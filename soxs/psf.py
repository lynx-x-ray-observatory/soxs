import numpy as np
from astropy.io import fits
from astropy.units import Quantity

from soxs.constants import sigma_to_fwhm
from soxs.lib.psf_cdf import eef_cdf, score_psf
from soxs.utils import convert_endian, get_data_file, image_pos, parse_prng

psf_model_registry = {}


class RegisteredPSFModel(type):
    def __init__(cls, name, b, d):
        type.__init__(cls, name, b, d)
        if hasattr(cls, "_psf_type"):
            psf_model_registry[cls._psf_type] = cls


class PSF(metaclass=RegisteredPSFModel):
    _psf_type = "null"

    def __init__(self, prng=None):
        self.prng = parse_prng(prng)

    def __str__(self):
        return self._psf_type


class GaussianPSF(PSF):
    _psf_type = "gaussian"

    def __init__(self, inst, prng=None):
        super().__init__(prng)
        fwhm = inst["psf"][1]
        plate_scale_arcsec = inst["fov"] / inst["num_pixels"] * 60.0
        self.sigma = fwhm / sigma_to_fwhm / plate_scale_arcsec

    def scatter(self, x, y, e):
        n_evt = x.size
        x += self.prng.normal(loc=0.0, scale=self.sigma, size=n_evt)
        y += self.prng.normal(loc=0.0, scale=self.sigma, size=n_evt)
        return x, y


class EEFPSF(PSF):
    _psf_type = "eef"

    def __init__(self, inst, prng=None):
        super().__init__(prng)
        eef_file = get_data_file(inst["psf"][1])
        hdu = inst["psf"][2]
        plate_scale_arcsec = inst["fov"] / inst["num_pixels"] * 60.0
        self.eefhdu = fits.open(get_data_file(eef_file))[hdu]
        unit = getattr(self.eefhdu.columns["psfrad"], "unit", "arcsec")
        self.scale = Quantity(1.0, unit).to_value("arcsec") / plate_scale_arcsec

    def scatter(self, x, y, e):
        n_evt = x.size
        # This returns image coordinates from the PSF
        # image
        rad = self.eefhdu.data["psfrad"] * self.scale
        cdf = self.eefhdu.data["eef"]
        randvec = self.prng.uniform(low=cdf[0], high=cdf[-1], size=n_evt)
        r = np.interp(randvec, cdf, rad)
        phi = self.prng.uniform(low=0.0, high=2.0 * np.pi, size=n_evt)
        return x + r * np.cos(phi), y + r * np.sin(phi)


class MultiEEFPSF(PSF):
    _psf_type = "multi_eef"

    def __init__(self, inst, prng=None):
        super().__init__(prng)
        self.eef_file = get_data_file(inst["psf"][1])
        self.eef_type = inst["psf"][2]
        self.det_ctr = np.array(inst["aimpt_coords"])
        plate_scale_arcmin = inst["fov"] / inst["num_pixels"]
        plate_scale_arcsec = plate_scale_arcmin * 60.0
        with fits.open(self.eef_file, lazy_load_hdus=True) as f:
            if self.eef_type == 1:
                eef_e = []
                eef_r = []
                eef_i = []
                eef_u = []
                for i, hdu in enumerate(f):
                    if hdu.header["XTENSION"] != "BINTABLE":
                        continue
                    eef_e.append(hdu.header["ENERGY"])
                    key = "THETA" if "OFFAXIS" not in hdu.header else "OFFAXIS"
                    eef_r.append(hdu.header[key])
                    eef_i.append(i)
                    eef_u.append(getattr(hdu.columns["psfrad"], "unit", "arcsec"))
                eef_e = np.array(eef_e)
                eef_r = np.array(eef_r)
                units = list(set(eef_u))
                if len(units) > 1:
                    raise RuntimeError("More than one psfrad unit detected!!")
                unit = units[0]
            elif self.eef_type == 2:
                hdu = f[1]
                unit = getattr(hdu.columns["Radius"], "unit", "arcsec")
                eef_e = convert_endian(hdu.data["Energy"])
                eef_r = Quantity(
                    convert_endian(hdu.data["Theta"]), hdu.columns["Theta"].unit
                ).to_value("arcmin")
                eef_i = np.arange(eef_e.size)
        eef_r = (eef_r / plate_scale_arcmin) ** 2
        if np.all(eef_e > 100.0):
            # this is probably in eV
            eef_e *= 1.0e-3
        self.eef_bins = np.array([eef_e, eef_r]).astype("float64")
        self.eef_i = eef_i
        self.num_eef = len(eef_e)
        self.eef_s = Quantity(1.0, unit).to_value("arcsec") / plate_scale_arcsec

    def scatter(self, x, y, e):
        r2 = (x - self.det_ctr[0]) ** 2 + (y - self.det_ctr[1]) ** 2
        idx_score = score_psf(self.eef_bins[0], self.eef_bins[1], e, r2)
        if self.eef_type == 1:
            x, y = self._scatter1(x, y, idx_score)
        elif self.eef_type == 2:
            x, y = self._scatter2(x, y, idx_score)
        return x, y

    def _scatter1(self, x, y, idx_score):
        n_in = x.size
        n_out = 0
        with fits.open(self.eef_file) as f:
            for j in range(self.num_eef):
                # This returns image coordinates from the EEF
                # curve
                idxs = np.where(idx_score == j)[0]
                n_out += idxs.size
                rad = f[self.eef_i[j]].data["psfrad"] * self.eef_s
                cdf = f[self.eef_i[j]].data["eef"]
                randvec = self.prng.uniform(low=cdf[0], high=cdf[-1], size=idxs.size)
                r = np.interp(randvec, cdf, rad)
                phi = self.prng.uniform(low=0.0, high=2.0 * np.pi, size=idxs.size)
                x[idxs] += r * np.cos(phi)
                y[idxs] += r * np.sin(phi)
        if n_in != n_out:
            raise ValueError(
                f"The number of photons scattered by the PSF "
                f"({n_out}) does not equal the input number "
                f"({n_in})!"
            )
        return x, y

    def _scatter2(self, x, y, idx_score):
        n_evt = x.size
        with fits.open(self.eef_file) as f:
            hdu = f[1]
            rad = convert_endian(hdu.data["Radius"]).astype("float64") * self.eef_s
            cdf = np.insert(
                convert_endian(hdu.data["EEF"]).astype("float64"), 0, 0.0, axis=1
            )
            cdf /= cdf[:, -1, np.newaxis]
        randvec = self.prng.uniform(low=0.0, high=1.0, size=n_evt)
        r = eef_cdf(idx_score, randvec, rad, cdf)
        phi = self.prng.uniform(low=0.0, high=2.0 * np.pi, size=n_evt)
        x += r * np.cos(phi)
        y += r * np.sin(phi)
        return x, y


class ImagePSF(PSF):
    _psf_type = "image"

    def __init__(self, inst, prng=None):
        super().__init__(prng)
        img_file = get_data_file(inst["psf"][1])
        hdu = inst["psf"][2]
        plate_scale_arcmin = inst["fov"] / inst["num_pixels"]
        plate_scale_deg = plate_scale_arcmin / 60.0
        plate_scale_mm = inst["focal_length"] * 1e3 * np.deg2rad(plate_scale_deg)
        self.imhdu = fits.open(get_data_file(img_file))[hdu]
        self.imctr = np.array(
            [self.imhdu.header["CRPIX1"], self.imhdu.header["CRPIX2"]]
        )
        unit = self.imhdu.header.get("CUNIT1", "mm")
        self.scale = Quantity(
            [self.imhdu.header["CDELT1"], self.imhdu.header["CDELT2"]], unit
        ).to_value("mm")
        self.scale /= plate_scale_mm

    def scatter(self, x, y, e):
        n_evt = x.size
        # This returns image coordinates from the PSF
        # image
        dx, dy = image_pos(self.imhdu.data, n_evt, self.prng)
        dx -= self.imctr[0]
        dy -= self.imctr[1]
        dx *= self.scale[0]
        dy *= self.scale[1]
        return x + dx, y + dy


class MultiImagePSF(PSF):
    _psf_type = "multi_image"

    def __init__(self, inst, prng=None):
        super().__init__(prng)
        self.img_file = get_data_file(inst["psf"][1])
        self.det_ctr = np.array(inst["aimpt_coords"])
        plate_scale_arcmin = inst["fov"] / inst["num_pixels"]
        plate_scale_deg = plate_scale_arcmin / 60.0
        plate_scale_mm = inst["focal_length"] * 1e3 * np.deg2rad(plate_scale_deg)
        img_e = []
        img_r = []
        img_i = []
        img_c = []
        img_s = []
        img_u = []
        with fits.open(self.img_file, lazy_load_hdus=True) as f:
            for i, hdu in enumerate(f):
                if not hdu.is_image or hdu.header["NAXIS"] != 2:
                    continue
                img_e.append(hdu.header["ENERGY"])
                key = "THETA" if "OFFAXIS" not in hdu.header else "OFFAXIS"
                img_r.append(hdu.header[key])
                img_c.append([hdu.header["CRPIX1"], hdu.header["CRPIX2"]])
                img_s.append([hdu.header["CDELT1"], hdu.header["CDELT2"]])
                img_i.append(i)
                img_u.append(hdu.header.get("CUNIT1", "mm"))
        img_e = np.array(img_e)
        img_r = (np.array(img_r) / plate_scale_arcmin) ** 2
        if np.all(img_e > 100.0):
            # this is probably in eV
            img_e *= 1.0e-3
        self.img_bins = np.array([img_e, img_r]).astype("float64")
        self.img_i = img_i
        self.num_images = len(img_e)
        self.img_c = np.array(img_c)
        unit = list(set(img_u))
        if len(unit) > 1:
            raise RuntimeError("More than one delta unit detected!!")
        self.img_s = Quantity(img_s, unit[0]).to_value("mm") / plate_scale_mm

    def scatter(self, x, y, e):
        r2 = (x - self.det_ctr[0]) ** 2 + (y - self.det_ctr[1]) ** 2
        idx_score = score_psf(self.img_bins[0], self.img_bins[1], e, r2)
        n_in = x.size
        n_out = 0
        with fits.open(self.img_file) as f:
            for j in range(self.num_images):
                # This returns image coordinates from the PSF
                # image
                idxs = np.where(idx_score == j)[0]
                n_out += idxs.size
                dx, dy = image_pos(f[self.img_i[j]].data, idxs.size, self.prng)
                dx -= self.img_c[j][0]
                dy -= self.img_c[j][1]
                dx *= self.img_s[j][0]
                dy *= self.img_s[j][1]
                x[idxs] += dx
                y[idxs] += dy
        if n_in != n_out:
            raise ValueError(
                f"The number of photons scattered by the PSF "
                f"({n_out}) does not equal the input number "
                f"({n_in})!"
            )
        return x, y
