import numpy as np
import astropy.io.fits as pyfits

from soxs.constants import sigma_to_fwhm
from soxs.utils import parse_prng, get_data_file, \
    image_pos, find_nearest


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


class GaussianPSF(PSF):
    _psf_type = "gaussian"

    def __init__(self, inst, prng=None):
        super().__init__(prng)
        fwhm = inst["psf"][1]
        plate_scale_arcsec = inst['fov']/inst['num_pixels']*60.0
        self.sigma = fwhm / sigma_to_fwhm / plate_scale_arcsec

    def scatter(self, x, y, e):
        n_evt = x.size
        x += self.prng.normal(loc=0.0, scale=self.sigma, size=n_evt)
        y += self.prng.normal(loc=0.0, scale=self.sigma, size=n_evt)
        return x, y


class ImagePSF(PSF):
    _psf_type = "image"

    def __init__(self, inst, prng=None):
        super().__init__(prng)
        img_file = get_data_file(inst['psf'][1])
        hdu = inst['psf'][2]
        plate_scale_rad = np.deg2rad(inst['fov']/inst['num_pixels']/60.0)
        plate_scale_mm = inst['focal_length']*1000.0*plate_scale_rad
        self.imhdu = pyfits.open(get_data_file(img_file))[hdu]
        self.imctr = np.array([self.imhdu.header["CRPIX1"],
                               self.imhdu.header["CRPIX2"]])
        self.scale = np.array([self.imhdu.header["CDELT1"],
                               self.imhdu.header["CDELT2"]])/plate_scale_mm

    def scatter(self, x, y, e):
        n_evt = x.size
        # This returns image coordinates from the PSF
        # image
        dx, dy = image_pos(self.imhdu.data, n_evt, self.prng)
        dx -= self.imctr[0]
        dy -= self.imctr[1]
        dx *= self.scale[0]
        dy *= self.scale[1]
        return x+dx, y+dy


class MultiImagePSF(PSF):
    _psf_type = "multi_image"

    def __init__(self, inst, prng=None):
        super().__init__(prng)
        self.img_file = get_data_file(inst['psf'][1])
        self.det_ctr = np.array(inst['aimpt_coords'])
        plate_scale_arcmin = inst['fov']/inst['num_pixels']
        plate_scale_deg = plate_scale_arcmin/60.0
        plate_scale_mm = inst['focal_length']*1e3*np.deg2rad(plate_scale_deg)
        img_e = []
        img_r = []
        img_i = []
        img_c = []
        img_s = []
        with pyfits.open(self.img_file, lazy_load_hdus=True) as f:
            for i, hdu in enumerate(f):
                if not isinstance(hdu, pyfits.ImageHDU):
                    continue
                img_e.append(hdu.header["ENERGY"])
                img_r.append(hdu.header["OFFAXIS"])
                img_c.append([hdu.header["CRPIX1"],
                              hdu.header["CRPIX2"]])
                img_s.append([hdu.header["CDELT1"],
                              hdu.header["CDELT2"]])
                img_i.append(i)
        self.img_e, ie = np.unique(img_e, return_inverse=True)
        self.img_r2, ir = np.unique(img_r, return_inverse=True)
        self.img_i = {j: (i, ie[j], ir[j]) for j, i in enumerate(img_i)}
        self.num_images = self.img_e.size
        self.img_r2 = (self.img_r2/plate_scale_arcmin)**2
        self.img_c = np.array(img_c)
        self.img_s = np.array(img_s)/plate_scale_mm

    def scatter(self, x, y, e):
        r2 = (x-self.det_ctr[0])**2 + (y-self.det_ctr[1])**2
        idx_e = find_nearest(self.img_e, e)
        idx_r = find_nearest(self.img_r2, r2)
        with pyfits.open(self.img_file) as f:
            for j in range(self.num_images):
                i, ie, ir = self.img_i[j]
                # This returns image coordinates from the PSF
                # image
                idxs = np.where((idx_e == ie) & (idx_r == ir))[0]
                dx, dy = image_pos(f[i].data, idxs.size, self.prng)
                dx -= self.img_c[j][0]
                dy -= self.img_c[j][1]
                dx *= self.img_s[j][0]
                dy *= self.img_s[j][1]
                x[idxs] += dx
                y[idxs] += dy
        return x, y
