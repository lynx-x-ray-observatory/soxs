from numbers import Number

import numpy as np
from astropy import wcs
from astropy.io import fits
from pathlib import Path

from soxs.events.utils import _combine_events, _region_filter, wcs_from_header
from soxs.utils import parse_value


def filter_events(
    evtfile,
    newfile,
    region=None,
    emin=None,
    emax=None,
    tmin=None,
    tmax=None,
    format="ds9",
    exclude=False,
    overwrite=False,
):
    r"""

    Parameters
    ----------
    evtfile : string
        The input events file to be read in.
    newfile : string
        The new event file that will be written.
    region : string, Region, or Regions, optional
        The region(s) to be used for the filtering. Default: None
    emin : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
        The minimum energy of the events to be included, in keV.
        Default is the lowest energy available.
    emax : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
        The maximum energy of the events to be included, in keV.
        Default is the highest energy available.
    tmin : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
        The minimum energy of the events to be included, in seconds.
        Default is the earliest time available.
    tmax : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
        The maximum energy of the events to be included, in seconds.
        Default is the latest time available.
    format : string, optional
        The file format specifier for the region, if supplied. "ds9",
        "crtf", "fits", etc. Default: "ds9"
    exclude : boolean, optional
        If True, the events in a specified *region* will be excluded instead of
        included. Default: False
    overwrite : boolean, optional
        Whether to overwrite an existing file with the
        same name. Default: False
    """
    with fits.open(evtfile) as f:
        hdu = f["EVENTS"]
        evt_mask = np.ones(hdu.data["ENERGY"].size, dtype="bool")
        if region is not None:
            evt_mask &= _region_filter(hdu, region, format=format, exclude=exclude)
        if emin is not None:
            emin = parse_value(emin, "keV") * 1000.0
            evt_mask &= hdu.data["ENERGY"] > emin
        if emax is not None:
            emax = parse_value(emax, "keV") * 1000.0
            evt_mask &= hdu.data["ENERGY"] < emax
        if tmin is not None:
            tmin = parse_value(tmin, "s")
            evt_mask &= hdu.data["TIME"] > tmin
        else:
            tmin = 0.0
        if tmax is not None:
            tmax = parse_value(tmax, "s")
            evt_mask &= hdu.data["TIME"] < tmax
        else:
            tmax = hdu.header["EXPOSURE"]
        hdu.data = hdu.data[evt_mask]
        hdu.header["EXPOSURE"] = tmax - tmin
        gtihdu = f["STDGTI"]
        gtihdu.data["START"][0] = tmin
        gtihdu.data["STOP"][0] = tmax
        gtihdu.header["TSTART"] = tmin
        gtihdu.header["TSTOP"] = tmax
        f.writeto(newfile, overwrite=overwrite)


def write_radial_profile(
    evt_file,
    out_file,
    ctr,
    rmin,
    rmax,
    nbins,
    ctr_type="celestial",
    emin=None,
    emax=None,
    expmap_file=None,
    overwrite=False,
):
    r"""
    Bin up events into a radial profile and write them to a FITS
    table.

    Parameters
    ----------
    evt_file : string
        Input event file.
    out_file : string
        The output file to write the profile to.
    ctr : array-like
        The central coordinate of the profile. Can either be in
        celestial coordinates (the default) or "physical" pixel
        coordinates. If the former, the ``ctr_type`` keyword
        argument must be explicity set to "physical".
    rmin : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The minimum radius of the profile, in arcseconds.
    rmax : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The maximum radius of the profile, in arcseconds.
    nbins : integer
        The number of bins in the profile.
    ctr_type : string, optional
        The type of center coordinate. Either "celestial" for
        (RA, Dec) coordinates (the default), or "physical" for
        pixel coordinates.
    emin : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
        The minimum energy of the events to be binned in keV.
        Default is the lowest energy available.
    emax : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
        The maximum energy of the events to be binned in keV.
        Default is the highest energy available.
    overwrite : boolean, optional
        Whether to overwrite an existing file with the
        same name. Default: False
    expmap_file : string, optional
        Supply an exposure map file to determine fluxes.
        Default: None
    """
    rmin = parse_value(rmin, "arcsec")
    rmax = parse_value(rmax, "arcsec")
    f = fits.open(evt_file)
    hdu = f["EVENTS"]
    orig_dx = hdu.header["TCDLT3"]
    e = hdu.data["ENERGY"]
    if emin is None:
        emin = e.min()
    else:
        emin = parse_value(emin, "keV")
        emin *= 1000.0
    if emax is None:
        emax = e.max()
    else:
        emax = parse_value(emax, "keV")
        emax *= 1000.0
    idxs = np.logical_and(e > emin, e < emax)
    x = hdu.data["X"][idxs]
    y = hdu.data["Y"][idxs]
    exp_time = hdu.header["EXPOSURE"]
    w = wcs_from_header(hdu.header)
    dtheta = np.abs(w.wcs.cdelt[1]) * 3600.0
    f.close()

    if ctr_type == "celestial":
        ctr = w.all_world2pix(ctr[0], ctr[1], 1)

    r = np.sqrt((x - ctr[0]) ** 2 + (y - ctr[1]) ** 2)
    rr = np.linspace(rmin / dtheta, rmax / dtheta, nbins + 1)
    C = np.histogram(r, bins=rr)[0]
    rbin = rr * dtheta
    rmid = 0.5 * (rbin[1:] + rbin[:-1])

    A = np.pi * (rbin[1:] ** 2 - rbin[:-1] ** 2)

    Cerr = np.sqrt(C)

    R = C / exp_time
    Rerr = Cerr / exp_time

    S = R / A
    Serr = Rerr / A

    col1 = fits.Column(name="RLO", format="D", unit="arcsec", array=rbin[:-1])
    col2 = fits.Column(name="RHI", format="D", unit="arcsec", array=rbin[1:])
    col3 = fits.Column(name="RMID", format="D", unit="arcsec", array=rmid)
    col4 = fits.Column(name="AREA", format="D", unit="arcsec**2", array=A)
    col5 = fits.Column(name="NET_COUNTS", format="D", unit="count", array=C)
    col6 = fits.Column(name="NET_ERR", format="D", unit="count", array=Cerr)
    col7 = fits.Column(name="NET_RATE", format="D", unit="count/s", array=R)
    col8 = fits.Column(name="ERR_RATE", format="D", unit="count/s", array=Rerr)
    col9 = fits.Column(name="SUR_BRI", format="D", unit="count/s/arcsec**2", array=S)
    col10 = fits.Column(
        name="SUR_BRI_ERR", format="1D", unit="count/s/arcsec**2", array=Serr
    )

    coldefs = [col1, col2, col3, col4, col5, col6, col7, col8, col9, col10]

    if expmap_file is not None:
        f = fits.open(expmap_file)
        ehdu = f["EXPMAP"]
        wexp = wcs.WCS(header=ehdu.header)
        cel = w.all_pix2world(ctr[0], ctr[1], 1)
        ectr = wexp.all_world2pix(cel[0], cel[1], 1)
        exp = ehdu.data[:, :]
        nx, ny = exp.shape
        reblock = ehdu.header["CDELT2"] / orig_dx
        x, y = np.mgrid[1 : nx + 1, 1 : ny + 1]
        r = np.sqrt((x - ectr[0]) ** 2 + (y - ectr[1]) ** 2)
        f.close()
        E = (
            np.histogram(r, bins=rr / reblock, weights=exp)[0]
            / np.histogram(r, bins=rr / reblock)[0]
        )
        with np.errstate(invalid="ignore", divide="ignore"):
            F = R / E
            Ferr = Rerr / E
        SF = F / A
        SFerr = Ferr / A
        col11 = fits.Column(name="MEAN_SRC_EXP", format="D", unit="cm**2", array=E)
        col12 = fits.Column(name="NET_FLUX", format="D", unit="count/s/cm**2", array=F)
        col13 = fits.Column(
            name="NET_FLUX_ERR", format="D", unit="count/s/cm**2", array=Ferr
        )
        col14 = fits.Column(
            name="SUR_FLUX", format="D", unit="count/s/cm**2/arcsec**2", array=SF
        )
        col15 = fits.Column(
            name="SUR_FLUX_ERR", format="D", unit="count/s/cm**2/arcsec**2", array=SFerr
        )
        coldefs += [col11, col12, col13, col14, col15]

    tbhdu = fits.BinTableHDU.from_columns(fits.ColDefs(coldefs))
    tbhdu.name = "PROFILE"

    hdulist = fits.HDUList([fits.PrimaryHDU(), tbhdu])

    hdulist.writeto(out_file, overwrite=overwrite)


def fill_regions(
    in_img, out_img, src_reg, bkg_val, median=False, overwrite=False, format="ds9"
):
    """
    Fill regions which have been removed from an image (e.g. by
    wavdetect) with a Poisson distribution of counts, either
    from a single background region, a number of background regions
    equal to the number of source regions, or from a single background
    value.

    Parameters
    ----------
    in_img : string
        The file path to the input image.
    out_img : string
        The file path to the image to be written with the filled
        regions.
    src_reg : string, Region, or Regions
        The region(s) which will be filled.
    bkg_val : float, string, Region, or Regions
        The background value(s) to use for the Poisson distribution.
        Can be a single value, a region, or a list of regions.
    median : boolean
        If True, then the median value of the counts within the
        region will be used as the mean of the Poisson distribution
    overwrite : boolean, optional
        Whether to overwrite an existing file with
        the same name. Default: False
    format : string, optional
        The file format specifier for the region. "ds9",
        "crtf", "fits", etc. Default: "ds9"
    """
    from regions import PixelRegion, Region, Regions, SkyRegion

    fill_op = np.median if median else np.mean
    if isinstance(src_reg, str):
        if Path(src_reg).exists():
            src_reg = Regions.read(src_reg, format=format)
        else:
            src_reg = Regions.parse(src_reg, format=format)
    elif not isinstance(src_reg, (Region, Regions)):
        raise RuntimeError("'src_reg' argument is not valid!")
    if isinstance(src_reg, Region):
        src_reg = [src_reg]
    if not isinstance(bkg_val, Number):
        if isinstance(bkg_val, str):
            if Path(bkg_val).exists():
                bkg_val = Regions.read(bkg_val, format=format)
            else:
                bkg_val = Regions.parse(bkg_val, format=format)
        elif not isinstance(bkg_val, (Region, Regions)):
            raise RuntimeError("'bkg_val' argument is not valid!")
    if isinstance(bkg_val, (Number, Region)):
        bkg_val = [bkg_val] * len(src_reg)
    if len(src_reg) != len(bkg_val):
        raise ValueError(
            "The number of background regions must either "
            "be the same as the number of the source regions, "
            "or there must be one background region!"
        )
    with fits.open(in_img) as f:
        if f[0].is_image and f[0].header["NAXIS"] == 2:
            hdu = f[0]
        else:
            hdu = f[1]
        for src_r, bkg_v in zip(src_reg, bkg_val):
            if isinstance(src_r, PixelRegion):
                src_mask = src_r.to_mask(mode="exact").astype("bool")
            elif isinstance(src_r, SkyRegion):
                w = wcs.WCS(header=hdu.header)
                src_mask = src_r.to_pixel(w).to_mask(mode="exact")
            else:
                raise NotImplementedError
            src_mask = src_mask.to_image(hdu.shape).astype("bool")
            if isinstance(bkg_v, Number):
                lam = bkg_v
            else:
                if isinstance(bkg_v, PixelRegion):
                    bkg_mask = bkg_v.to_mask().astype("bool")
                elif isinstance(bkg_v, SkyRegion):
                    w = wcs.WCS(header=hdu.header)
                    bkg_mask = bkg_v.to_pixel(w).to_mask()
                else:
                    raise NotImplementedError
                bkg_mask = bkg_mask.to_image(hdu.shape).astype("bool")
                lam = fill_op(hdu.data[bkg_mask])
            n_src = src_mask.sum()
            hdu.data[src_mask] = np.random.poisson(lam=lam, size=n_src)
        f.writeto(out_img, overwrite=overwrite)


def merge_event_files(input_files, output_file, overwrite=False):
    """
    Merge SOXS-generated event files together. The WCS used for the final
    file will be based on the first source file, so the original celestial
    coordinates for the other files will be discarded.

    Parameters
    ----------
    input_files : list of strings
        The events file(s) containing the background events.
    output_file : string
        The merged events file.
    overwrite : boolean, optional
        Whether to overwrite an existing file with
        the same name. Default: False
    """
    with fits.open(input_files[0], memmap=True) as f:
        wcs_out = wcs_from_header(f["EVENTS"].header)
        shape_out = 2.0 * wcs_out.wcs.crpix - 1
    _combine_events(input_files, wcs_out, shape_out, output_file, overwrite=overwrite)
