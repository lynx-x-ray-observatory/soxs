import numpy as np
from astropy import wcs
from astropy.io import fits
from pathlib import Path, PurePath
from tqdm.auto import tqdm

from soxs.instrument_registry import instrument_registry
from soxs.utils import create_region, get_rot_mat, mylog, parse_value


def wcs_from_header(h):
    w = wcs.WCS(naxis=2)
    w.wcs.crval = [h["TCRVL2"], h["TCRVL3"]]
    w.wcs.crpix = [h["TCRPX2"], h["TCRPX3"]]
    w.wcs.cdelt = [h["TCDLT2"], h["TCDLT3"]]
    w.wcs.ctype = [h["TCTYP2"], h["TCTYP3"]]
    w.wcs.cunit = [h["TCUNI2"], h["TCUNI3"]]
    return w


def write_event_file(events, parameters, filename, overwrite=False):
    from astropy.time import Time, TimeDelta

    mylog.info("Writing events to file %s.", filename)

    t_begin = Time.now()
    dt = TimeDelta(parameters["exposure_time"], format="sec")
    t_end = t_begin + dt

    col_x = fits.Column(name="X", format="D", unit="pixel", array=events["xpix"])
    col_y = fits.Column(name="Y", format="D", unit="pixel", array=events["ypix"])
    col_e = fits.Column(
        name="ENERGY", format="E", unit="eV", array=events["energy"] * 1000.0
    )
    col_dx = fits.Column(name="DETX", format="D", unit="pixel", array=events["detx"])
    col_dy = fits.Column(name="DETY", format="D", unit="pixel", array=events["dety"])
    col_id = fits.Column(
        name="CCD_ID", format="J", unit="pixel", array=events["ccd_id"]
    )

    chantype = parameters["channel_type"].upper()
    if chantype == "PHA":
        cunit = "adu"
    elif chantype == "PI":
        cunit = "Chan"
    col_ch = fits.Column(name=chantype, format="1J", unit=cunit, array=events[chantype])

    col_t = fits.Column(name="TIME", format="1D", unit="s", array=events["time"])

    cols = [col_e, col_x, col_y, col_ch, col_t, col_dx, col_dy, col_id]

    coldefs = fits.ColDefs(cols)
    tbhdu = fits.BinTableHDU.from_columns(coldefs)
    tbhdu.name = "EVENTS"

    tbhdu.header["MTYPE1"] = "sky"
    tbhdu.header["MFORM1"] = "x,y"
    tbhdu.header["MTYPE2"] = "EQPOS"
    tbhdu.header["MFORM2"] = "RA,DEC"
    tbhdu.header["TCTYP2"] = "RA---TAN"
    tbhdu.header["TCTYP3"] = "DEC--TAN"
    tbhdu.header["TCRVL2"] = parameters["sky_center"][0]
    tbhdu.header["TCRVL3"] = parameters["sky_center"][1]
    tbhdu.header["TCDLT2"] = -parameters["plate_scale"]
    tbhdu.header["TCDLT3"] = parameters["plate_scale"]
    tbhdu.header["TCRPX2"] = parameters["pix_center"][0]
    tbhdu.header["TCRPX3"] = parameters["pix_center"][1]
    tbhdu.header["TCUNI2"] = "deg"
    tbhdu.header["TCUNI3"] = "deg"
    tbhdu.header["TLMIN2"] = 0.5
    tbhdu.header["TLMIN3"] = 0.5
    tbhdu.header["TLMAX2"] = 2.0 * parameters["num_pixels"] + 0.5
    tbhdu.header["TLMAX3"] = 2.0 * parameters["num_pixels"] + 0.5
    tbhdu.header["TLMIN4"] = parameters["chan_lim"][0]
    tbhdu.header["TLMAX4"] = parameters["chan_lim"][1]
    tbhdu.header["TLMIN6"] = -0.5 * parameters["num_pixels"]
    tbhdu.header["TLMAX6"] = 0.5 * parameters["num_pixels"]
    tbhdu.header["TLMIN7"] = -0.5 * parameters["num_pixels"]
    tbhdu.header["TLMAX7"] = 0.5 * parameters["num_pixels"]
    tbhdu.header["EXPOSURE"] = parameters["exposure_time"]
    tbhdu.header["TSTART"] = 0.0
    tbhdu.header["TSTOP"] = parameters["exposure_time"]
    tbhdu.header["HDUVERS"] = "1.1.0"
    tbhdu.header["RADECSYS"] = "FK5"
    tbhdu.header["EQUINOX"] = 2000.0
    tbhdu.header["HDUCLASS"] = "OGIP"
    tbhdu.header["HDUCLAS1"] = "EVENTS"
    tbhdu.header["HDUCLAS2"] = "ACCEPTED"
    tbhdu.header["DATE"] = t_begin.tt.isot
    tbhdu.header["DATE-OBS"] = t_begin.tt.isot
    tbhdu.header["DATE-END"] = t_end.tt.isot
    tbhdu.header["RESPFILE"] = PurePath(parameters["rmf"]).parts[-1]
    tbhdu.header["PHA_BINS"] = parameters["nchan"]
    tbhdu.header["ANCRFILE"] = PurePath(parameters["arf"]).parts[-1]
    tbhdu.header["CHANTYPE"] = parameters["channel_type"]
    tbhdu.header["MISSION"] = parameters["mission"]
    tbhdu.header["TELESCOP"] = parameters["telescope"]
    tbhdu.header["INSTRUME"] = parameters["instrument"]
    tbhdu.header["RA_PNT"] = parameters["sky_center"][0]
    tbhdu.header["DEC_PNT"] = parameters["sky_center"][1]
    tbhdu.header["ROLL_PNT"] = parameters["roll_angle"]
    tbhdu.header["AIMPT_X"] = parameters["aimpt_coords"][0]
    tbhdu.header["AIMPT_Y"] = parameters["aimpt_coords"][1]
    tbhdu.header["AIMPT_DX"] = parameters["aimpt_shift"][0]
    tbhdu.header["AIMPT_DY"] = parameters["aimpt_shift"][1]
    if parameters["dither_params"]["dither_on"]:
        tbhdu.header["DITHXAMP"] = parameters["dither_params"]["x_amp"]
        tbhdu.header["DITHYAMP"] = parameters["dither_params"]["y_amp"]
        tbhdu.header["DITHXPER"] = parameters["dither_params"]["x_period"]
        tbhdu.header["DITHYPER"] = parameters["dither_params"]["y_period"]

    start = fits.Column(name="START", format="1D", unit="s", array=np.array([0.0]))
    stop = fits.Column(
        name="STOP",
        format="1D",
        unit="s",
        array=np.array([parameters["exposure_time"]]),
    )

    tbhdu_gti = fits.BinTableHDU.from_columns([start, stop])
    tbhdu_gti.name = "STDGTI"
    tbhdu_gti.header["TSTART"] = 0.0
    tbhdu_gti.header["TSTOP"] = parameters["exposure_time"]
    tbhdu_gti.header["HDUCLASS"] = "OGIP"
    tbhdu_gti.header["HDUCLAS1"] = "GTI"
    tbhdu_gti.header["HDUCLAS2"] = "STANDARD"
    tbhdu_gti.header["RADECSYS"] = "FK5"
    tbhdu_gti.header["EQUINOX"] = 2000.0
    tbhdu_gti.header["DATE"] = t_begin.tt.isot
    tbhdu_gti.header["DATE-OBS"] = t_begin.tt.isot
    tbhdu_gti.header["DATE-END"] = t_end.tt.isot

    hdulist = [fits.PrimaryHDU(), tbhdu, tbhdu_gti]

    fits.HDUList(hdulist).writeto(filename, overwrite=overwrite)


def make_exposure_map(
    event_file,
    expmap_file,
    energy,
    weights=None,
    asol_file=None,
    normalize=True,
    overwrite=False,
    reblock=1,
    nhistx=16,
    nhisty=16,
):
    """
    Make an exposure map for a SOXS event file, and optionally write
    an aspect solution file. The exposure map will be created by
    binning an aspect histogram over the range of the aspect solution.

    Parameters
    ----------
    event_file : string
        The path to the event file to use for making the exposure map.
    expmap_file : string
        The path to write the exposure map file to.
    energy : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, or NumPy array
        The energy in keV to use when computing the exposure map, or
        a set of energies to be used with the *weights* parameter. If
        providing a set, it must be in keV.
    weights : array-like, optional
        The weights to use with a set of energies given in the
        *energy* parameter. Used to create a more accurate exposure
        map weighted by a range of energies. Default: None
    asol_file : string, optional
        The path to write the aspect solution file to, if desired.
        Default: None
    normalize : boolean, optional
        If True, the exposure map will be divided by the exposure time
        so that the map's units are cm**2. Default: True
    overwrite : boolean, optional
        Whether to overwrite an existing file. Default: False
    reblock : integer, optional
        Supply an integer power of 2 here to make an exposure map
        with a different binning. Default: 1
    nhistx : integer, optional
        The number of bins in the aspect histogram in the DETX
        direction. Default: 16
    nhisty : integer, optional
        The number of bins in the aspect histogram in the DETY
        direction. Default: 16
    order : integer, optional
        The interpolation order to use when making the exposure map.
        Default: 1
    """
    from scipy.ndimage import rotate

    from soxs.instrument import perform_dither
    from soxs.response import AuxiliaryResponseFile

    if isinstance(energy, np.ndarray) and weights is None:
        raise RuntimeError(
            "Must supply a single value for the energy if " "you do not supply weights!"
        )
    if not isinstance(energy, np.ndarray):
        energy = parse_value(energy, "keV")
    f_evt = fits.open(event_file)
    hdu = f_evt["EVENTS"]
    arf = AuxiliaryResponseFile(hdu.header["ANCRFILE"])
    exp_time = hdu.header["EXPOSURE"]
    nx = int(hdu.header["TLMAX2"] - 0.5) // 2
    ny = int(hdu.header["TLMAX3"] - 0.5) // 2
    ra0 = hdu.header["TCRVL2"]
    dec0 = hdu.header["TCRVL3"]
    xdel = hdu.header["TCDLT2"]
    ydel = hdu.header["TCDLT3"]
    x0 = hdu.header["TCRPX2"]
    y0 = hdu.header["TCRPX3"]
    xaim = hdu.header.get("AIMPT_X", 0.0)
    yaim = hdu.header.get("AIMPT_Y", 0.0)
    xaim += hdu.header.get("AIMPT_DX", 0.0)
    yaim += hdu.header.get("AIMPT_DY", 0.0)
    roll = hdu.header["ROLL_PNT"]
    instr = instrument_registry[hdu.header["INSTRUME"].lower()]
    dither_params = {}
    if "DITHXAMP" in hdu.header:
        dither_params["x_amp"] = hdu.header["DITHXAMP"]
        dither_params["y_amp"] = hdu.header["DITHYAMP"]
        dither_params["x_period"] = hdu.header["DITHXPER"]
        dither_params["y_period"] = hdu.header["DITHYPER"]
        dither_params["plate_scale"] = ydel * 3600.0
        dither_params["dither_on"] = True
    else:
        dither_params["dither_on"] = False
    f_evt.close()

    # Create time array for aspect solution
    dt = 1.0  # Seconds
    t = np.arange(0.0, exp_time + dt, dt)

    # Construct WCS
    w = wcs.WCS(naxis=2)
    w.wcs.crval = [ra0, dec0]
    w.wcs.crpix = [x0, y0]
    w.wcs.cdelt = [xdel, ydel]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.cunit = ["deg"] * 2

    # Create aspect solution if we had dithering.
    # otherwise just set the offsets to zero
    if dither_params["dither_on"]:
        x_off, y_off = perform_dither(t, dither_params)
        # Make the aspect histogram
        x_amp = dither_params["x_amp"] / dither_params["plate_scale"]
        y_amp = dither_params["y_amp"] / dither_params["plate_scale"]
        x_edges = np.linspace(-x_amp, x_amp, nhistx + 1, endpoint=True)
        y_edges = np.linspace(-y_amp, y_amp, nhisty + 1, endpoint=True)
        asphist = np.histogram2d(x_off, y_off, (x_edges, y_edges))[0]
        asphist *= dt
        x_mid = 0.5 * (x_edges[1:] + x_edges[:-1]) / reblock
        y_mid = 0.5 * (y_edges[1:] + y_edges[:-1]) / reblock
    else:
        asphist = exp_time * np.ones((1, 1))

    # Determine the effective area
    eff_area = arf.interpolate_area(energy).value
    if weights is not None:
        eff_area = np.average(eff_area, weights=weights)

    rtypes = []
    args = []
    for chip in instr["chips"]:
        rtypes.append(chip[0])
        args.append(np.array(chip[1:]) / reblock)

    xdet0 = 0.5 * (2 * nx // reblock + 1)
    ydet0 = 0.5 * (2 * ny // reblock + 1)
    xaim //= reblock
    yaim //= reblock
    dx = xdet0 - xaim - 1.0
    dy = ydet0 - yaim - 1.0

    if dither_params["dither_on"]:
        niterx = nhistx
        nitery = nhisty
    else:
        niterx = 1
        nitery = 1

    expmap = np.zeros((2 * nx // reblock, 2 * ny // reblock))
    niter = niterx * nitery
    pbar = tqdm(leave=True, total=niter, desc="Creating exposure map ")
    for i in range(niterx):
        for j in range(nitery):
            chips, _ = create_region(rtypes[0], args[0], dx + x_mid[i], dy + y_mid[j])
            for rtype, arg in zip(rtypes[1:], args[1:]):
                r, _ = create_region(rtype, arg, dx + x_mid[i], dy + y_mid[j])
                chips = chips | r
            dexp = chips.to_mask().to_image(expmap.shape).astype("float64")
            expmap += dexp * asphist[i, j]
        pbar.update(nitery)
    pbar.close()

    expmap *= eff_area
    if normalize:
        expmap /= exp_time

    if roll != 0.0:
        rotate(expmap, roll, output=expmap, reshape=False)

    expmap[expmap < 0.0] = 0.0

    map_header = {
        "EXPOSURE": exp_time,
        "MTYPE1": "EQPOS",
        "MFORM1": "RA,DEC",
        "CTYPE1": "RA---TAN",
        "CTYPE2": "DEC--TAN",
        "CRVAL1": ra0,
        "CRVAL2": dec0,
        "CUNIT1": "deg",
        "CUNIT2": "deg",
        "CDELT1": xdel * reblock,
        "CDELT2": ydel * reblock,
        "CRPIX1": 0.5 * (2.0 * nx // reblock + 1),
        "CRPIX2": 0.5 * (2.0 * ny // reblock + 1),
    }

    map_hdu = fits.ImageHDU(expmap, header=fits.Header(map_header))
    map_hdu.name = "EXPMAP"
    map_hdu.writeto(expmap_file, overwrite=overwrite)

    if asol_file is not None:

        if dither_params["dither_on"]:

            det = np.array([x_off, y_off])

            pix = np.dot(get_rot_mat(roll).T, det)

            ra, dec = w.wcs_pix2world(pix[0, :] + x0, pix[1, :] + y0, 1)

            col_t = fits.Column(name="time", format="D", unit="s", array=t)
            col_ra = fits.Column(name="ra", format="D", unit="deg", array=ra)
            col_dec = fits.Column(name="dec", format="D", unit="deg", array=dec)

            coldefs = fits.ColDefs([col_t, col_ra, col_dec])
            tbhdu = fits.BinTableHDU.from_columns(coldefs)
            tbhdu.name = "ASPSOL"
            tbhdu.header["EXPOSURE"] = exp_time

            hdulist = [fits.PrimaryHDU(), tbhdu]

            fits.HDUList(hdulist).writeto(asol_file, overwrite=overwrite)

        else:

            mylog.warning(
                "Refusing to write an aspect solution file because "
                "there was no dithering."
            )


def _write_spectrum(
    bins, spec, exp_time, spectype, parameters, specfile, overwrite=False
):

    col1 = fits.Column(name="CHANNEL", format="1J", array=bins)
    col2 = fits.Column(name=spectype.upper(), format="1D", array=bins.astype("float64"))
    col3 = fits.Column(name="COUNTS", format="1J", array=spec.astype("int32"))
    col4 = fits.Column(name="COUNT_RATE", format="1D", array=spec / exp_time)

    coldefs = fits.ColDefs([col1, col2, col3, col4])

    tbhdu = fits.BinTableHDU.from_columns(coldefs)
    tbhdu.name = "SPECTRUM"

    tbhdu.header["DETCHANS"] = spec.size
    tbhdu.header["TOTCTS"] = spec.sum()
    tbhdu.header["EXPOSURE"] = exp_time
    tbhdu.header["LIVETIME"] = exp_time
    tbhdu.header["CONTENT"] = spectype
    tbhdu.header["HDUCLASS"] = "OGIP"
    tbhdu.header["HDUCLAS1"] = "SPECTRUM"
    tbhdu.header["HDUCLAS2"] = "TOTAL"
    tbhdu.header["HDUCLAS3"] = "TYPE:I"
    tbhdu.header["HDUCLAS4"] = "COUNT"
    tbhdu.header["HDUVERS"] = "1.1.0"
    tbhdu.header["HDUVERS1"] = "1.1.0"
    tbhdu.header["CHANTYPE"] = spectype
    tbhdu.header["BACKFILE"] = "none"
    tbhdu.header["CORRFILE"] = "none"
    tbhdu.header["POISSERR"] = True
    for key in ["RESPFILE", "ANCRFILE", "MISSION", "TELESCOP", "INSTRUME"]:
        tbhdu.header[key] = parameters[key]
    tbhdu.header["AREASCAL"] = 1.0
    tbhdu.header["CORRSCAL"] = 0.0
    tbhdu.header["BACKSCAL"] = 1.0

    hdulist = fits.HDUList([fits.PrimaryHDU(), tbhdu])

    hdulist.writeto(specfile, overwrite=overwrite)


def _region_filter(hdu, region, format="ds9"):
    from regions import PixCoord, PixelRegion, Region, Regions, SkyRegion

    if isinstance(region, str):
        if Path(region).exists():
            region = Regions.read(region, format=format)[0]
        else:
            region = Regions.parse(region, format=format)[0]
    elif not isinstance(region, Region):
        raise RuntimeError("'region' argument is not valid!")
    pixcoords = PixCoord(hdu.data["X"], hdu.data["Y"])
    if isinstance(region, PixelRegion):
        evt_mask = region.contains(pixcoords)
    elif isinstance(region, SkyRegion):
        w = wcs_from_header(hdu.header)
        skycoords = pixcoords.to_sky(w, origin=1)
        evt_mask = region.contains(skycoords, w)
    else:
        raise NotImplementedError
    return evt_mask


def filter_events(
    evtfile, newfile, region=None, emin=None, emax=None, format="ds9", overwrite=False
):
    r"""

    Parameters
    ----------
    evtfile : string
        The input events file to be read in.
    newfile : string
        The new event file that will be written.
    region : string or Region, optional
        The region to be used for the filtering. Default: None
    emin : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
        The minimum energy of the events to be binned in keV.
        Default is the lowest energy available.
    emax : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
        The maximum energy of the events to be binned in keV.
        Default is the highest energy available.
    format : string, optional
        The file format specifier for the region. "ds9",
        "crtf", "fits", etc. Default: "ds9"
    overwrite : boolean, optional
        Whether to overwrite an existing file with the
        same name. Default: False
    """
    with fits.open(evtfile) as f:
        hdu = f["EVENTS"]
        evt_mask = np.ones(hdu.data["ENERGY"].size, dtype="bool")
        if region is not None:
            evt_mask &= _region_filter(hdu, region, format=format)
        if emin is not None:
            emin = parse_value(emin, "keV") * 1000.0
            evt_mask &= hdu.data["ENERGY"] > emin
        if emax is not None:
            emax = parse_value(emax, "keV") * 1000.0
            evt_mask &= hdu.data["ENERGY"] < emax
        hdu.data = hdu.data[evt_mask]
        f.writeto(newfile, overwrite=overwrite)


def write_spectrum(evtfile, specfile, region=None, format="ds9", overwrite=False):
    r"""
    Bin event energies into a spectrum and write it to
    a FITS binary table. Does not do any grouping of
    channels, and will automatically determine PI or PHA.

    Parameters
    ----------
    evtfile : string
        The name of the event file to read the events from.
    specfile : string
        The name of the spectrum file to be written.
    region : string or Region, optional
        The region to be used for the filtering. Default: None
    format : string, optional
        The file format specifier for the region. "ds9",
        "crtf", "fits", etc. Default: "ds9"
    overwrite : boolean, optional
        Whether to overwrite an existing file with
        the same name. Default: False
    """
    from soxs.response import RedistributionMatrixFile

    parameters = {}
    if isinstance(evtfile, str):
        with fits.open(evtfile) as f:
            hdu = f["EVENTS"]
            evt_mask = np.ones(hdu.data["ENERGY"].size, dtype="bool")
            if region is not None:
                evt_mask &= _region_filter(hdu, region, format=format)
            spectype = hdu.header["CHANTYPE"]
            rmf = hdu.header["RESPFILE"]
            p = hdu.data[spectype][evt_mask]
            exp_time = hdu.header["EXPOSURE"]
            for key in ["RESPFILE", "ANCRFILE", "MISSION", "TELESCOP", "INSTRUME"]:
                parameters[key] = hdu.header[key]
    else:
        rmf = evtfile["rmf"]
        spectype = evtfile["channel_type"]
        p = evtfile[spectype]
        parameters["RESPFILE"] = PurePath(rmf).parts[-1]
        parameters["ANCRFILE"] = PurePath(evtfile["arf"]).parts[-1]
        parameters["TELESCOP"] = evtfile["telescope"]
        parameters["INSTRUME"] = evtfile["instrument"]
        parameters["MISSION"] = evtfile["mission"]
        exp_time = evtfile["exposure_time"]

    rmf = RedistributionMatrixFile(rmf)
    minlength = rmf.n_ch
    if rmf.cmin == 1:
        minlength += 1
    spec = np.bincount(p, minlength=minlength)
    if rmf.cmin == 1:
        spec = spec[1:]
    bins = (np.arange(rmf.n_ch) + rmf.cmin).astype("int32")

    _write_spectrum(
        bins, spec, exp_time, spectype, parameters, specfile, overwrite=overwrite
    )


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


coord_types = {"sky": ("X", "Y", 2, 3), "det": ("DETX", "DETY", 6, 7)}


def write_image(
    evt_file,
    out_file,
    coord_type="sky",
    emin=None,
    emax=None,
    overwrite=False,
    expmap_file=None,
    reblock=1,
):
    r"""
    Generate a image by binning X-ray counts and write
    it to a FITS file.

    Parameters
    ----------
    evt_file : string
        The name of the input event file to read.
    out_file : string
        The name of the image file to write.
    coord_type : string, optional
        The type of coordinate to bin into an image.
        Can be "sky" or "det". Default: "sky"
    emin : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
        The minimum energy of the photons to put in the image, in keV.
    emax : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
        The maximum energy of the photons to put in the image, in keV.
    overwrite : boolean, optional
        Whether to overwrite an existing file with
        the same name. Default: False
    expmap_file : string, optional
        Supply an exposure map file to divide this image by
        to get a flux map. Default: None
    reblock : integer, optional
        Change this value to reblock the image to larger
        pixel sizes (reblock >= 1). Only supported for
        sky coordinates. Default: 1
    """
    if emin is None:
        emin = 0.0
    else:
        emin = parse_value(emin, "keV")
    emin *= 1000.0
    if emax is None:
        emax = 100.0
    else:
        emax = parse_value(emax, "keV")
    emax *= 1000.0
    if coord_type == "det" and reblock > 1:
        raise RuntimeError(
            "Reblocking images is not supported " "for detector coordinates!"
        )
    f = fits.open(evt_file)
    e = f["EVENTS"].data["ENERGY"]
    idxs = np.logical_and(e > emin, e < emax)
    xcoord, ycoord, xcol, ycol = coord_types[coord_type]
    x = f["EVENTS"].data[xcoord][idxs]
    y = f["EVENTS"].data[ycoord][idxs]
    exp_time = f["EVENTS"].header["EXPOSURE"]
    xmin = f["EVENTS"].header[f"TLMIN{xcol}"]
    ymin = f["EVENTS"].header[f"TLMIN{ycol}"]
    xmax = f["EVENTS"].header[f"TLMAX{xcol}"]
    ymax = f["EVENTS"].header[f"TLMAX{ycol}"]
    if coord_type == "sky":
        xctr = f["EVENTS"].header[f"TCRVL{xcol}"]
        yctr = f["EVENTS"].header[f"TCRVL{ycol}"]
        xdel = f["EVENTS"].header[f"TCDLT{xcol}"] * reblock
        ydel = f["EVENTS"].header[f"TCDLT{ycol}"] * reblock
    f.close()

    nx = int(xmax - xmin) // reblock
    ny = int(ymax - ymin) // reblock

    xbins = np.linspace(xmin, xmax, nx + 1, endpoint=True)
    ybins = np.linspace(ymin, ymax, ny + 1, endpoint=True)

    H, xedges, yedges = np.histogram2d(x, y, bins=[xbins, ybins])

    if expmap_file is not None:
        if coord_type == "det":
            raise RuntimeError(
                "Cannot divide by an exposure map for images "
                "binned in detector coordinates!"
            )
        f = fits.open(expmap_file)
        if f["EXPMAP"].shape != (nx, ny):
            raise RuntimeError("Exposure map and image do not have the same shape!!")
        with np.errstate(invalid="ignore", divide="ignore"):
            H /= f["EXPMAP"].data.T
        H[np.isinf(H)] = 0.0
        H = np.nan_to_num(H)
        H[H < 0.0] = 0.0
        f.close()

    hdu = fits.PrimaryHDU(H.T)

    if coord_type == "sky":
        hdu.header["MTYPE1"] = "EQPOS"
        hdu.header["MFORM1"] = "RA,DEC"
        hdu.header["CTYPE1"] = "RA---TAN"
        hdu.header["CTYPE2"] = "DEC--TAN"
        hdu.header["CRVAL1"] = xctr
        hdu.header["CRVAL2"] = yctr
        hdu.header["CUNIT1"] = "deg"
        hdu.header["CUNIT2"] = "deg"
        hdu.header["CDELT1"] = xdel
        hdu.header["CDELT2"] = ydel
        hdu.header["CRPIX1"] = 0.5 * (nx + 1)
        hdu.header["CRPIX2"] = 0.5 * (ny + 1)
    else:
        hdu.header["CUNIT1"] = "pixel"
        hdu.header["CUNIT2"] = "pixel"

    hdu.header["EXPOSURE"] = exp_time
    hdu.name = "IMAGE"

    hdu.writeto(out_file, overwrite=overwrite)


def plot_spectrum(
    specfile,
    plot_energy=True,
    ebins=None,
    lw=2,
    xmin=None,
    xmax=None,
    ymin=None,
    ymax=None,
    xscale=None,
    yscale=None,
    label=None,
    fontsize=18,
    fig=None,
    ax=None,
    plot_counts=False,
    noerr=False,
    plot_used=False,
    **kwargs,
):
    """
    Make a quick Matplotlib plot of a convolved spectrum
    from a file. A Matplotlib figure and axis is returned.

    Parameters
    ----------
    specfile : string
        The file to be opened for plotting.
    plot_energy : boolean, optional
        Whether to plot in energy or channel space. Default is
        to plot in energy, unless the RMF for the spectrum
        cannot be found.
    ebins : NumPy array, optional
        If set, these are the energy bin edges in which the spectrum
        will be binned. If not set, the counts will be binned according
        to channel. Default: None
    lw : float, optional
        The width of the lines in the plots. Default: 2.0 px.
    xmin : float, optional
        The left-most energy (in keV) or channel to plot. Default is the
        minimum value in the spectrum.
    xmax : float, optional
        The right-most energy (in keV) or channel to plot. Default is the
        maximum value in the spectrum.
    ymin : float, optional
        The lower extent of the y-axis. By default it is set automatically.
    ymax : float, optional
        The upper extent of the y-axis. By default it is set automatically.
    xscale : string, optional
        The scaling of the x-axis of the plot. Default: "log"
    yscale : string, optional
        The scaling of the y-axis of the plot. Default: "log"
    label : string, optional
        The label of the spectrum. Default: None
    fontsize : int
        Font size for labels and axes. Default: 18
    fig : :class:`~matplotlib.figure.Figure`, optional
        A Figure instance to plot in. Default: None, one will be
        created if not provided.
    ax : :class:`~matplotlib.axes.Axes`, optional
        An Axes instance to plot in. Default: None, one will be
        created if not provided.
    plot_counts : boolean, optional
        If set to True, the counts instead of the count rate will
        be plotted. Default: False
    noerr : boolean, optional
        If True, the spectrum will be plotted without errorbars.
        Default: False
    plot_used : boolean, optional
        If set to True, only the bins which contain more than 0
        counts will be plotted. Default: False

    Returns
    -------
    A tuple of the :class:`~matplotlib.figure.Figure` and the
    :class:`~matplotlib.axes.Axes` objects.
    """
    import matplotlib.pyplot as plt

    from soxs.instrument import RedistributionMatrixFile

    f = fits.open(specfile)
    hdu = f["SPECTRUM"]
    chantype = hdu.header["CHANTYPE"]
    y = hdu.data["COUNTS"].astype("float64")
    if plot_energy:
        rmf = hdu.header.get("RESPFILE", None)
        if rmf is not None:
            rmf = RedistributionMatrixFile(rmf)
            e = 0.5 * (rmf.ebounds_data["E_MIN"] + rmf.ebounds_data["E_MAX"])
            if ebins is None:
                xmid = e
                xerr = 0.5 * (rmf.ebounds_data["E_MAX"] - rmf.ebounds_data["E_MIN"])
            else:
                xmid = 0.5 * (ebins[1:] + ebins[:-1])
                xerr = 0.5 * np.diff(ebins)
                y = np.histogram(e, ebins, weights=y)[0].astype("float64")
            xlabel = "Energy (keV)"
        else:
            raise RuntimeError(
                "Cannot find the RMF associated with this "
                "spectrum, so I cannot plot in energy!"
            )
    else:
        xmid = hdu.data[chantype]
        xerr = 0.5
        xlabel = f"Channel ({chantype})"
    dx = 2.0 * xerr
    yerr = np.sqrt(y)
    if not plot_counts:
        y /= hdu.header["EXPOSURE"]
        yerr /= hdu.header["EXPOSURE"]
    if plot_energy:
        yunit = "keV"
        y /= dx
        yerr /= dx
    else:
        yunit = "bin"
    f.close()
    if fig is None:
        fig = plt.figure(figsize=(10, 10))
    if xscale is None:
        if ax is None:
            xscale = "log"
        else:
            xscale = ax.get_xscale()
    if yscale is None:
        if ax is None:
            yscale = "log"
        else:
            yscale = ax.get_yscale()
    if ax is None:
        ax = fig.add_subplot(111)
    if plot_used:
        used = y > 0
        xmid = xmid[used]
        y = y[used]
        xerr = xerr[used]
        yerr = yerr[used]
    if noerr:
        ax.plot(xmid, y, lw=lw, label=label, **kwargs)
    else:
        ax.errorbar(xmid, y, yerr=yerr, xerr=xerr, lw=lw, label=label, **kwargs)
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel(xlabel, fontsize=fontsize)
    if plot_counts:
        ylabel = "Counts (counts/{0})"
    else:
        ylabel = "Count Rate (counts/s/{0})"
    ax.set_ylabel(ylabel.format(yunit), fontsize=fontsize)
    ax.tick_params(axis="both", labelsize=fontsize)
    return fig, ax


def plot_image(
    img_file,
    hdu="IMAGE",
    stretch="linear",
    vmin=None,
    vmax=None,
    facecolor="black",
    center=None,
    width=None,
    figsize=(10, 10),
    cmap=None,
):
    """
    Plot a FITS image created by SOXS using Matplotlib.

    Parameters
    ----------
    img_file : str
        The on-disk FITS image to plot.
    hdu : str or int, optional
        The image extension to plot. Default is "IMAGE"
    stretch : str, optional
        The stretch to apply to the colorbar scale. Options are "linear",
        "log", and "sqrt". Default: "linear"
    vmin : float, optional
        The minimum value of the colorbar. If not set, it will be the minimum
        value in the image.
    vmax : float, optional
        The maximum value of the colorbar. If not set, it will be the maximum
        value in the image.
    facecolor : str, optional
        The color of zero-valued pixels. Default: "black"
    center : array-like
        A 2-element object giving an (RA, Dec) coordinate for the center
        in degrees. If not set, the reference pixel of the image (usually
        the center) is used.
    width : float, optional
        The width of the image in degrees. If not set, the width of the
        entire image will be used.
    figsize : tuple, optional
        A 2-tuple giving the size of the image in inches, e.g. (12, 15).
        Default: (10,10)
    cmap : str, optional
        The colormap to be used. If not set, the default Matplotlib
        colormap will be used.

    Returns
    -------
    A tuple of the :class:`~matplotlib.figure.Figure` and the
    :class:`~matplotlib.axes.Axes` objects.
    """
    import matplotlib.pyplot as plt
    from astropy.visualization.wcsaxes import WCSAxes
    from astropy.wcs.utils import proj_plane_pixel_scales
    from matplotlib.colors import LogNorm, Normalize, PowerNorm

    if stretch == "linear":
        norm = Normalize(vmin=vmin, vmax=vmax)
    elif stretch == "log":
        norm = LogNorm(vmin=vmin, vmax=vmax)
    elif stretch == "sqrt":
        norm = PowerNorm(0.5, vmin=vmin, vmax=vmax)
    else:
        raise RuntimeError(f"'{stretch}' is not a valid stretch!")
    with fits.open(img_file) as f:
        hdu = f[hdu]
        w = wcs.WCS(hdu.header)
        pix_scale = proj_plane_pixel_scales(w)
        if center is None:
            center = w.wcs.crpix
        else:
            center = w.wcs_world2pix(center[0], center[1], 0)
        if width is None:
            dx_pix = 0.5 * hdu.shape[0]
            dy_pix = 0.5 * hdu.shape[1]
        else:
            dx_pix = width / pix_scale[0]
            dy_pix = width / pix_scale[1]
        fig = plt.figure(figsize=figsize)
        ax = WCSAxes(fig, [0.15, 0.1, 0.8, 0.8], wcs=w)
        fig.add_axes(ax)
        im = ax.imshow(hdu.data, norm=norm, cmap=cmap)
        ax.set_xlim(center[0] - 0.5 * dx_pix, center[0] + 0.5 * dx_pix)
        ax.set_ylim(center[1] - 0.5 * dy_pix, center[1] + 0.5 * dy_pix)
        ax.set_facecolor(facecolor)
        plt.colorbar(im)
    return fig, ax
