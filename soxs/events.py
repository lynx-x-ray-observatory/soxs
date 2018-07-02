import numpy as np
import astropy.io.fits as pyfits
import astropy.wcs as pywcs
import os
from six import string_types
from soxs.utils import mylog, parse_value, get_rot_mat, \
    downsample
from soxs.instrument_registry import instrument_registry
from tqdm import tqdm


def wcs_from_event_file(f):
    h = f["EVENTS"].header
    w = pywcs.WCS(naxis=2)
    w.wcs.crval = [h["TCRVL2"], h["TCRVL3"]]
    w.wcs.crpix = [h["TCRPX2"], h["TCRPX3"]]
    w.wcs.cdelt = [h["TCDLT2"], h["TCDLT3"]]
    w.wcs.ctype = [h["TCTYP2"], h["TCTYP3"]]
    w.wcs.cunit = [h["TCUNI2"], h["TCUNI3"]]
    return w


def write_event_file(events, parameters, filename, overwrite=False):
    from astropy.time import Time, TimeDelta
    mylog.info("Writing events to file %s." % filename)

    t_begin = Time.now()
    dt = TimeDelta(parameters["exposure_time"], format='sec')
    t_end = t_begin + dt

    col_x = pyfits.Column(name='X', format='D', unit='pixel', array=events["xpix"])
    col_y = pyfits.Column(name='Y', format='D', unit='pixel', array=events["ypix"])
    col_e = pyfits.Column(name='ENERGY', format='E', unit='eV', array=events["energy"]*1000.)
    col_dx = pyfits.Column(name='DETX', format='D', unit='pixel', array=events["detx"])
    col_dy = pyfits.Column(name='DETY', format='D', unit='pixel', array=events["dety"])
    col_id = pyfits.Column(name='CCD_ID', format='D', unit='pixel', array=events["ccd_id"])

    chantype = parameters["channel_type"]
    if chantype == "PHA":
        cunit = "adu"
    elif chantype == "PI":
        cunit = "Chan"
    col_ch = pyfits.Column(name=chantype.upper(), format='1J', unit=cunit, array=events[chantype])

    col_t = pyfits.Column(name="TIME", format='1D', unit='s', array=events['time'])

    cols = [col_e, col_x, col_y, col_ch, col_t, col_dx, col_dy, col_id]

    coldefs = pyfits.ColDefs(cols)
    tbhdu = pyfits.BinTableHDU.from_columns(coldefs)
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
    tbhdu.header["TLMAX2"] = 2.0*parameters["num_pixels"]+0.5
    tbhdu.header["TLMAX3"] = 2.0*parameters["num_pixels"]+0.5
    tbhdu.header["TLMIN4"] = parameters["chan_lim"][0]
    tbhdu.header["TLMAX4"] = parameters["chan_lim"][1]
    tbhdu.header["TLMIN6"] = -0.5*parameters["num_pixels"]
    tbhdu.header["TLMAX6"] = 0.5*parameters["num_pixels"]
    tbhdu.header["TLMIN7"] = -0.5*parameters["num_pixels"]
    tbhdu.header["TLMAX7"] = 0.5*parameters["num_pixels"]
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
    tbhdu.header["RESPFILE"] = os.path.split(parameters["rmf"])[-1]
    tbhdu.header["PHA_BINS"] = parameters["nchan"]
    tbhdu.header["ANCRFILE"] = os.path.split(parameters["arf"])[-1]
    tbhdu.header["CHANTYPE"] = parameters["channel_type"]
    tbhdu.header["MISSION"] = parameters["mission"]
    tbhdu.header["TELESCOP"] = parameters["telescope"]
    tbhdu.header["INSTRUME"] = parameters["instrument"]
    tbhdu.header["RA_PNT"] = parameters["sky_center"][0]
    tbhdu.header["DEC_PNT"] = parameters["sky_center"][1]
    tbhdu.header["ROLL_PNT"] = parameters["roll_angle"]
    tbhdu.header["AIMPT_X"] = parameters["aimpt_coords"][0]
    tbhdu.header["AIMPT_Y"] = parameters["aimpt_coords"][1]
    if parameters["dither_params"]["dither_on"]:
        tbhdu.header["DITHXAMP"] = parameters["dither_params"]["x_amp"]
        tbhdu.header["DITHYAMP"] = parameters["dither_params"]["y_amp"]
        tbhdu.header["DITHXPER"] = parameters["dither_params"]["x_period"]
        tbhdu.header["DITHYPER"] = parameters["dither_params"]["y_period"]

    start = pyfits.Column(name='START', format='1D', unit='s',
                          array=np.array([0.0]))
    stop = pyfits.Column(name='STOP', format='1D', unit='s',
                         array=np.array([parameters["exposure_time"]]))

    tbhdu_gti = pyfits.BinTableHDU.from_columns([start,stop])
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

    hdulist = [pyfits.PrimaryHDU(), tbhdu, tbhdu_gti]

    pyfits.HDUList(hdulist).writeto(filename, overwrite=overwrite)


def parse_region_args(rtype, args, dx, dy):
    if rtype == "Box":
        xctr, yctr, xw, yw = args
        new_args = [xctr + dx, yctr + dy, xw, yw]
    elif rtype == "Circle":
        xctr, yctr, radius = args
        new_args = [xctr + dx, yctr + dx, radius]
    elif rtype == "Polygon":
        new_args = [[x + dx for x in args[0]],
                    [y + dy for y in args[1]]]
    else:
        raise NotImplementedError
    return new_args


def make_exposure_map(event_file, expmap_file, energy, weights=None,
                      asol_file=None, normalize=True, overwrite=False,
                      reblock=1, nhistx=16, nhisty=16, order=1):
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
        Whether or not to overwrite an existing file. Default: False
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
    import pyregion._region_filter as rfilter
    from scipy.ndimage.interpolation import rotate, shift
    from soxs.instrument import AuxiliaryResponseFile, perform_dither
    if isinstance(energy, np.ndarray) and weights is None:
        raise RuntimeError("Must supply a single value for the energy if "
                           "you do not supply weights!")
    if not isinstance(energy, np.ndarray):
        energy = parse_value(energy, "keV")
    f_evt = pyfits.open(event_file)
    hdu = f_evt["EVENTS"]
    arf = AuxiliaryResponseFile(hdu.header["ANCRFILE"])
    exp_time = hdu.header["EXPOSURE"]
    nx = int(hdu.header["TLMAX2"]-0.5)//2
    ny = int(hdu.header["TLMAX3"]-0.5)//2
    ra0 = hdu.header["TCRVL2"]
    dec0 = hdu.header["TCRVL3"]
    xdel = hdu.header["TCDLT2"]
    ydel = hdu.header["TCDLT3"]
    x0 = hdu.header["TCRPX2"]
    y0 = hdu.header["TCRPX3"]
    xdet0 = 0.5*(2*nx+1)
    ydet0 = 0.5*(2*ny+1)
    xaim = hdu.header.get("AIMPT_X", 0.0)
    yaim = hdu.header.get("AIMPT_Y", 0.0)
    roll = hdu.header["ROLL_PNT"]
    instr = instrument_registry[hdu.header["INSTRUME"].lower()]
    dither_params = {}
    if "DITHXAMP" in hdu.header:
        dither_params["x_amp"] = hdu.header["DITHXAMP"]
        dither_params["y_amp"] = hdu.header["DITHYAMP"]
        dither_params["x_period"] = hdu.header["DITHXPER"]
        dither_params["y_period"] = hdu.header["DITHYPER"]
        dither_params["plate_scale"] = ydel*3600.0
        dither_params["dither_on"] = True
    else:
        dither_params["dither_on"] = False
    f_evt.close()

    # Create time array for aspect solution
    dt = 1.0 # Seconds
    t = np.arange(0.0, exp_time+dt, dt)

    # Construct WCS
    w = pywcs.WCS(naxis=2)
    w.wcs.crval = [ra0, dec0]
    w.wcs.crpix = [x0, y0]
    w.wcs.cdelt = [xdel, ydel]
    w.wcs.ctype = ["RA---TAN","DEC--TAN"]
    w.wcs.cunit = ["deg"]*2

    # Create aspect solution if we had dithering.
    # otherwise just set the offsets to zero
    if dither_params["dither_on"]:
        x_off, y_off = perform_dither(t, dither_params)
        # Make the aspect histogram
        x_amp = dither_params["x_amp"]/dither_params["plate_scale"]
        y_amp = dither_params["y_amp"]/dither_params["plate_scale"]
        x_edges = np.linspace(-x_amp, x_amp, nhistx+1, endpoint=True)
        y_edges = np.linspace(-y_amp, y_amp, nhisty+1, endpoint=True)
        asphist = np.histogram2d(x_off, y_off, (x_edges, y_edges))[0]
        asphist *= dt
        x_mid = 0.5*(x_edges[1:]+x_edges[:-1])/reblock
        y_mid = 0.5*(y_edges[1:]+y_edges[:-1])/reblock

    # Determine the effective area
    eff_area = arf.interpolate_area(energy).value
    if weights is not None:
        eff_area = np.average(eff_area, weights=weights)

    if instr["chips"] is None:
        rtypes = ["Box"]
        args = [[0.0, 0.0, instr["num_pixels"], instr["num_pixels"]]]
    else:
        rtypes = []
        args = []
        for i, chip in enumerate(instr["chips"]):
            rtypes.append(chip[0])
            args.append(np.array(chip[1:]))

    tmpmap = np.zeros((2*nx, 2*ny))

    for rtype, arg in zip(rtypes, args):
        rfunc = getattr(rfilter, rtype)
        new_args = parse_region_args(rtype, arg, xdet0-xaim-1.0, ydet0-yaim-1.0)
        r = rfunc(*new_args)
        tmpmap += r.mask(tmpmap).astype("float64")

    tmpmap = downsample(tmpmap, reblock)

    if dither_params["dither_on"]:
        expmap = np.zeros(tmpmap.shape)
        niter = nhistx*nhisty
        pbar = tqdm(leave=True, total=niter, desc="Creating exposure map ")
        for i in range(nhistx):
            for j in range(nhisty):
                expmap += shift(tmpmap, (x_mid[i], y_mid[j]), order=order)*asphist[i, j]
            pbar.update(nhisty)
        pbar.close()
    else:
        expmap = tmpmap*exp_time

    expmap *= eff_area
    if normalize:
        expmap /= exp_time

    if roll != 0.0:
        rotate(expmap, roll, output=expmap, reshape=False)

    expmap[expmap < 0.0] = 0.0

    map_header = {"EXPOSURE": exp_time,
                  "MTYPE1": "EQPOS",
                  "MFORM1": "RA,DEC",
                  "CTYPE1": "RA---TAN",
                  "CTYPE2": "DEC--TAN",
                  "CRVAL1": ra0,
                  "CRVAL2": dec0,
                  "CUNIT1": "deg",
                  "CUNIT2": "deg",
                  "CDELT1": xdel*reblock,
                  "CDELT2": ydel*reblock,
                  "CRPIX1": 0.5*(2.0*nx//reblock+1),
                  "CRPIX2": 0.5*(2.0*ny//reblock+1)}

    map_hdu = pyfits.ImageHDU(expmap, header=pyfits.Header(map_header))
    map_hdu.name = "EXPMAP"
    map_hdu.writeto(expmap_file, overwrite=overwrite)

    if asol_file is not None:

        if dither_params["dither_on"]:

            det = np.array([x_off, y_off])

            pix = np.dot(get_rot_mat(roll).T, det)

            ra, dec = w.wcs_pix2world(pix[0,:]+x0, pix[1,:]+y0, 1)

            col_t = pyfits.Column(name='time', format='D', unit='s', array=t)
            col_ra = pyfits.Column(name='ra', format='D', unit='deg', array=ra)
            col_dec = pyfits.Column(name='dec', format='D', unit='deg', array=dec)

            coldefs = pyfits.ColDefs([col_t, col_ra, col_dec])
            tbhdu = pyfits.BinTableHDU.from_columns(coldefs)
            tbhdu.name = "ASPSOL"
            tbhdu.header["EXPOSURE"] = exp_time

            hdulist = [pyfits.PrimaryHDU(), tbhdu]

            pyfits.HDUList(hdulist).writeto(asol_file, overwrite=overwrite)

        else:

            mylog.warning("Refusing to write an aspect solution file because "
                          "there was no dithering.")


def _write_spectrum(bins, spec, exp_time, spectype, parameters,
                    specfile, overwrite=False):

    col1 = pyfits.Column(name='CHANNEL', format='1J', array=bins)
    col2 = pyfits.Column(name=spectype.upper(), format='1D', array=bins.astype("float64"))
    col3 = pyfits.Column(name='COUNTS', format='1J', array=spec.astype("int32"))
    col4 = pyfits.Column(name='COUNT_RATE', format='1D', array=spec/exp_time)

    coldefs = pyfits.ColDefs([col1, col2, col3, col4])

    tbhdu = pyfits.BinTableHDU.from_columns(coldefs)
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

    hdulist = pyfits.HDUList([pyfits.PrimaryHDU(), tbhdu])

    hdulist.writeto(specfile, overwrite=overwrite)


def write_spectrum(evtfile, specfile, overwrite=False):
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
    overwrite : boolean, optional
        Whether or not to overwrite an existing file with 
        the same name. Default: False
    """
    from soxs.instrument import RedistributionMatrixFile
    parameters = {}
    if isinstance(evtfile, string_types):
        f = pyfits.open(evtfile)
        spectype = f["EVENTS"].header["CHANTYPE"]
        rmf = f["EVENTS"].header["RESPFILE"]
        p = f["EVENTS"].data[spectype]
        exp_time = f["EVENTS"].header["EXPOSURE"]
        for key in ["RESPFILE", "ANCRFILE", "MISSION", "TELESCOP", "INSTRUME"]:
            parameters[key] = f["EVENTS"].header[key]
        f.close()
    else:
        rmf = evtfile["rmf"]
        spectype = evtfile["channel_type"]
        p = evtfile[spectype]
        parameters["RESPFILE"] = os.path.split(rmf)[-1]
        parameters["ANCRFILE"] = os.path.split(evtfile["arf"])[-1]
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
    bins = (np.arange(rmf.n_ch)+rmf.cmin).astype("int32")

    _write_spectrum(bins, spec, exp_time, spectype, parameters,
                    specfile, overwrite=overwrite)


def write_radial_profile(evt_file, out_file, ctr, rmin,
                         rmax, nbins, ctr_type="celestial",
                         emin=None, emax=None, expmap_file=None,
                         overwrite=False):
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
        Whether or not to overwrite an existing file with the 
        same name. Default: False
    expmap_file : string, optional
        Supply an exposure map file to determine fluxes. 
        Default: None
    """
    import astropy.wcs as pywcs
    rmin = parse_value(rmin, "arcsec")
    rmax = parse_value(rmax, "arcsec")
    f = pyfits.open(evt_file)
    hdu = f["EVENTS"]
    orig_dx = hdu.header["TCDLT3"]
    e = hdu.data["ENERGY"]
    if emin is None:
        emin = e.min()
    else:
        emin = parse_value(emin, "keV")
        emin *= 1000.
    if emax is None:
        emax = e.max()
    else:
        emax = parse_value(emax, "keV")
        emax *= 1000.
    idxs = np.logical_and(e > emin, e < emax)
    x = hdu.data["X"][idxs]
    y = hdu.data["Y"][idxs]
    exp_time = hdu.header["EXPOSURE"]
    w = wcs_from_event_file(f)
    dtheta = np.abs(w.wcs.cdelt[1])*3600.0
    f.close()

    if ctr_type == "celestial":
        ctr = w.all_world2pix(ctr[0], ctr[1], 1)

    r = np.sqrt((x-ctr[0])**2+(y-ctr[1])**2)
    rr = np.linspace(rmin/dtheta, rmax/dtheta, nbins+1)
    C = np.histogram(r, bins=rr)[0]
    rbin = rr*dtheta
    rmid = 0.5*(rbin[1:]+rbin[:-1])

    A = np.pi*(rbin[1:]**2-rbin[:-1]**2)

    Cerr = np.sqrt(C)

    R = C/exp_time
    Rerr = Cerr/exp_time

    S = R/A
    Serr = Rerr/A

    col1 = pyfits.Column(name='RLO', format='D', unit='arcsec', array=rbin[:-1])
    col2 = pyfits.Column(name='RHI', format='D', unit='arcsec', array=rbin[1:])
    col3 = pyfits.Column(name='RMID', format='D', unit='arcsec', array=rmid)
    col4 = pyfits.Column(name='AREA', format='D', unit='arcsec**2', array=A)
    col5 = pyfits.Column(name='NET_COUNTS', format='D', unit='count', array=C)
    col6 = pyfits.Column(name='NET_ERR', format='D', unit='count', array=Cerr)
    col7 = pyfits.Column(name='NET_RATE', format='D', unit='count/s', array=R)
    col8 = pyfits.Column(name='ERR_RATE', format='D', unit='count/s', array=Rerr)
    col9 = pyfits.Column(name='SUR_BRI', format='D', unit='count/s/arcsec**2', array=S)
    col10 = pyfits.Column(name='SUR_BRI_ERR', format='1D', unit='count/s/arcsec**2', array=Serr)

    coldefs = [col1, col2, col3, col4, col5, col6, col7, col8, col9, col10]

    if expmap_file is not None:
        f = pyfits.open(expmap_file)
        ehdu = f["EXPMAP"]
        wexp = pywcs.WCS(header=ehdu.header)
        cel = w.all_pix2world(ctr[0], ctr[1], 1)
        ectr = wexp.all_world2pix(cel[0], cel[1], 1)
        exp = ehdu.data[:,:]
        nx, ny = exp.shape
        reblock = ehdu.header["CDELT2"]/orig_dx
        x, y = np.mgrid[1:nx+1,1:ny+1]
        r = np.sqrt((x-ectr[0])**2 + (y-ectr[1])**2)
        f.close()
        E = np.histogram(r, bins=rr/reblock, weights=exp)[0] / np.histogram(r, bins=rr/reblock)[0]
        with np.errstate(invalid='ignore', divide='ignore'):
            F = R/E
            Ferr = Rerr/E
        SF = F/A
        SFerr = Ferr/A
        col11 = pyfits.Column(name='MEAN_SRC_EXP', format='D', unit='cm**2', array=E)
        col12 = pyfits.Column(name='NET_FLUX', format='D', unit='count/s/cm**2', array=F)
        col13 = pyfits.Column(name='NET_FLUX_ERR', format='D', unit='count/s/cm**2', array=Ferr)
        col14 = pyfits.Column(name='SUR_FLUX', format='D', unit='count/s/cm**2/arcsec**2', array=SF)
        col15 = pyfits.Column(name='SUR_FLUX_ERR', format='D', unit='count/s/cm**2/arcsec**2', array=SFerr)
        coldefs += [col11, col12, col13, col14, col15]

    tbhdu = pyfits.BinTableHDU.from_columns(pyfits.ColDefs(coldefs))
    tbhdu.name = "PROFILE"

    hdulist = pyfits.HDUList([pyfits.PrimaryHDU(), tbhdu])

    hdulist.writeto(out_file, overwrite=overwrite)

coord_types = {"sky": ("X", "Y", 2, 3),
               "det": ("DETX", "DETY", 6, 7)}


def write_image(evt_file, out_file, coord_type='sky', emin=None, emax=None,
                overwrite=False, expmap_file=None, reblock=1):
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
        The minimum energy of the photons to put in the 
        image, in keV.
    emax : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
        The maximum energy of the photons to put in the 
        image, in keV.
    overwrite : boolean, optional
        Whether or not to overwrite an existing file with 
        the same name. Default: False
    expmap_file : string, optional
        Supply an exposure map file to divide this image by
        to get a flux map. Default: None
    reblock : integer, optional
        Change this value to reblock the image to larger 
        pixel sizes (reblock >= 1). Only supported for
        sky coordinates. Default: 1
    """
    if coord_type == "det" and reblock > 1:
        raise RuntimeError("Reblocking images is not supported "
                           "for detector coordinates!")
    f = pyfits.open(evt_file)
    e = f["EVENTS"].data["ENERGY"]
    if emin is None:
        emin = e.min()
    else:
        emin = parse_value(emin, "keV")
        emin *= 1000.
    if emax is None:
        emax = e.max()
    else:
        emax = parse_value(emax, "keV")
        emax *= 1000.
    idxs = np.logical_and(e > emin, e < emax)
    xcoord, ycoord, xcol, ycol = coord_types[coord_type]
    x = f["EVENTS"].data[xcoord][idxs]
    y = f["EVENTS"].data[ycoord][idxs]
    exp_time = f["EVENTS"].header["EXPOSURE"]
    xmin = f["EVENTS"].header["TLMIN%d" % xcol]
    ymin = f["EVENTS"].header["TLMIN%d" % ycol]
    xmax = f["EVENTS"].header["TLMAX%d" % xcol]
    ymax = f["EVENTS"].header["TLMAX%d" % ycol]
    if coord_type == 'sky':
        xctr = f["EVENTS"].header["TCRVL%d" % xcol]
        yctr = f["EVENTS"].header["TCRVL%d" % ycol]
        xdel = f["EVENTS"].header["TCDLT%d" % xcol]*reblock
        ydel = f["EVENTS"].header["TCDLT%d" % ycol]*reblock
    f.close()

    nx = int(xmax-xmin)//reblock
    ny = int(ymax-ymin)//reblock

    xbins = np.linspace(xmin, xmax, nx+1, endpoint=True)
    ybins = np.linspace(ymin, ymax, ny+1, endpoint=True)

    H, xedges, yedges = np.histogram2d(x, y, bins=[xbins, ybins])

    if expmap_file is not None:
        if coord_type == "det":
            raise RuntimeError("Cannot divide by an exposure map for images "
                               "binned in detector coordinates!")
        f = pyfits.open(expmap_file)
        if f["EXPMAP"].shape != (nx, ny):
            raise RuntimeError("Exposure map and image do not have the same shape!!")
        with np.errstate(invalid='ignore', divide='ignore'):
            H /= f["EXPMAP"].data.T
        H[np.isinf(H)] = 0.0
        H = np.nan_to_num(H)
        H[H < 0.0] = 0.0
        f.close()

    hdu = pyfits.PrimaryHDU(H.T)

    if coord_type == 'sky':
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
        hdu.header["CRPIX1"] = 0.5*(nx+1)
        hdu.header["CRPIX2"] = 0.5*(ny+1)
    else:
        hdu.header["CUNIT1"] = "pixel"
        hdu.header["CUNIT2"] = "pixel"

    hdu.header["EXPOSURE"] = exp_time

    hdu.writeto(out_file, overwrite=overwrite)


def plot_spectrum(specfile, plot_energy=True, lw=2, xmin=None, xmax=None,
                  ymin=None, ymax=None, xscale=None, yscale=None, 
                  label=None, fontsize=18, fig=None, ax=None, 
                  plot_counts=False, **kwargs):
    """
    Make a quick Matplotlib plot of a convolved spectrum
    from a file. A Matplotlib figure and axis is returned.

    Parameters
    ----------
    specfile : string
        The file to be opened for plotting.
    figsize : tuple of integers, optional
        The size of the figure on both sides in inches.
        Default: (10,10)
    plot_energy : boolean, optional
        Whether to plot in energy or channel space. Default is
        to plot in energy, unless the RMF for the spectrum
        cannot be found. 
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

    Returns
    -------
    A tuple of the :class:`~matplotlib.figure.Figure` and the :class:`~matplotlib.axes.Axes` objects.
    """
    import matplotlib.pyplot as plt
    from soxs.instrument import RedistributionMatrixFile
    f = pyfits.open(specfile)
    hdu = f["SPECTRUM"]
    chantype = hdu.header["CHANTYPE"]
    rmf = hdu.header.get("RESPFILE", None)
    xerr = None
    if plot_energy:
        if rmf is not None:
            rmf = RedistributionMatrixFile(rmf)
            x = 0.5*(rmf.ebounds_data["E_MIN"]+rmf.ebounds_data["E_MAX"])
            xerr = 0.5*(rmf.ebounds_data["E_MAX"]-rmf.ebounds_data["E_MIN"])
            xlabel = "Energy (keV)"
        else:
            raise RuntimeError("Cannot find the RMF associated with this "
                               "spectrum, so I cannot plot in energy!")
    else:
        x = hdu.data[chantype]
        xlabel = "Channel (%s)" % chantype
    if plot_counts:
        y = hdu.data["COUNTS"].astype("float64")
        yerr = np.sqrt(y)
    else:
        if "COUNT_RATE" in hdu.columns.names:
            y = hdu.data["COUNT_RATE"]
        else:
            y = hdu.data["COUNTS"]/hdu.header["EXPOSURE"]
        yerr = np.sqrt(hdu.data["COUNTS"])/hdu.header["EXPOSURE"]
    if plot_energy:
        yunit = "keV"
        y /= 2.0*xerr
        yerr /= 2.0*xerr
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
    ax.errorbar(x, y, yerr=yerr, xerr=xerr, lw=lw, label=label, **kwargs)
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel(xlabel, fontsize=fontsize)
    if plot_counts:
        ylabel = "Counts (counts/%s)"
    else:
        ylabel = "Count Rate (counts/s/%s)"
    ax.set_ylabel(ylabel % yunit, fontsize=fontsize)
    ax.tick_params(axis='both', labelsize=fontsize)
    return fig, ax