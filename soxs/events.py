import numpy as np
import astropy.io.fits as pyfits
import astropy.wcs as pywcs
import os

def wcs_from_event_file(f):
    h = f["EVENTS"].header
    w = pywcs.WCS(naxis=2)
    w.wcs.crval = [h["TCRVL2"], h["TCRVL3"]]
    w.wcs.crpix = [h["TCRPX2"], h["TCRPX3"]]
    w.wcs.cdelt = [h["TCDLT2"], h["TCDLT3"]]
    w.wcs.ctype = [h["TCTYP2"], h["TCTYP3"]]
    w.wcs.cunit = [h["TCUNI2"], h["TCUNI3"]]
    return w

def write_event_file(events, parameters, filename, clobber=False):
    from astropy.time import Time, TimeDelta

    t_begin = Time.now()
    dt = TimeDelta(parameters["exposure_time"], format='sec')
    t_end = t_begin + dt

    col_x = pyfits.Column(name='X', format='D', unit='pixel', array=events["xpix"])
    col_y = pyfits.Column(name='Y', format='D', unit='pixel', array=events["ypix"])
    col_e = pyfits.Column(name='ENERGY', format='E', unit='eV', array=events["energy"]*1000.)
    col_cx = pyfits.Column(name='CHIPX', format='D', unit='pixel', array=events["chipx"])
    col_cy = pyfits.Column(name='CHIPY', format='D', unit='pixel', array=events["chipy"])
    col_dx = pyfits.Column(name='DETX', format='D', unit='pixel', array=events["detx"])
    col_dy = pyfits.Column(name='DETY', format='D', unit='pixel', array=events["dety"])

    chantype = parameters["channel_type"]
    if chantype == "PHA":
        cunit = "adu"
    elif chantype == "PI":
        cunit = "Chan"
    col_ch = pyfits.Column(name=chantype.upper(), format='1J', unit=cunit, array=events[chantype])

    col_t = pyfits.Column(name="TIME", format='1D', unit='s', array=events['time'])

    cols = [col_e, col_x, col_y, col_ch, col_t, col_cx, col_cy, col_dx, col_dy]

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
    tbhdu.header["TLMAX2"] = parameters["num_pixels"]+0.5
    tbhdu.header["TLMAX3"] = parameters["num_pixels"]+0.5
    tbhdu.header["TLMIN4"] = parameters["chan_lim"][0]
    tbhdu.header["TLMAX4"] = parameters["chan_lim"][1]
    tbhdu.header["TLMIN6"] = 0.5
    tbhdu.header["TLMAX6"] = parameters["num_pixels"]+0.5
    tbhdu.header["TLMIN7"] = 0.5
    tbhdu.header["TLMAX7"] = parameters["num_pixels"]+0.5
    tbhdu.header["TLMIN8"] = 1.0-parameters["pix_center"][0]
    tbhdu.header["TLMAX8"] = parameters["num_pixels"]-parameters["pix_center"][0]
    tbhdu.header["TLMIN9"] = 1.0-parameters["pix_center"][1]
    tbhdu.header["TLMAX9"] = parameters["num_pixels"]-parameters["pix_center"][1]
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

    pyfits.HDUList(hdulist).writeto(filename, clobber=clobber)

def write_spectrum(evtfile, specfile, clobber=False):
    r"""
    Bin event energies into a spectrum and write it to a FITS binary table.
    Does not do any grouping of channels, and will automatically determine
    PI or PHA. 

    Parameters
    ----------
    evtfile : string
        The name of the event file to read the events from. 
    specfile : string
        The name of the spectrum file to be written.
    clobber : boolean, optional
        Whether or not to clobber an existing file with the same name.
        Default: False
    """
    from soxs.instrument import RedistributionMatrixFile
    f = pyfits.open(evtfile)
    spectype = f["EVENTS"].header["CHANTYPE"]
    rmf = RedistributionMatrixFile(f["EVENTS"].header["RESPFILE"])
    minlength = rmf.n_ch
    if rmf.cmin == 1:
        minlength += 1
    spec = np.bincount(f["EVENTS"].data[spectype], minlength=minlength)
    if rmf.cmin == 1:
        spec = spec[1:]
    bins = (np.arange(rmf.n_ch)+rmf.cmin).astype("int32")

    col1 = pyfits.Column(name='CHANNEL', format='1J', array=bins)
    col2 = pyfits.Column(name=spectype.upper(), format='1D', array=bins.astype("float64"))
    col3 = pyfits.Column(name='COUNTS', format='1J', array=spec.astype("int32"))
    col4 = pyfits.Column(name='COUNT_RATE', format='1D', array=spec/f["EVENTS"].header["EXPOSURE"])

    coldefs = pyfits.ColDefs([col1, col2, col3, col4])

    tbhdu = pyfits.BinTableHDU.from_columns(coldefs)
    tbhdu.name = "SPECTRUM"

    tbhdu.header["DETCHANS"] = spec.size
    tbhdu.header["TOTCTS"] = spec.sum()
    tbhdu.header["EXPOSURE"] = f["EVENTS"].header["EXPOSURE"]
    tbhdu.header["LIVETIME"] = f["EVENTS"].header["EXPOSURE"]
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
        tbhdu.header[key] = f["EVENTS"].header[key]
    tbhdu.header["AREASCAL"] = 1.0
    tbhdu.header["CORRSCAL"] = 0.0
    tbhdu.header["BACKSCAL"] = 1.0

    f.close()

    hdulist = pyfits.HDUList([pyfits.PrimaryHDU(), tbhdu])

    hdulist.writeto(specfile, clobber=clobber)

def write_radial_profile(evt_file, out_file, ctr, rmin, 
                         rmax, nbins, ctr_type="celestial", 
                         emin=None, emax=None, clobber=False):
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
    rmin : float
        The minimum radius of the profile, in arcseconds. 
    rmax : float
        The maximum radius of the profile, in arcseconds.
    nbins : integer
        The number of bins in the profile.
    ctr_type : string, optional
        The type of center coordinate. Either "celestial" for (RA, Dec)
        coordinates (the default), or "physical" for pixel coordinates.
    emin : float
        The minimum energy of the events to be binned in keV. Default
        is the lowest energy available.
    emax : float
        The maximum energy of the events to be binned in keV. Default
        is the highest energy available.
    clobber : boolean, optional
        Whether or not to clobber an existing file with the same name.
        Default: False
    """
    f = pyfits.open(evt_file)
    e = f["EVENTS"].data["ENERGY"]
    if emin is None:
        emin = e.min()
    else:
        emin *= 1000.
    if emax is None:
        emax = e.max()
    else:
        emax *= 1000.
    idxs = np.logical_and(e > emin, e < emax)
    x = f["EVENTS"].data["X"][idxs]
    y = f["EVENTS"].data["Y"][idxs]
    exp_time = f["EVENTS"].header["EXPOSURE"]
    w = wcs_from_event_file(f)
    dtheta = np.abs(w.wcs.cdelt[1])*3600.0
    f.close()

    if ctr_type == "celestial":
        ctr = w.all_world2pix(ctr[0], ctr[1], 1)

    r = np.sqrt((x-ctr[0])**2+(y-ctr[1])**2)
    rbin = np.linspace(rmin/dtheta, rmax/dtheta, nbins+1)
    C, _ = np.histogram(r, bins=rbin)
    rbin *= dtheta
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

    coldefs = pyfits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8, col9, col10])

    tbhdu = pyfits.BinTableHDU.from_columns(coldefs)
    tbhdu.name = "PROFILE"

    hdulist = pyfits.HDUList([pyfits.PrimaryHDU(), tbhdu])

    hdulist.writeto(out_file, clobber=clobber)

coord_types = {"sky": ("X", "Y", 2, 3),
               "chip": ("CHIPX", "CHIPY", 6, 7),
               "det": ("DETX", "DETY", 8, 9)}

def write_image(evt_file, out_file, coord_type='sky', emin=None, emax=None, 
                clobber=False):
    r"""
    Generate a image by binning X-ray counts and write it to a FITS file.

    Parameters
    ----------
    evt_file : string
        The name of the input event file to read.
    out_file : string
        The name of the image file to write.
    coord_type : string, optional
        The type of coordinate to bin into an image. Can be "sky", "det",
        or "chip". Default: "sky"
    emin : float, optional
        The minimum energy of the photons to put in the image, in keV.
    emax : float, optional
        The maximum energy of the photons to put in the image, in keV.
    clobber : boolean, optional
        Whether or not to clobber an existing file with the same name.
        Default: False
    """
    f = pyfits.open(evt_file)
    e = f["EVENTS"].data["ENERGY"]
    if emin is None:
        emin = e.min()
    else:
        emin *= 1000.
    if emax is None:
        emax = e.max()
    else:
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
        xdel = f["EVENTS"].header["TCDLT%d" % xcol]
        ydel = f["EVENTS"].header["TCDLT%d" % ycol]
        xpix = f["EVENTS"].header["TCRPX%d" % xcol]
        ypix = f["EVENTS"].header["TCRPX%d" % ycol]
    f.close()

    nx = int(xmax-xmin)
    ny = int(ymax-ymin)

    xbins = np.linspace(xmin, xmax, nx+1, endpoint=True)
    ybins = np.linspace(ymin, ymax, ny+1, endpoint=True)

    H, xedges, yedges = np.histogram2d(x, y, bins=[xbins, ybins])

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
        hdu.header["CRPIX1"] = xpix
        hdu.header["CRPIX2"] = ypix
    else:
        hdu.header["CUNIT1"] = "pixel"
        hdu.header["CUNIT2"] = "pixel"

    hdu.header["EXPOSURE"] = exp_time

    hdu.writeto(out_file, clobber=clobber)
