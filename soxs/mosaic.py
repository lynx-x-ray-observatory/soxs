from soxs.instrument import instrument_simulator
from soxs.events import write_image, make_exposure_map, wcs_from_header
from soxs.utils import mylog
import numpy as np
from astropy.io import fits, ascii
from astropy.table import Table


def make_mosaic_events(pointing_list, input_source, out_prefix, exp_time, 
                       instrument, overwrite=False, instr_bkgnd=True,
                       foreground=True, ptsrc_bkgnd=True, bkgnd_file=None, 
                       no_dither=False, dither_params=None, subpixel_res=False, 
                       aimpt_shift=None, prng=None):
    """
    Observe a source from many different pointings.

    Parameters
    ----------
    pointing_list : list of tuples or str
        Either a list of tuples or a two-column ASCII table, containing 
        RA and Dec pointings for each mock observation.
    input_source : string
        The path to the SIMPUT catalog file which contains the input
        source(s).
    out_prefix : string
        The prefix for the event files which will be generated. 
    exp_time : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The exposure time in seconds.
    instrument : string
        The name of the instrument to use, which picks an instrument
        specification from the instrument registry.
    overwrite : boolean, optional
        Whether to overwrite an existing file with the same name.
        Default: False
    instr_bkgnd : boolean, optional
        Whether to include the instrumental/particle background. 
        Default: True
    foreground : boolean, optional
        Whether to include the local foreground. 
        Default: True
    ptsrc_bkgnd : boolean, optional
        Whether to include the point-source background. 
        Default: True
    bkgnd_file : string, optional
        If set, backgrounds will be loaded from this file and not generated
        on the fly. Default: None
    no_dither : boolean, optional
        If True, turn off dithering entirely. Default: False
    dither_params : array-like of floats, optional
        The parameters to use to control the size and period of the dither
        pattern. The first two numbers are the dither amplitude in x and y
        detector coordinates in arcseconds, and the second two numbers are
        the dither period in x and y detector coordinates in seconds. 
        Default: [8.0, 8.0, 1000.0, 707.0].
    subpixel_res: boolean, optional
        If True, event positions are not randomized within the pixels 
        within which they are detected. Default: False
    aimpt_shift : array-like, optional
        A two-float array-like object which shifts the aimpoint on the 
        detector from the nominal position. Units are in arcseconds.
        Default: None, which results in no shift from the nominal aimpoint. 
    prng : :class:`~numpy.random.RandomState` object, integer, or None
        A pseudo-random number generator. Typically will only 
        be specified if you have a reason to generate the same 
        set of random numbers, such as for a test. Default is None, 
        which sets the seed based on the system time. 
    """
    if isinstance(pointing_list, str):
        t = ascii.read(pointing_list, format='commented_header', guess=False,
                       header_start=0, delimiter="\t")
    elif not isinstance(pointing_list, Table):
        t = Table(np.array(pointing_list), names=["ra", "dec"])
    out_list = []
    for i, row in enumerate(t):
        out_file = f"{out_prefix}_{i}_evt.fits"
        out_list.append(out_file)
        instrument_simulator(input_source, out_file, exp_time, instrument,
                             (row["ra"], row["dec"]), overwrite=overwrite,
                             instr_bkgnd=instr_bkgnd, foreground=foreground,
                             ptsrc_bkgnd=ptsrc_bkgnd, bkgnd_file=bkgnd_file,
                             no_dither=no_dither, dither_params=dither_params,
                             subpixel_res=subpixel_res, aimpt_shift=aimpt_shift,
                             prng=prng)
    t["evtfile"] = out_list
    outfile = f"{out_prefix}_event_mosaic.dat"
    mylog.info(f"Writing mosaic information to {outfile}.")
    t.write(outfile, overwrite=overwrite, delimiter="\t",
            format='ascii.commented_header')
    return outfile


def _combine_events(eventfiles, wcs_out, shape_out, outfile, overwrite=False):
    x = []
    y = []
    e = []
    ch = []
    t = []
    chantype = ""
    for eventfile in eventfiles:
        with fits.open(eventfile, memmap=True) as f:
            hdu = f["EVENTS"]
            wcs_in = wcs_from_header(hdu.header)
            ra, dec = wcs_in.wcs_pix2world(hdu.data["X"], hdu.data["Y"], 1)
            xx, yy = wcs_out.wcs_world2pix(ra, dec, 1)
            x.append(xx)
            y.append(yy)
            e.append(hdu.data["ENERGY"])
            ch.append(hdu.data[hdu.header["CHANTYPE"]])
            t.append(hdu.data["TIME"])
            chantype = hdu.header["CHANTYPE"]

    x = np.concatenate(x)
    y = np.concatenate(y)
    e = np.concatenate(e)
    ch = np.concatenate(ch)
    t = np.concatenate(t)

    if chantype == "PHA":
        cunit = "adu"
    elif chantype == "PI":
        cunit = "Chan"

    col_x = fits.Column(name='X', format='D', unit='pixel', array=x)
    col_y = fits.Column(name='Y', format='D', unit='pixel', array=y)
    col_e = fits.Column(name='ENERGY', format='E', unit='eV', array=e)
    col_ch = fits.Column(name=chantype, format='1J', unit=cunit, array=ch)
    col_t = fits.Column(name='TIME', format='1D', unit='s', array=t)

    coldefs = fits.ColDefs([col_e, col_x, col_y, col_ch, col_t])
    tbhdu = fits.BinTableHDU.from_columns(coldefs)
    tbhdu.name = "EVENTS"

    tbhdu.header["HDUVERS"] = "1.1.0"
    tbhdu.header["RADECSYS"] = "FK5"
    tbhdu.header["EQUINOX"] = 2000.0
    tbhdu.header["HDUCLASS"] = "OGIP"
    tbhdu.header["HDUCLAS1"] = "EVENTS"
    tbhdu.header["HDUCLAS2"] = "ACCEPTED"
    tbhdu.header["MTYPE1"] = "sky"
    tbhdu.header["MFORM1"] = "x,y"
    tbhdu.header["MTYPE2"] = "EQPOS"
    tbhdu.header["MFORM2"] = "RA,DEC"
    tbhdu.header["TCTYP2"] = "RA---TAN"
    tbhdu.header["TCTYP3"] = "DEC--TAN"
    tbhdu.header["TCRVL2"] = wcs_out.wcs.crval[0]
    tbhdu.header["TCRVL3"] = wcs_out.wcs.crval[1]
    tbhdu.header["TCDLT2"] = wcs_out.wcs.cdelt[0]
    tbhdu.header["TCDLT3"] = wcs_out.wcs.cdelt[1]
    tbhdu.header["TCRPX2"] = wcs_out.wcs.crpix[0]
    tbhdu.header["TCRPX3"] = wcs_out.wcs.crpix[1]
    tbhdu.header["TCUNI2"] = "deg"
    tbhdu.header["TCUNI3"] = "deg"
    tbhdu.header["TLMIN2"] = 0.5
    tbhdu.header["TLMIN3"] = 0.5
    tbhdu.header["TLMAX2"] = shape_out[0]+0.5
    tbhdu.header["TLMAX3"] = shape_out[1]+0.5
    tbhdu.header["CHANTYPE"] = chantype

    with fits.open(eventfiles[0], memmap=True) as f:
        for key in ["TELESCOP", "INSTRUME", "ANCRFILE", "RESPFILE",
                    "EXPOSURE", "DATE", "DATE-OBS", "DATE-END", "MISSION",
                    "TSTART", "TSTOP", "TLMIN4", "TLMAX4", "PHA_BINS"]:
            tbhdu.header[key] = f["EVENTS"].header[key]

    start = fits.Column(name='START', format='1D', unit='s',
                          array=np.array([0.0]))
    stop = fits.Column(name='STOP', format='1D', unit='s',
                         array=np.array([tbhdu.header["EXPOSURE"]]))

    tbhdu_gti = fits.BinTableHDU.from_columns([start,stop])
    tbhdu_gti.name = "STDGTI"
    tbhdu_gti.header["TSTART"] = 0.0
    tbhdu_gti.header["TSTOP"] = tbhdu.header["EXPOSURE"]
    tbhdu_gti.header["HDUCLASS"] = "OGIP"
    tbhdu_gti.header["HDUCLAS1"] = "GTI"
    tbhdu_gti.header["HDUCLAS2"] = "STANDARD"
    tbhdu_gti.header["RADECSYS"] = "FK5"
    tbhdu_gti.header["EQUINOX"] = 2000.0
    tbhdu_gti.header["DATE"] = tbhdu.header["DATE"]
    tbhdu_gti.header["DATE-OBS"] = tbhdu.header["DATE-OBS"]
    tbhdu_gti.header["DATE-END"] = tbhdu.header["DATE-END"]

    hdulist = [fits.PrimaryHDU(), tbhdu, tbhdu_gti]

    fits.HDUList(hdulist).writeto(outfile, overwrite=overwrite)


def make_mosaic_image(evtfile_list, image_file, evt_file=None, emin=None, emax=None,
                      reblock=1, use_expmap=False, expmap_energy=None, 
                      expmap_weights=None, normalize=True, nhistx=16,
                      nhisty=16, overwrite=False):
    """
    Make a single FITS image from a grid of observations. Optionally,
    an exposure map can be computed and a flux image may be generated.

    Parameters
    ----------
    evtfile_list : filename
        The ASCII table produced by :meth:`~soxs.grid.observe_grid_source`
        containing the information about the event files and their
        locations on the sky.
    image_file : filename
        The name of the FITS image file to be written. This name will
        also be used for the exposure map, event, and flux files if they are
        written.
    evt_file : filename, optional
        The name of the mosaicked FITS event file to be written, if desired. Default
        is None, which writes no file.
    emin : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
        The minimum energy of the photons to put in the image, in keV.
    emax : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
        The maximum energy of the photons to put in the image, in keV.
    reblock : integer, optional
        Supply an integer power of 2 here to make an exposure map 
        with a different binning. Default: 1
    use_expmap : boolean, optional
        Whether or not to use (and potentially generate) an exposure map
        and a flux map. Default: False
    expmap_energy : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, or NumPy array, optional
        The energy in keV to use when computing the exposure map, or 
        a set of energies to be used with the *weights* parameter. If
        providing a set, it must be in keV.
    expmap_weights : array-like, optional
        The weights to use with a set of energies given in the
        *energy* parameter. Used to create a more accurate exposure
        map weighted by a range of energies. Default: None
    overwrite : boolean, optional
        Whether to overwrite an existing file with the same name.
        Default: False
    """
    try:
        from reproject.mosaicking import find_optimal_celestial_wcs, \
            reproject_and_coadd
        from reproject import reproject_interp
    except ImportError:
        raise ImportError("The mosaic functionality of SOXS requires the "
                          "'reproject' package to be installed!")
    t = ascii.read(evtfile_list, format='commented_header',
                   guess=False, header_start=0, delimiter="\t")

    files = []
    for row in t:
        evtfile = row["evtfile"]
        img_file = evtfile.replace("evt", "img")
        if use_expmap:
            emap_file = evtfile.replace("evt", "expmap")
            make_exposure_map(evtfile, emap_file, energy=expmap_energy,
                              weights=expmap_weights, normalize=normalize,
                              overwrite=overwrite, reblock=reblock, 
                              nhistx=nhistx, nhisty=nhisty)
        else:
            emap_file = None
        write_image(evtfile, img_file, emin=emin, emax=emax,
                    overwrite=overwrite, reblock=reblock)
        files.append([img_file, emap_file])

    img_hdus = [fits.open(fns[0], memmap=True)[0] for fns in files]
    wcs_out, shape_out = find_optimal_celestial_wcs(img_hdus)
    if evt_file is not None:
        _combine_events(t["evtfile"], wcs_out, (shape_out[1], shape_out[0]), evt_file,
                        overwrite=overwrite)
    img, footprint = reproject_and_coadd(img_hdus, wcs_out, shape_out=shape_out,
                                         reproject_function=reproject_interp,
                                         combine_function='sum')
    hdu = fits.PrimaryHDU(img, header=wcs_out.to_header())
    hdu.writeto(image_file, overwrite=overwrite)

    if use_expmap:
        if expmap_energy is None:
            raise RuntimeError("The 'expmap_energy' argument must be set if "
                               "making a mosaicked exposure map!")
        emap_hdus = [fits.open(fns[1], memmap=True)[1] for fns in files]
        emap, footprint = reproject_and_coadd(
            emap_hdus, wcs_out, shape_out=shape_out,
            reproject_function=reproject_interp, combine_function='sum')
        hdu = fits.PrimaryHDU(emap, header=wcs_out.to_header())
        expmap_file = image_file.replace("fits", "expmap")
        hdu.writeto(expmap_file, overwrite=overwrite)

        with np.errstate(invalid='ignore', divide='ignore'):
            flux = img / emap
        flux[np.isinf(flux)] = 0.0
        flux = np.nan_to_num(flux)
        flux[flux < 0.0] = 0.0
        hdu = fits.PrimaryHDU(flux, header=wcs_out.to_header())
        flux_file = image_file.replace("fits", "flux")
        hdu.writeto(flux_file, overwrite=overwrite)
