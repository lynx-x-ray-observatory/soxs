import astropy.io.fits as pyfits
import numpy as np
import astropy.wcs as pywcs

from xrs_tools.simput import read_simput_catalog
from xrs_tools.utils import mylog
from xrs_tools.instrument import instrument_registry, \
    AuxiliaryResponseFile, RedistributionMatrixFile
from xrs_tools.constants import erg_per_keV

def write_event_file(events, parameters, filename, clobber=False):
    from astropy.time import Time, TimeDelta

    t_begin = Time.now()
    dt = TimeDelta(parameters["exposure_time"], format='sec')
    t_end = t_begin + dt

    num_events = len(events["energy"])

    col_x = pyfits.Column(name='X', format='D', unit='pixel', array=events["xpix"])
    col_y = pyfits.Column(name='Y', format='D', unit='pixel', array=events["ypix"])
    col_e = pyfits.Column(name='ENERGY', format='E', unit='eV', array=events["energy"]*1000.)

    chantype = parameters["channel_type"]
    if chantype == "PHA":
        cunit = "adu"
    elif chantype == "PI":
        cunit = "Chan"
    col_ch = pyfits.Column(name=chantype.upper(), format='1J', unit=cunit, array=events[chantype])

    time = np.random.uniform(size=num_events, low=0.0, high=parameters["exposure_time"])
    col_t = pyfits.Column(name="TIME", format='1D', unit='s', array=time)

    cols = [col_e, col_x, col_y, col_ch, col_t]

    coldefs = pyfits.ColDefs(cols)
    tbhdu = pyfits.BinTableHDU.from_columns(coldefs)
    tbhdu.update_ext_name("EVENTS")

    tbhdu.header["MTYPE1"] = "sky"
    tbhdu.header["MFORM1"] = "x,y"
    tbhdu.header["MTYPE2"] = "EQPOS"
    tbhdu.header["MFORM2"] = "RA,DEC"
    tbhdu.header["TCTYP2"] = "RA---TAN"
    tbhdu.header["TCTYP3"] = "DEC--TAN"
    tbhdu.header["TCRVL2"] = parameters["sky_center"][0]
    tbhdu.header["TCRVL3"] = parameters["sky_center"][1]
    tbhdu.header["TCDLT2"] = -parameters["dtheta"]
    tbhdu.header["TCDLT3"] = parameters["dtheta"]
    tbhdu.header["TCRPX2"] = parameters["pix_center"][0]
    tbhdu.header["TCRPX3"] = parameters["pix_center"][1]
    tbhdu.header["TLMIN2"] = 0.5
    tbhdu.header["TLMIN3"] = 0.5
    tbhdu.header["TLMAX2"] = parameters["num_pixels"]+0.5
    tbhdu.header["TLMAX3"] = parameters["num_pixels"]+0.5
    tbhdu.header["TLMIN4"] = parameters["chan_lim"][0]
    tbhdu.header["TLMAX4"] = parameters["chan_lim"][1]
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
    tbhdu.header["RESPFILE"] = parameters["rmf"]
    tbhdu.header["PHA_BINS"] = parameters["nchan"]
    tbhdu.header["ANCRFILE"] = parameters["arf"]
    tbhdu.header["CHANTYPE"] = parameters["channel_type"]
    tbhdu.header["MISSION"] = parameters["mission"]
    tbhdu.header["TELESCOP"] = parameters["telescope"]
    tbhdu.header["INSTRUME"] = parameters["instrument"]

    start = pyfits.Column(name='START', format='1D', unit='s', array=np.array([0.0]))
    stop = pyfits.Column(name='STOP', format='1D', unit='s', array=np.array([parameters["exposure_time"]]))

    tbhdu_gti = pyfits.BinTableHDU.from_columns([start,stop])
    tbhdu_gti.update_ext_name("STDGTI")
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

def make_event_file(simput_file, out_file, exp_time, instrument,
                    sky_center, clobber=False, prng=np.random):
    """
    Take unconvolved events in a SIMPUT catalog and create an event
    file from them. This function does the following:

    1. Convolves the events with an ARF and RMF
    2. Pixelizes the events
    3. Writes the events to a file

    PSF effects and dithering are handled separately, in 

    Parameters
    ----------
    simput_file : string
        The SIMPUT catalog file to be used as input.
    out_file : string
        The name of the event file to be written.
    exp_time : float
        The exposure time to use, in seconds. 
    instrument : string
        The name of the instrument to use, which picks an instrument
        specification from the instrument registry.
    sky_center : array, tuple, or list
        The center RA, Dec coordinates of the observation, in degrees.
    clobber : boolean, optional
        Whether or not to clobber an existing file with the same name.
        Default: False
    prng : :class:`~numpy.random.RandomState` object or :mod:`~numpy.random`, optional
        A pseudo-random number generator. Typically will only be specified
        if you have a reason to generate the same set of random numbers, such as for a
        test. Default is the :mod:`numpy.random` module.

    Examples
    --------
    >>> make_event_file("sloshing_simput.fits", "sloshing_evt.fits", "hdxi",
    ...                 [30., 45.], clobber=True)
    """
    event_list, parameters = read_simput_catalog(simput_file)

    try:
        instrument_spec = instrument_registry[instrument]
    except KeyError:
        raise KeyError("Instrument %s is not in the instrument registry!" % instrument)
    arf_file = instrument_spec["arf"]
    rmf_file = instrument_spec["rmf"]
    arf = AuxiliaryResponseFile(arf_file)
    rmf = RedistributionMatrixFile(rmf_file)

    nx = instrument_spec["num_pixels"]
    dtheta = instrument_spec["dtheta"]/3600. # deg to arcsec

    event_params = {}
    event_params["exposure_time"] = exp_time
    event_params["arf"] = arf.filename
    event_params["sky_center"] = sky_center
    event_params["pix_center"] = np.array([0.5*(nx+1)]*2)
    event_params["num_pixels"] = nx
    event_params["dtheta"] = dtheta
    event_params["rmf"] = rmf.filename
    event_params["channel_type"] = rmf.header["CHANTYPE"]
    event_params["telescope"] = rmf.header["TELESCOP"]
    event_params["instrument"] = rmf.header["INSTRUME"]
    event_params["mission"] = rmf.header.get("MISSION", "")
    event_params["nchan"] = rmf.ebounds_header["DETCHANS"]
    num = 0
    for i in range(1, rmf.num_mat_columns+1):
        if rmf.header["TTYPE%d" % i] == "F_CHAN":
            num = i
            break
    event_params["chan_lim"] = [rmf.header["TLMIN%d" % num], 
                                rmf.header["TLMAX%d" % num]]

    w = pywcs.WCS(naxis=2)
    w.wcs.crval = event_params["sky_center"]
    w.wcs.crpix = event_params["pix_center"]
    w.wcs.cdelt = [-dtheta, dtheta]
    w.wcs.ctype = ["RA---TAN","DEC--TAN"]
    w.wcs.cunit = ["deg"]*2

    all_events = {}

    for i, events in enumerate(event_list):

        mylog.info("Creating events for source %d" % (i+1))

        # Step 1: Use ARF to determine which photons are observed

        mylog.info("Applying energy-dependent effective area from %s. This may take a minute." % arf_file)
        refband = [parameters["emin"][i], parameters["emax"][i]]
        events = arf.detect_events(events, exp_time, parameters["flux"][i],
                                   refband, prng=prng)
        if events["energy"].size == 0:
            mylog.warning("No events were observed for this source!!!")

        # Step 2: Assign pixel coordinates to events and clip events that don't fall
        # within the detection region

        if events["energy"].size > 0:
            mylog.info("Pixeling events.")

            xpix, ypix = w.wcs_world2pix(events["ra"], events["dec"], 1)

            events["xpix"] = xpix
            events["ypix"] = ypix

            keepx = np.logical_and(events["xpix"] >= 0.5, events["xpix"] <= nx+0.5)
            keepy = np.logical_and(events["ypix"] >= 0.5, events["ypix"] <= nx+0.5)
            keep = np.logical_and(keepx, keepy)
            if keep.sum() == 0:
                mylog.warning("No events are within the field of view for this source!!!")

            for key in events:
                events[key] = events[key][keep]

        # Step 3: Scatter energies with RMF

        if events["energy"].size > 0:
            mylog.info("Scattering energies with RMF.")
            events = rmf.scatter_energies(events, prng=prng)

        for key in events:
            if i == 0:
                all_events[key] = events[key]
            else:
                all_events[key] = np.append(all_events[key], events[key])

    if all_events["energy"].size == 0:
        raise RuntimeError("No events were detected!!!")

    write_event_file(all_events, event_params, out_file, clobber=clobber)

def add_background_events(bkgnd_spectrum, event_file, flat_response=False, 
                          prng=np.random):
    """
    Add background events to an event file. This is for astrophysical
    backgrounds as well as instrumental backgrounds. For the latter, 
    setting ``flat_response=True`` may be appropriate.

    Parameters
    ----------
    bkgnd_spectrum : :class:`~xrs_tools.spectra.Spectrum`
        A ``Spectrum`` object for the background.
    event_file : string
        The event file to add the new events to. 
    flat_response : boolean, optional
        If True, a flat ARF is assumed, which may be appropriate
        for instrumental backgrounds. Default: False
    prng : :class:`~numpy.random.RandomState` object or :mod:`~numpy.random`, optional
        A pseudo-random number generator. Typically will only be specified
        if you have a reason to generate the same set of random numbers, such as for a
        test. Default is the :mod:`numpy.random` module.

    """
    f = pyfits.open(event_file, mode='update')
    xmax = f["EVENTS"].header["TLMAX2"]
    ymax = f["EVENTS"].header["TLMAX3"]
    exp_time = f["EVENTS"].header["EXPOSURE"]
    if flat_response:
        area = 0.0
    else:
        arf = AuxiliaryResponseFile(f["EVENTS"].header["ANCRFILE"])
        area = arf.max_area
    bkg_events = {}
    bkg_events["energy"] = bkgnd_spectrum.generate_energies(exp_time, 
                                                            area, prng=prng)
    flux = bkg_events["energy"].sum()*erg_per_keV/exp_time/area
    if not flat_response:
        refband = [bkg_events["energy"].min(), bkg_events["energy"].max()]
        arf.detect_events(bkg_events, exp_time, flux, refband, prng=prng)
    n_events = bkg_events["energy"].size
    bkg_events['x'] = prng.uniform(low=0.5, high=xmax, size=n_events)
    bkg_events['y'] = prng.uniform(low=0.5, high=ymax, size=n_events)
    rmf = RedistributionMatrixFile(f["EVENTS"].header["RESPFILE"])
    bkg_events = rmf.scatter_energies(bkg_events, prng=prng)
    chantype = rmf.header["CHANTYPE"].upper()
    f["EVENTS"].data[chantype] = np.concatenate(f["EVENTS"].data[chantype][:],
                                                bkg_events[rmf.header["CHANTYPE"]])
    f["EVENTS"].data["X"] = np.concatenate(f["EVENTS"].data["X"][:], 
                                           bkg_events["x"])
    f["EVENTS"].data["Y"] = np.concatenate(f["EVENTS"].data["Y"][:],
                                           bkg_events["y"])
    f["EVENTS"].data["ENERGY"] = np.concatenate(f["EVENTS"].data["ENERGY"][:], 
                                                bkg_events["energy"])
    f["EVENTS"].data["Y"] = np.concatenate(f["EVENTS"].data["Y"][:],
                                           bkg_events["y"])
    t = np.random.uniform(size=n_events, low=0.0, high=exp_time)
    f["EVENTS"].data["TIME"] = np.concatenate(f["EVENTS"].data["TIME"][:], t)
    f.flush()
    f.close()
