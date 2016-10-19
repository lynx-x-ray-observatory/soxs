import astropy.io.fits as pyfits
import numpy as np
import astropy.wcs as pywcs
import os

from xrs_tools.simput import read_simput_catalog
from xrs_tools.utils import mylog, check_file_location
from xrs_tools.instrument import instrument_registry, \
    AuxiliaryResponseFile, RedistributionMatrixFile, \
    add_instrument_to_registry
from xrs_tools.constants import erg_per_keV

sigma_to_fwhm = 2.*np.sqrt(2.*np.log(2.))

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
    tbhdu.header["RESPFILE"] = parameters["rmf"]
    tbhdu.header["PHA_BINS"] = parameters["nchan"]
    tbhdu.header["ANCRFILE"] = parameters["arf"]
    tbhdu.header["CHANTYPE"] = parameters["channel_type"]
    tbhdu.header["MISSION"] = parameters["mission"]
    tbhdu.header["TELESCOP"] = parameters["telescope"]
    tbhdu.header["INSTRUME"] = parameters["instrument"]
    tbhdu.header["RA_PNT"] = parameters["sky_center"][0]
    tbhdu.header["DEC_PNT"] = parameters["sky_center"][1]
    tbhdu.header["ROLL_PNT"] = parameters["roll_angle"]

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
                    sky_center, clobber=False, dither_shape="square",
                    dither_size=16.0, roll_angle=0.0, prng=np.random):
    """
    Take unconvolved events in a SIMPUT catalog and create an event
    file from them. This function does the following:

    1. Convolves the events with an ARF and RMF
    2. Pixelizes the events, applying PSF effects and dithering
    3. Writes the events to a file

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
        specification from the instrument registry. Can also be a JSON
        file with a new instrument specification. If this is the case,
        it will be loaded into the instrument registry. 
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

    if instrument not in instrument_registry and os.path.exists(instrument):
        instrument = add_instrument_to_registry(instrument)
    try:
        instrument_spec = instrument_registry[instrument]
    except KeyError:
        raise KeyError("Instrument %s is not in the instrument registry!" % instrument)
    arf_file = check_file_location(instrument_spec["arf"], "files")
    rmf_file = check_file_location(instrument_spec["rmf"], "files")
    arf = AuxiliaryResponseFile(arf_file)
    rmf = RedistributionMatrixFile(rmf_file)

    nx = instrument_spec["num_pixels"]
    dtheta = instrument_spec["dtheta"]/3600. # deg to arcsec

    event_params = {}
    event_params["exposure_time"] = exp_time
    event_params["arf"] = os.path.split(arf.filename)[-1]
    event_params["sky_center"] = sky_center
    event_params["pix_center"] = np.array([0.5*(nx+1)]*2)
    event_params["num_pixels"] = nx
    event_params["dtheta"] = dtheta
    event_params["rmf"] = os.path.split(rmf.filename)[-1]
    event_params["channel_type"] = rmf.header["CHANTYPE"]
    event_params["telescope"] = rmf.header["TELESCOP"]
    event_params["instrument"] = rmf.header["INSTRUME"]
    event_params["mission"] = rmf.header.get("MISSION", "")
    event_params["nchan"] = rmf.ebounds_header["DETCHANS"]
    event_params["roll_angle"] = roll_angle
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

        mylog.info("Applying energy-dependent effective area from %s. " % event_params["arf"] +
                   "This may take a minute.")
        refband = [parameters["emin"][i], parameters["emax"][i]]
        events = arf.detect_events(events, exp_time, parameters["flux"][i],
                                   refband, prng=prng)
        if events["energy"].size == 0:
            mylog.warning("No events were observed for this source!!!")

        # Step 2: Assign pixel coordinates to events. Apply dithering and
        # PSF. Clip events that don't fall within the detection region.

        if events["energy"].size > 0:
            mylog.info("Pixeling events.")

            # Convert RA, Dec to pixel coordinates
            xpix, ypix = w.wcs_world2pix(events["ra"], events["dec"], 1)

            xpix -= event_params["pix_center"][0]
            ypix -= event_params["pix_center"][1]

            events.pop("ra")
            events.pop("dec")

            n_evt = xpix.size

            # Dither pixel coordinates

            x_offset = 0.0
            y_offset = 0.0
            if dither_shape == "circle":
                r = dither_size*prng.uniform(size=n_evt)
                theta = 2.*np.pi*prng.uniform(size=n_evt)
                x_offset = r*np.cos(theta)
                y_offset = r*np.sin(theta)
            elif dither_shape == "square":
                x_offset = dither_size*prng.uniform(low=-0.5, high=0.5, size=n_evt)
                y_offset = dither_size*prng.uniform(low=-0.5, high=0.5, size=n_evt)

            xpix -= x_offset
            ypix -= y_offset

            roll_angle = np.deg2rad(roll_angle)

            # Rotate physical coordinates to detector coordinates

            rot_mat = np.array([[np.sin(roll_angle), -np.cos(roll_angle)],
                                [-np.cos(roll_angle), -np.sin(roll_angle)]])

            det = np.dot(rot_mat, np.array([xpix, ypix]))
            detx = det[0,:]
            dety = det[1,:]

            # PSF scattering of detector coordinates

            sigma = instrument_spec["psf_fwhm"]/sigma_to_fwhm/instrument_spec["dtheta"]

            detx += prng.normal(loc=0.0, scale=sigma, size=n_evt)
            dety += prng.normal(loc=0.0, scale=sigma, size=n_evt)

            # Convert detector coordinates to chip coordinates

            events["chipx"] = np.round(detx + event_params['pix_center'][0])
            events["chipy"] = np.round(dety + event_params['pix_center'][1])

            # Throw out events that don't fall on the chip

            keepx = np.logical_and(events["chipx"] >= 1.0, events["chipx"] <= nx)
            keepy = np.logical_and(events["chipy"] >= 1.0, events["chipy"] <= nx)
            keep = np.logical_and(keepx, keepy)
            if keep.sum() == 0:
                mylog.warning("No events are within the field of view for this source!!!")

            for key in events:
                events[key] = events[key][keep]

            n_evt = events["energy"].shape

            # Convert chip coordinates back to detector coordinates

            events["detx"] = np.round(events["chipx"] - event_params['pix_center'][0] +
                                      prng.uniform(low=-0.5, high=0.5, size=n_evt))
            events["dety"] = np.round(events["chipy"] - event_params['pix_center'][1] +
                                      prng.uniform(low=-0.5, high=0.5, size=n_evt))

            # Convert detector coordinates back to pixel coordinates

            pix = np.dot(rot_mat, np.array([events["detx"], events["dety"]]))

            events["xpix"] = pix[0,:] + event_params['pix_center'][0] + x_offset
            events["ypix"] = pix[1,:] + event_params['pix_center'][1] + y_offset

    # Step 3: Scatter energies with RMF

        if events["energy"].size > 0:
            mylog.info("Scattering energies with RMF %s." % event_params['rmf'])
            events = rmf.scatter_energies(events, prng=prng)

        for key in events:
            if i == 0:
                all_events[key] = events[key]
            else:
                all_events[key] = np.append(all_events[key], events[key])

    if all_events["energy"].size == 0:
        raise RuntimeError("No events were detected!!!")

    all_events['time'] = np.random.uniform(size=all_events["energy"].size, low=0.0, 
                                           high=parameters["exposure_time"])

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
    cxmax = f["EVENTS"].header["TLMAX2"]-0.5
    cymax = f["EVENTS"].header["TLMAX3"]-0.5
    xc = f["EVENTS"].header["TCRPX2"]
    yc = f["EVENTS"].header["TCRPX3"]
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
    bkg_events['chipx'] = np.round(prng.uniform(low=1.0, high=cxmax, size=n_events))
    bkg_events['chipy'] = np.round(prng.uniform(low=1.0, high=cymax, size=n_events))
    roll_angle = np.deg2rad(f["EVENTS"].header["ROLL_PNT"])
    rot_mat = np.array([[np.sin(roll_angle), -np.cos(roll_angle)],
                        [-np.cos(roll_angle), -np.sin(roll_angle)]])
    bkg_events["detx"] = np.round(bkg_events["chipx"] - xc +
                                  prng.uniform(low=-0.5, high=0.5, size=n_events))
    bkg_events["dety"] = np.round(bkg_events["chipy"] - yc +
                                  prng.uniform(low=-0.5, high=0.5, size=n_events))
    pix = np.dot(rot_mat, np.array([bkg_events["detx"], bkg_events["dety"]]))
    bkg_events["xpix"] = pix[0,:]
    bkg_events["ypix"] = pix[1,:]
    rmf = RedistributionMatrixFile(f["EVENTS"].header["RESPFILE"])
    bkg_events = rmf.scatter_energies(bkg_events, prng=prng)
    chantype = rmf.header["CHANTYPE"].upper()
    f["EVENTS"].data[chantype] = np.concatenate(f["EVENTS"].data[chantype][:],
                                                bkg_events[rmf.header["CHANTYPE"]])
    f["EVENTS"].data["X"] = np.concatenate(f["EVENTS"].data["X"][:], 
                                           bkg_events["xpix"])
    f["EVENTS"].data["Y"] = np.concatenate(f["EVENTS"].data["Y"][:],
                                           bkg_events["ypix"])
    f["EVENTS"].data["ENERGY"] = np.concatenate(f["EVENTS"].data["ENERGY"][:], 
                                                bkg_events["energy"])
    f["EVENTS"].data[chantype] = np.concatenate(f["EVENTS"].data[chantype][:],
                                                bkg_events[chantype])
    t = np.random.uniform(size=n_events, low=0.0, high=exp_time)
    f["EVENTS"].data["TIME"] = np.concatenate(f["EVENTS"].data["TIME"][:], t)
    f.flush()
    f.close()
