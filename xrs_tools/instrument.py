from __future__ import print_function

import numpy as np
import astropy.wcs as pywcs
import astropy.io.fits as pyfits
from astropy.utils.console import ProgressBar

from xrs_tools.simput import read_simput_phlist
from xrs_tools.utils import mylog, iterable, ensure_numpy_array

sigma_to_fwhm = 2.*np.sqrt(2.*np.log(2.))

class AuxiliaryResponseFile(object):
    r"""
    A class for auxiliary response files (ARFs).

    Parameters
    ----------
    filename : string
        The filename of the ARF to be read.

    Examples
    --------
    >>> arf = AuxiliaryResponseFile("xrs_calorimeter.arf")
    """
    def __init__(self, filename):
        self.filename = filename
        f = pyfits.open(self.filename)
        self.elo = f["SPECRESP"].data.field("ENERG_LO"), "keV"
        self.ehi = f["SPECRESP"].data.field("ENERG_HI")
        self.emid = 0.5*(self.elo+self.ehi)
        self.eff_area = np.nan_to_num(f["SPECRESP"].data.field("SPECRESP"))
        f.close()

    def __str__(self):
        return self.filename

    def detect_events(self, energy, area, prng=None):
        """
        Use the ARF to determine a subset of photons which will be
        detected. Returns a boolean NumPy array which is the same
        is the same size as the number of photons, wherever it is
        "true" means those photons have been detected.

        Parameters
        ----------
        energy : np.ndarray
            The energies of the photons to attempt to detect, in keV.
        area : float
            The collecting area in cm^2 associated with the event energies.
        prng : :class:`~numpy.random.RandomState` object or :mod:`~numpy.random`, optional
            A pseudo-random number generator. Typically will only be specified
            if you have a reason to generate the same set of random numbers, such as for a
            test. Default is the :mod:`~numpy.random` module.
        """
        if prng is None:
            prng = np.random
        earea = np.interp(energy, self.emid, self.eff_area, left=0.0, right=0.0)
        randvec = area*prng.uniform(size=energy.shape)
        return randvec < earea

class RedistributionMatrixFile(object):
    r"""
    A class for redistribution matrix files (RMFs).

    Parameters
    ----------
    filename : string
        The filename of the RMF to be read.

    Examples
    --------
    >>> rmf = RedistributionMatrixFile("xrs_hdxi.rmf")
    """
    def __init__(self, filename):
        self.filename = filename
        self.handle = pyfits.open(self.filename)
        if "MATRIX" in self.handle:
            self.mat_key = "MATRIX"
        elif "SPECRESP MATRIX" in self.handle:
            self.mat_key = "SPECRESP MATRIX"
        else:
            raise RuntimeError("Cannot find the response matrix in the RMF "
                               "file %s! " % filename+"It should be named "
                                                      "\"MATRIX\" or \"SPECRESP MATRIX\".")
        self.data = self.handle[self.mat_key].data
        self.header = self.handle[self.mat_key].header
        self.num_mat_columns = len(self.handle[self.mat_key].columns)
        self.ebounds = self.handle["EBOUNDS"].data
        self.ebounds_header = self.handle["EBOUNDS"].header
        self.weights = np.array([w.sum() for w in self.data["MATRIX"]])

    def __str__(self):
        return self.filename

def write_event_file(events, parameters, filename, clobber=False):
    from astropy.time import Time, TimeDelta

    t_begin = Time.now()
    dt = TimeDelta(parameters["exposure_time"], format='sec')
    t_end = t_begin + dt

    num_events = len(events["energy"])

    col_e = pyfits.Column(name='ENERGY', format='E', unit='eV', array=events["energy"]*1000.)
    col_x = pyfits.Column(name='X', format='D', unit='pixel', array=events["xpix"])
    col_y = pyfits.Column(name='Y', format='D', unit='pixel', array=events["ypix"])

    chantype = parameters["channel_type"]
    if chantype == "PHA":
        cunit = "adu"
    elif chantype == "PI":
        cunit = "Chan"
    col_ch = pyfits.Column(name=chantype.upper(), format='1J', unit=cunit, array=events[chantype])

    time = np.random.uniform(size=num_events, low=0.0, high=parameters["exposure_time"])
    col_t = pyfits.Column(name="TIME", format='1D', unit='s', array=time)

    cols = [col_x, col_y, col_e, col_ch, col_t]

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

instrument_registry = {}
instrument_registry["xcal"] = {"arf": "xrs_calorimeter.arf",
                               "rmf": "xrs_calorimeter.rmf",
                               "num_pixel": 300,
                               "dtheta": 1.0}
instrument_registry["hdxi"] = {"arf": "xrs_hdxi.arf",
                               "rmf": "xrs_hdxi.rmf",
                               "num_pixel": 4096,
                               "dtheta": 1./3.}

def add_instrument_to_registry(filename):
    """
    Add an instrument specification to the registry contained
    in a JSON file. The JSON file must have this structure (the
    order is not important):

    {'name': 'hdxi', # The short name of the instrument
     'arf': 'xrs_hdxi.arf', # The file containing the ARF
     'rmf': 'xrs_hdxi.rmf' # The file containing the RMF
     'dtheta': 0.33333333333, # The central pixel scale in arcsec
     'num_pixel': 4096} # The number of pixels on a side in the FOV

    Parameters
    ----------
    filename : string
        The JSON file containing the instrument specification.
    """
    import json
    f = open(filename)
    inst = json.load(f)
    f.close()
    name = inst.pop("name")
    instrument_registry[name] = inst
    return name

def show_instrument_registry():
    """
    Print the contents of the instrument registry.
    """
    for name, spec in instrument_registry.items():
        print("Instrument: %s" % name.upper())
        for k, v in spec.items():
            print("    %s: %s" % (k, v))

def make_event_file(simput_file, out_file, exp_time, instrument, 
                    sky_center, clobber=False, prng=np.random):
    """
    Take unconvolved events in a SIMPUT file and create an event
    file from them. This function does the following:

    1. Convolves the events with an ARF and RMF
    2. Pixelizes the events
    3. Writes the events to a file

    PSF effects and dithering are handled separately, in 

    Parameters
    ----------
    simput_file : string
        The SIMPUT file to be used as input.
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
    events, parameters = read_simput_phlist(simput_file)

    try:
        instrument_spec = instrument_registry[instrument]
    except KeyError:
        raise KeyError("Instrument %s is not in the instrument registry!" % instrument)
    arf_file = instrument_spec["arf"]
    rmf_file = instrument_spec["rmf"]
    nx = instrument_spec["num_pixels"]
    dtheta = instrument_spec["dtheta"]/3600. # deg to arcsec

    if exp_time > parameters["exposure_time"]:
        raise ValueError("Specified exposure time %g s cannot be larger " % exp_time +
                         "than maximum exposure time %s!" % parameters["exposure_time"])
    area = parameters["exposure_time"]/parameters["flux"]
    parameters["exposure_time"] = exp_time

    # Step 1: Use ARF to determine which photons are observed

    mylog.info("Applying energy-dependent effective area from %s." % arf_file)
    arf = AuxiliaryResponseFile(arf_file)
    detected = arf.detect_events(events["energy"], area, prng=prng)
    mylog.info("%s events detected." % detected.sum())
    for key in events:
        events[key] = events[key][detected]
    parameters["arf"] = arf.filename

    # Step 2: Assign pixel coordinates to events and clip events that don't fall
    # within the detection region

    mylog.info("Pixeling image.")

    parameters["sky_center"] = sky_center
    parameters["pix_center"] = np.array([0.5*(nx+1)]*2)
    parameters["num_pixels"] = nx
    parameters["dtheta"] = dtheta

    w = pywcs.WCS(naxis=2)
    w.wcs.crval = parameters["sky_center"]
    w.wcs.crpix = parameters["pix_center"]
    w.wcs.cdelt = [-dtheta, dtheta]
    w.wcs.ctype = ["RA---TAN","DEC--TAN"]
    w.wcs.cunit = ["deg"]*2

    events["xpix"], events["ypix"] = w.wcs_world2pix(events["xsky"], events["ysky"], 1)

    keepx = np.logical_and(events["xpix"] < 0.5, events["xpix"] > nx+0.5)
    keepy = np.logical_and(events["ypix"] < 0.5, events["ypix"] > nx+0.5)
    keep = np.logical_and(keepx, keepy)
    for key in events:
        events[key] = events[key][keep]

    # Step 3: Scatter energies with RMF

    mylog.info("Reading response matrix file (RMF): %s" % rmf_file)
    rmf = RedistributionMatrixFile(rmf_file)

    elo = rmf.data["ENERG_LO"]
    ehi = rmf.data["ENERG_HI"]
    n_de = elo.shape[0]
    mylog.info("Number of energy bins in RMF: %d" % n_de)
    mylog.info("Energy limits: %g %g" % (min(elo), max(ehi)))

    n_ch = len(rmf.ebounds["CHANNEL"])
    mylog.info("Number of channels in RMF: %d" % n_ch)

    eidxs = np.argsort(events["energy"])
    sorted_e = events["energy"][eidxs]

    detectedChannels = []

    # run through all photon energies and find which bin they go in
    fcurr = 0
    last = sorted_e.shape[0]

    mylog.info("Scattering energies with RMF.")

    with ProgressBar(last) as pbar:
        for (k, low), high in zip(enumerate(elo), ehi):
            # weight function for probabilities from RMF
            weights = np.nan_to_num(np.float64(rmf.data["MATRIX"][k]))
            weights /= weights.sum()
            # build channel number list associated to array value,
            # there are groups of channels in rmfs with nonzero probabilities
            trueChannel = []
            f_chan = ensure_numpy_array(np.nan_to_num(rmf.data["F_CHAN"][k]))
            n_chan = ensure_numpy_array(np.nan_to_num(rmf.data["N_CHAN"][k]))
            if not iterable(f_chan):
                f_chan = [f_chan]
                n_chan = [n_chan]
            for start, nchan in zip(f_chan, n_chan):
                if nchan == 0:
                    trueChannel.append(start)
                else:
                    trueChannel += list(range(start, start+nchan))
            if len(trueChannel) > 0:
                for q in range(fcurr, last):
                    if low <= sorted_e[q] < high:
                        channelInd = prng.choice(len(weights), p=weights)
                        fcurr += 1
                        pbar.update(fcurr)
                        detectedChannels.append(trueChannel[channelInd])
                    else:
                        break

    for key in events:
        events[key] = events[key][eidxs]

    events[rmf.header["CHANTYPE"]] = np.array(detectedChannels, dtype="int")

    parameters["rmf"] = rmf.filename
    parameters["channel_type"] = rmf.header["CHANTYPE"]
    parameters["telescope"] = rmf.header["TELESCOP"]
    parameters["instrument"] = rmf.header["INSTRUME"]
    parameters["mission"] = rmf.header["MISSION"]
    parameters["nchan"] = rmf.ebounds_header["DETCHANS"]
    num = 0
    for i in range(1, rmf.num_mat_columns+1):
        if rmf.header["TTYPE%d" % i] == "F_CHAN":
            num = i
            break
    parameters["chan_lim"] = [rmf.header["TLMIN%d" % num], rmf.header["TLMAX%d" % num]]

    write_event_file(events, parameters, out_file, clobber=clobber)
