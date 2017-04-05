import numpy as np
import astropy.io.fits as pyfits
import astropy.wcs as pywcs
import astropy.units as u
import os
from collections import defaultdict

from soxs.constants import erg_per_keV
from soxs.simput import read_simput_catalog
from soxs.utils import mylog, check_file_location, \
    ensure_numpy_array, parse_prng
from soxs.events import write_event_file
from soxs.background import make_instrument_background, \
    make_foreground, add_background_from_file, \
    make_ptsrc_background
from soxs.instrument_registry import instrument_registry
from six import string_types
from tqdm import tqdm

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
    >>> arf = AuxiliaryResponseFile("xrs_mucal_3x10.arf")
    """
    def __init__(self, filename):
        self.filename = check_file_location(filename, "files")
        f = pyfits.open(self.filename)
        self.elo = f["SPECRESP"].data.field("ENERG_LO")
        self.ehi = f["SPECRESP"].data.field("ENERG_HI")
        self.emid = 0.5*(self.elo+self.ehi)
        self.eff_area = np.nan_to_num(f["SPECRESP"].data.field("SPECRESP")).astype("float64")
        self.max_area = self.eff_area.max()
        f.close()

    def __str__(self):
        return self.filename

    def interpolate_area(self, energy):
        """
        Interpolate the effective area to the energies 
        provided  by the supplied *energy* array.
        """
        earea = np.interp(energy, self.emid, self.eff_area, left=0.0, right=0.0)
        return u.Quantity(earea, "cm**2")

    def detect_events(self, events, exp_time, flux, refband, prng=None):
        """
        Use the ARF to determine a subset of photons which 
        will be detected. Returns a boolean NumPy array 
        which is the same is the same size as the number 
        of photons, wherever it is "true" means those photons 
        have been detected.

        Parameters
        ----------
        events : dict of np.ndarrays
            The energies and positions of the photons. 
        exp_time : float
            The exposure time in seconds.
        flux : float
            The total flux of the photons in erg/s/cm^2. 
        refband : array_like
            A two-element array or list containing the limits 
            of the energy band which the flux was computed in. 
        prng : :class:`~numpy.random.RandomState` object, integer, or None
            A pseudo-random number generator. Typically will only 
            be specified if you have a reason to generate the same 
            set of random numbers, such as for a test. Default is None, 
            which sets the seed based on the system time. 
        """
        prng = parse_prng(prng)
        energy = events["energy"]
        if energy.size == 0:
            return events
        earea = self.interpolate_area(energy).value
        idxs = np.logical_and(energy >= refband[0], energy <= refband[1])
        rate = flux/(energy[idxs].sum()*erg_per_keV)*earea[idxs].sum()
        n_ph = np.modf(rate*exp_time)
        n_ph = np.int64(n_ph[1]) + np.int64(n_ph[0] >= prng.uniform())
        fak = float(n_ph)/energy.size
        if fak > 1.0:
            mylog.error("Number of events in sample: %d, Number of events wanted: %d" % (energy.size, n_ph))
            raise ValueError("This combination of exposure time and effective area "
                             "will result in more photons being drawn than are available "
                             "in the sample!!!")
        w = earea / self.max_area
        randvec = prng.uniform(size=energy.size)
        eidxs = prng.permutation(np.where(randvec < w)[0])[:n_ph].astype("int64")
        mylog.info("%s events detected." % n_ph)
        for key in events:
            events[key] = events[key][eidxs]
        return events

class FlatResponse(AuxiliaryResponseFile):
    """
    A flat effective area response. 

    Parameters
    ----------
    emin : float
        The minimum energy of the response in keV.
    emax : float
        The maximum energy of the response in keV. 
    area : float
        The effective area in cm**2. 
    nbins : integer
        The number of bins in the response file. 

    Examples
    --------
    >>> arf = FlatResponse(0.1, 10.0, 3000.0, 10000)
    """
    def __init__(self, emin, emax, area, nbins):
        self.filename = "flat_response"
        de = (emax-emin)/nbins
        self.elo = np.arange(nbins)*de + emin
        self.ehi = self.elo + de
        self.emid = 0.5*(self.elo+self.ehi)
        self.eff_area = area*np.ones(nbins)
        self.max_area = area

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
        self.filename = check_file_location(filename, "files")
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
        self.elo = self.data["ENERG_LO"]
        self.ehi = self.data["ENERG_HI"]
        self.emid = 0.5*(self.elo+self.ehi)
        self.de = self.ehi-self.elo
        self.n_de = self.elo.size
        self.n_ch = len(self.ebounds["CHANNEL"])
        num = 0
        for i in range(1, self.num_mat_columns+1):
            if self.header["TTYPE%d" % i] == "F_CHAN":
                num = i
                break
        self.cmin = self.header["TLMIN%d" % num]
        self.cmax = self.header["TLMAX%d" % num]

    def __str__(self):
        return self.filename

    def scatter_energies(self, events, prng=None):
        """
        Scatter photon energies with the RMF and produce the 
        corresponding channel values.

        Parameters
        ----------
        events : dict of np.ndarrays
            The energies and positions of the photons. 
        prng : :class:`~numpy.random.RandomState` object, integer, or None
            A pseudo-random number generator. Typically will only 
            be specified if you have a reason to generate the same 
            set of random numbers, such as for a test. Default is None, 
            which sets the seed based on the system time. 
        """
        prng = parse_prng(prng)
        eidxs = np.argsort(events["energy"])
        sorted_e = events["energy"][eidxs]

        detectedChannels = []

        # run through all photon energies and find which bin they go in
        fcurr = 0
        last = sorted_e.shape[0]

        pbar = tqdm(leave=True, total=last, desc="Scattering energies ")
        for (k, low), high in zip(enumerate(self.elo), self.ehi):
            # weight function for probabilities from RMF
            weights = np.nan_to_num(np.float64(self.data["MATRIX"][k]))
            weights /= weights.sum()
            # build channel number list associated to array value,
            # there are groups of channels in rmfs with nonzero probabilities
            trueChannel = []
            f_chan = ensure_numpy_array(np.nan_to_num(self.data["F_CHAN"][k]))
            n_chan = ensure_numpy_array(np.nan_to_num(self.data["N_CHAN"][k]))
            for start, nchan in zip(f_chan, n_chan):
                if nchan == 0:
                    trueChannel.append(start)
                else:
                    trueChannel += list(range(start, start+nchan))
            trueChannel = np.array(trueChannel)
            if len(trueChannel) > 0:
                e = sorted_e[fcurr:last]
                nn = np.logical_and(low <= e, e < high).sum()
                channelInd = prng.choice(len(weights), size=nn, p=weights)
                detectedChannels.append(trueChannel[channelInd])
                fcurr += nn
                pbar.update(nn)

        pbar.close()

        for key in events:
            events[key] = events[key][eidxs]
        events[self.header["CHANTYPE"]] = np.concatenate(detectedChannels)

        return events


def generate_events(input_events, exp_time, instrument, sky_center, 
                    dither_shape="square", dither_size=16.0, roll_angle=0.0, 
                    prng=None):
    """
    Take unconvolved events and convolve them with instrumental responses. This 
    function does the following:

    1. Determines which events are observed using the ARF
    2. Pixelizes the events, applying PSF effects and dithering
    3. Determines energy channels using the RMF

    This function is not meant to be called by the end-user but is used by
    the :func:`~soxs.instrument.instrument_simulator` function.

    Parameters
    ----------
    input_events : string, dict, or None
        The unconvolved events to be used as input. Can be one of the
        following:
        1. The name of a SIMPUT catalog file.
        2. A Python dictionary containing the following items:
        "ra": A NumPy array of right ascension values in degrees.
        "dec": A NumPy array of declination values in degrees.
        "energy": A NumPy array of energy values in keV.
        "flux": The flux of the entire source, in units of erg/cm**2/s.
    out_file : string
        The name of the event file to be written.
    exp_time : float
        The exposure time to use, in seconds. 
    instrument : string
        The name of the instrument to use, which picks an instrument
        specification from the instrument registry. 
    sky_center : array, tuple, or list
        The center RA, Dec coordinates of the observation, in degrees.
    dither_shape : string
        The shape of the dither. Currently "circle" or "square" 
        Default: "square"
    dither_size : float
        The size of the dither in arcseconds. Width of square or radius
        of circle. Default: 16.0
    roll_angle : float
        The roll angle of the observation in degrees. Default: 0.0
    prng : :class:`~numpy.random.RandomState` object, integer, or None
        A pseudo-random number generator. Typically will only 
        be specified if you have a reason to generate the same 
        set of random numbers, such as for a test. Default is None, 
        which sets the seed based on the system time. 
    """
    prng = parse_prng(prng)
    if isinstance(input_events, dict):
        parameters = {}
        for key in ["flux", "emin", "emax", "sources"]:
            parameters[key] = input_events[key]
        event_list = []
        for i in range(len(parameters["flux"])):
            edict = {}
            for key in ["ra", "dec", "energy"]:
                edict[key] = input_events[key][i]
            event_list.append(edict)
    elif isinstance(input_events, string_types):
        # Assume this is a SIMPUT catalog
        event_list, parameters = read_simput_catalog(input_events)

    try:
        instrument_spec = instrument_registry[instrument]
    except KeyError:
        raise KeyError("Instrument %s is not in the instrument registry!" % instrument)
    arf_file = check_file_location(instrument_spec["arf"], "files")
    rmf_file = check_file_location(instrument_spec["rmf"], "files")
    arf = AuxiliaryResponseFile(arf_file)
    rmf = RedistributionMatrixFile(rmf_file)

    nx = instrument_spec["num_pixels"]
    plate_scale = instrument_spec["fov"]/nx/60. # arcmin to deg
    plate_scale_arcsec = plate_scale * 3600.0
    dsize = dither_size/plate_scale_arcsec

    event_params = {}
    event_params["exposure_time"] = exp_time
    event_params["arf"] = arf.filename
    event_params["sky_center"] = sky_center
    event_params["pix_center"] = np.array([0.5*(nx+1)]*2)
    event_params["num_pixels"] = nx
    event_params["plate_scale"] = plate_scale
    event_params["rmf"] = rmf.filename
    event_params["channel_type"] = rmf.header["CHANTYPE"]
    event_params["telescope"] = rmf.header["TELESCOP"]
    event_params["instrument"] = rmf.header["INSTRUME"]
    event_params["mission"] = rmf.header.get("MISSION", "")
    event_params["nchan"] = rmf.ebounds_header["DETCHANS"]
    event_params["roll_angle"] = roll_angle
    event_params["fov"] = instrument_spec["fov"]
    event_params["chan_lim"] = [rmf.cmin, rmf.cmax]

    w = pywcs.WCS(naxis=2)
    w.wcs.crval = event_params["sky_center"]
    w.wcs.crpix = event_params["pix_center"]
    w.wcs.cdelt = [-plate_scale, plate_scale]
    w.wcs.ctype = ["RA---TAN","DEC--TAN"]
    w.wcs.cunit = ["deg"]*2

    roll_angle = np.deg2rad(roll_angle)
    rot_mat = np.array([[np.sin(roll_angle), -np.cos(roll_angle)],
                        [-np.cos(roll_angle), -np.sin(roll_angle)]])

    all_events = defaultdict(list)

    for i, evts in enumerate(event_list):

        mylog.info("Detecting events from source %s." % parameters["sources"][i])

        # Step 1: Use ARF to determine which photons are observed

        mylog.info("Applying energy-dependent effective area from %s." % os.path.split(arf.filename)[-1])
        refband = [parameters["emin"][i], parameters["emax"][i]]
        events = arf.detect_events(evts, exp_time, parameters["flux"][i], refband, prng=prng)

        n_evt = events["energy"].size

        if n_evt == 0:
            mylog.warning("No events were observed for this source!!!")
        else:

            # Step 2: Assign pixel coordinates to events. Apply dithering and
            # PSF. Clip events that don't fall within the detection region.

            mylog.info("Pixeling events.")

            # Convert RA, Dec to pixel coordinates
            xpix, ypix = w.wcs_world2pix(events["ra"], events["dec"], 1)

            xpix -= event_params["pix_center"][0]
            ypix -= event_params["pix_center"][1]

            events.pop("ra")
            events.pop("dec")

            n_evt = xpix.size

            # Dither pixel coordinates

            x_offset = np.zeros(n_evt)
            y_offset = np.zeros(n_evt)

            if instrument_spec["dither"]:
                if dither_shape == "circle":
                    r = dsize*np.sqrt(prng.uniform(size=n_evt))
                    theta = 2.*np.pi*prng.uniform(size=n_evt)
                    x_offset = r*np.cos(theta)
                    y_offset = r*np.sin(theta)
                elif dither_shape == "square":
                    x_offset = dsize*prng.uniform(low=-0.5, high=0.5, size=n_evt)
                    y_offset = dsize*prng.uniform(low=-0.5, high=0.5, size=n_evt)

            xpix -= x_offset
            ypix -= y_offset

            # Rotate physical coordinates to detector coordinates

            det = np.dot(rot_mat, np.array([xpix, ypix]))
            detx = det[0,:]
            dety = det[1,:]

            # PSF scattering of detector coordinates

            if instrument_spec["psf"] is not None:
                psf_type, psf_spec = instrument_spec["psf"]
                if psf_type == "gaussian":
                    sigma = psf_spec/sigma_to_fwhm/plate_scale_arcsec
                    detx += prng.normal(loc=0.0, scale=sigma, size=n_evt)
                    dety += prng.normal(loc=0.0, scale=sigma, size=n_evt)
                else:
                    raise NotImplementedError("PSF type %s not implemented!" % psf_type)

            # Convert detector coordinates to chip coordinates

            events["chipx"] = np.round(detx + event_params['pix_center'][0])
            events["chipy"] = np.round(dety + event_params['pix_center'][1])

            # Throw out events that don't fall on the chip

            keepx = np.logical_and(events["chipx"] >= 1.0, events["chipx"] <= nx)
            keepy = np.logical_and(events["chipy"] >= 1.0, events["chipy"] <= nx)
            keep = np.logical_and(keepx, keepy)
            mylog.info("%d events were rejected because " % (n_evt-keep.sum()) +
                       "they fall outside the field of view.")
            n_evt = keep.sum()

            if n_evt == 0:
                mylog.warning("No events are within the field of view for this source!!!")
            else:

                for key in events:
                    events[key] = events[key][keep]

                # Convert chip coordinates back to detector coordinates

                events["detx"] = events["chipx"] - event_params['pix_center'][0] + \
                    prng.uniform(low=-0.5, high=0.5, size=n_evt)
                events["dety"] = events["chipy"] - event_params['pix_center'][1] + \
                    prng.uniform(low=-0.5, high=0.5, size=n_evt)

                # Convert detector coordinates back to pixel coordinates

                pix = np.dot(rot_mat, np.array([events["detx"], events["dety"]]))

                events["xpix"] = pix[0,:] + event_params['pix_center'][0] + x_offset[keep]
                events["ypix"] = pix[1,:] + event_params['pix_center'][1] + y_offset[keep]

        if n_evt > 0:
            for key in events:
                all_events[key] = np.concatenate([all_events[key], events[key]])

    if len(all_events["energy"]) == 0:
        mylog.warning("No events from any of the sources in the catalog were detected!")
        for key in ["xpix", "ypix", "chipx", "chipy", "detx", "dety", "time", 
                    event_params["channel_type"]]:
            all_events[key] = np.array([])
    else:
        # Step 4: Scatter energies with RMF
        mylog.info("Scattering energies with RMF %s." % os.path.split(rmf.filename)[-1])
        all_events = rmf.scatter_energies(all_events, prng=prng)

        # Step 5: Add times to events
        all_events['time'] = prng.uniform(size=all_events["energy"].size, low=0.0,
                                          high=event_params["exposure_time"])

    return all_events, event_params


def make_background(exp_time, instrument, sky_center, foreground=True, 
                    ptsrc_bkgnd=True, instr_bkgnd=True, dither_shape="square", 
                    dither_size=16.0, roll_angle=0.0, prng=None):
    """
    Make background events. 

    Parameters
    ----------
    exp_time : float
        The exposure time to use, in seconds. 
    instrument : string
        The name of the instrument to use, which picks an instrument
        specification from the instrument registry. 
    sky_center : array, tuple, or list
        The center RA, Dec coordinates of the observation, in degrees.
    foreground : boolean, optional
        Whether or not to include the Galactic foreground. Default: True
    instr_bkgnd : boolean, optional
        Whether or not to include the instrumental background. Default: True
    ptsrc_bkgnd : boolean, optional
        Whether or not to include the point-source background. Default: True
        Default: 0.05
    dither_shape : string
        The shape of the dither. Currently "circle" or "square" 
        Default: "square"
    dither_size : float
        The size of the dither in arcseconds. Width of square or radius
        of circle. Default: 16.0
    roll_angle : float
        The roll angle of the observation in degrees. Default: 0.0
    prng : :class:`~numpy.random.RandomState` object, integer, or None
        A pseudo-random number generator. Typically will only 
        be specified if you have a reason to generate the same 
        set of random numbers, such as for a test. Default is None, 
        which sets the seed based on the system time. 
    """
    prng = parse_prng(prng)
    try:
        instrument_spec = instrument_registry[instrument]
    except KeyError:
        raise KeyError("Instrument %s is not in the instrument registry!" % instrument)
    fov = instrument_spec["fov"]

    input_events = defaultdict(list)

    arf_file = check_file_location(instrument_spec["arf"], "files")
    arf = AuxiliaryResponseFile(arf_file)
    rmf_file = check_file_location(instrument_spec["rmf"], "files")
    rmf = RedistributionMatrixFile(rmf_file)

    if ptsrc_bkgnd:
        mylog.info("Adding in point-source background.")
        ptsrc_events = make_ptsrc_background(exp_time, fov, sky_center, 
                                             prng=prng)
        for key in ["ra", "dec", "energy"]:
            input_events[key].append(ptsrc_events[key])
        input_events["flux"].append(ptsrc_events["flux"])
        input_events["emin"].append(ptsrc_events["energy"].min())
        input_events["emax"].append(ptsrc_events["energy"].max())
        input_events["sources"].append("ptsrc_bkgnd")
        events, event_params = generate_events(input_events, exp_time,
                                               instrument, sky_center,
                                               dither_shape=dither_shape, 
                                               dither_size=dither_size, 
                                               roll_angle=roll_angle, prng=prng)
        mylog.info("Generated %d photons from the point-source background." % len(events["energy"]))
    else:
        nx = instrument_spec["num_pixels"]
        events = defaultdict(list)
        event_params = {"exposure_time": exp_time, 
                        "fov": instrument_spec["fov"],
                        "num_pixels": nx,
                        "pix_center": np.array([0.5*(nx+1)]*2),
                        "channel_type": rmf.header["CHANTYPE"],
                        "sky_center": sky_center,
                        "plate_scale": instrument_spec["fov"]/nx/60.,
                        "chan_lim": [rmf.cmin, rmf.cmax],
                        "rmf": rmf_file, "arf": arf_file,
                        "telescope": rmf.header["TELESCOP"],
                        "instrument": rmf.header["INSTRUME"],
                        "mission": rmf.header.get("MISSION", ""),
                        "nchan": rmf.ebounds_header["DETCHANS"],
                        "roll_angle": 0.0}

    if foreground:
        mylog.info("Adding in astrophysical foreground.")
        bkg_events = make_foreground(event_params, arf, rmf, prng=prng)
        for key in bkg_events:
            events[key] = np.concatenate([events[key], bkg_events[key]])
    if instr_bkgnd:
        mylog.info("Adding in instrumental background.")
        bkg_events = make_instrument_background(instrument_spec["bkgnd"], 
                                                event_params,
                                                instrument_spec["focal_length"], 
                                                rmf, prng=prng)
        for key in bkg_events:
            events[key] = np.concatenate([events[key], bkg_events[key]])

    return events, event_params

def make_background_file(out_file, exp_time, instrument, sky_center, 
                         overwrite=False, foreground=True, instr_bkgnd=True,
                         ptsrc_bkgnd=True, dither_shape="square", 
                         dither_size=16.0, prng=None):
    """
    Make an event file consisting entirely of background events. This will be 
    useful for creating backgrounds that can be added to simulations of sources.

    Parameters
    ----------
    exp_time : float
        The exposure time to use, in seconds. 
    instrument : string
        The name of the instrument to use, which picks an instrument
        specification from the instrument registry. 
    sky_center : array, tuple, or list
        The center RA, Dec coordinates of the observation, in degrees.
    overwrite : boolean, optional
        Whether or not to overwrite an existing file with the same name.
        Default: False
    foreground : boolean, optional
        Whether or not to include the Galactic foreground. Default: True
    cosmo_bkgnd : boolean, optional
        Whether or not to include the cosmological halo background. Default:
        True
    instr_bkgnd : boolean, optional
        Whether or not to include the instrumental background. Default: True
    ptsrc_bkgnd : boolean, optional
        Whether or not to include the point-source background. Default: True
    dither_shape : string
        The shape of the dither. Currently "circle" or "square" 
        Default: "square"
    dither_size : float
        The size of the dither in arcseconds. Width of square or radius
        of circle. Default: 16.0
    prng : :class:`~numpy.random.RandomState` object, integer, or None
        A pseudo-random number generator. Typically will only 
        be specified if you have a reason to generate the same 
        set of random numbers, such as for a test. Default is None, 
        which sets the seed based on the system time. 
    """
    prng = parse_prng(prng)
    events, event_params = make_background(exp_time, instrument, sky_center, 
                                           ptsrc_bkgnd=ptsrc_bkgnd, 
                                           foreground=foreground, 
                                           instr_bkgnd=instr_bkgnd,
                                           dither_shape=dither_shape, 
                                           dither_size=dither_size, prng=prng)
    write_event_file(events, event_params, out_file, overwrite=overwrite)

def instrument_simulator(input_events, out_file, exp_time, instrument,
                         sky_center, overwrite=False, instr_bkgnd=True, 
                         foreground=True, ptsrc_bkgnd=True, 
                         bkgnd_file=None, dither_shape="square", 
                         dither_size=16.0, roll_angle=0.0, prng=None):
    """
    Take unconvolved events and create an event file from them. This
    function calls generate_events to do the following:

    1. Determines which events are observed using the ARF
    2. Pixelizes the events, applying PSF effects and dithering
    3. Determines energy channels using the RMF

    and then calls make_background to add instrumental and astrophysical
    backgrounds, unless a background file is provided, in which case
    the background events are read from this file. The events are
    then written out to a file.

    Parameters
    ----------
    input_events : string, dict, or None
        The unconvolved events to be used as input. Can be one of the
        following:
        1. The name of a SIMPUT catalog file.
        2. A Python dictionary containing the following items:
        "ra": A NumPy array of right ascension values in degrees.
        "dec": A NumPy array of declination values in degrees.
        "energy": A NumPy array of energy values in keV.
        "flux": The flux of the entire source, in units of erg/cm**2/s.
    out_file : string
        The name of the event file to be written.
    exp_time : float
        The exposure time to use, in seconds. 
    instrument : string
        The name of the instrument to use, which picks an instrument
        specification from the instrument registry. 
    sky_center : array, tuple, or list
        The center RA, Dec coordinates of the observation, in degrees.
    overwrite : boolean, optional
        Whether or not to overwrite an existing file with the same name.
        Default: False
    instr_bkgnd : boolean, optional
        Whether or not to include the instrumental/particle background. 
        Default: True
    foreground : boolean, optional
        Whether or not to include the local foreground. 
        Default: True
    ptsrc_bkgnd : boolean, optional
        Whether or not to include the point-source background. 
        Default: True
    bkgnd_file : string, optional
        If set, backgrounds will be loaded from this file and not generated
        on the fly. Default: None
    dither_shape : string
        The shape of the dither. Currently "circle" or "square" 
        Default: "square"
    dither_size : float
        The size of the dither in arcseconds. Width of square or radius
        of circle. Default: 16.0
    roll_angle : float
        The roll angle of the observation in degrees. Default: 0.0
    prng : :class:`~numpy.random.RandomState` object, integer, or None
        A pseudo-random number generator. Typically will only 
        be specified if you have a reason to generate the same 
        set of random numbers, such as for a test. Default is None, 
        which sets the seed based on the system time. 

    Examples
    --------
    >>> instrument_simulator("sloshing_simput.fits", "sloshing_evt.fits", 
    ...                      300000.0, "hdxi_3x10", [30., 45.], overwrite=True)
    """
    if not out_file.endswith(".fits"):
        out_file += ".fits"
    mylog.info("Making observation of source in %s." % out_file)
    # Make the source first
    events, event_params = generate_events(input_events, exp_time, instrument, sky_center,
                                           dither_shape=dither_shape, dither_size=dither_size, 
                                           roll_angle=roll_angle, prng=prng)
    # If the user wants backgrounds, either make the background or add an already existing
    # background event file. It may be necessary to reproject events to a new coordinate system.
    if bkgnd_file is None:
        if not instr_bkgnd and not ptsrc_bkgnd and not foreground:
            mylog.info("No backgrounds will be added to this observation.")
        else:
            mylog.info("Adding background events.")
            bkg_events, _ = make_background(exp_time, instrument, sky_center,
                                            foreground=foreground, instr_bkgnd=instr_bkgnd, 
                                            dither_shape=dither_shape, dither_size=dither_size, 
                                            ptsrc_bkgnd=ptsrc_bkgnd, prng=prng, 
                                            roll_angle=roll_angle)
            for key in events:
                events[key] = np.concatenate([events[key], bkg_events[key]])
    else:
        mylog.info("Adding background events from the file %s." % bkgnd_file)
        if not os.path.exists(bkgnd_file):
            raise IOError("Cannot find the background event file %s!" % bkgnd_file)
        events = add_background_from_file(events, event_params, bkgnd_file)
    write_event_file(events, event_params, out_file, overwrite=overwrite)
    mylog.info("Observation complete.")

def simulate_spectrum(spec, instrument, exp_time, out_file, overwrite=False,
                      prng=None):
    """
    Generate a PI or PHA spectrum from a :class:`~soxs.spectra.Spectrum`
    by convolving it with responses. To be used if one wants to 
    create a spectrum without worrying about spatial response. Similar
    to XSPEC's "fakeit". 

    Parameters
    ----------
    spec : :class:`~soxs.spectra.Spectrum`
        The spectrum to be convolved.
    instrument : string
        The name of the instrument to use, which picks an instrument
        specification from the instrument registry. 
    exp_time : float
        The exposure time in seconds.
    out_file : string
        The file to write the spectrum to.
    overwrite : boolean, optional
        Whether or not to overwrite an existing file. Default: False
    prng : :class:`~numpy.random.RandomState` object, integer, or None
        A pseudo-random number generator. Typically will only 
        be specified if you have a reason to generate the same 
        set of random numbers, such as for a test. Default is None, 
        which sets the seed based on the system time. 

    Examples
    --------
    >>> spec = soxs.Spectrum.from_file("my_spectrum.txt")
    >>> soxs.simulate_spectrum(spec, "mucal", 100000.0, 
    ...                        "my_spec.pi", overwrite=True)
    """
    from soxs.events import write_spectrum
    from soxs.instrument import RedistributionMatrixFile, \
        AuxiliaryResponseFile
    from soxs.spectra import ConvolvedSpectrum
    prng = parse_prng(prng)
    try:
        instrument_spec = instrument_registry[instrument]
    except KeyError:
        raise KeyError("Instrument %s is not in the instrument registry!" % instrument)
    arf_file = check_file_location(instrument_spec["arf"], "files")
    rmf_file = check_file_location(instrument_spec["rmf"], "files")
    arf = AuxiliaryResponseFile(arf_file)
    rmf = RedistributionMatrixFile(rmf_file)
    cspec = ConvolvedSpectrum(spec, arf)
    events = {}
    events["energy"] = cspec.generate_energies(exp_time, prng=prng).value
    events = rmf.scatter_energies(events, prng=prng)
    events["arf"] = arf.filename
    events["rmf"] = rmf.filename
    events["exposure_time"] = exp_time
    events["channel_type"] = rmf.header["CHANTYPE"]
    events["telescope"] = rmf.header["TELESCOP"]
    events["instrument"] = rmf.header["INSTRUME"]
    events["mission"] = rmf.header.get("MISSION", "")
    write_spectrum(events, out_file, overwrite=overwrite)

