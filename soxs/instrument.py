from __future__ import print_function

import json
import numpy as np
import astropy.io.fits as pyfits
import astropy.wcs as pywcs
import os

from astropy.utils.console import ProgressBar
from soxs.constants import erg_per_keV
from soxs.simput import read_simput_catalog
from soxs.utils import mylog, check_file_location, \
    ensure_numpy_array, write_event_file
from soxs.background import background_registry

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
        Interpolate the effective area to the energies provided by the supplied *energy* array.
        """
        earea = np.interp(energy, self.emid, self.eff_area, left=0.0, right=0.0)
        return earea

    def detect_events(self, events, exp_time, flux, refband, prng=None):
        """
        Use the ARF to determine a subset of photons which will be
        detected. Returns a boolean NumPy array which is the same
        is the same size as the number of photons, wherever it is
        "true" means those photons have been detected.

        Parameters
        ----------
        events : dict of np.ndarrays
            The energies and positions of the photons. 
        exp_time : float
            The exposure time in seconds.
        flux : float
            The total flux of the photons in erg/s/cm^2. 
        refband : array_like
            A two-element array or list containing the limits of the energy band
            which the flux was computed in. 
        prng : :class:`~numpy.random.RandomState` object or :mod:`~numpy.random`, optional
            A pseudo-random number generator. Typically will only be specified
            if you have a reason to generate the same set of random numbers, such as for a
            test. Default is the :mod:`~numpy.random` module.
        """
        if prng is None:
            prng = np.random
        energy = events["energy"]
        earea = self.interpolate_area(energy)
        idxs = np.logical_and(energy >= refband[0], energy <= refband[1])
        rate = flux/(energy[idxs].sum()*erg_per_keV)*earea[idxs].sum()
        n_ph = np.uint64(rate*exp_time)
        fak = float(n_ph)/energy.size
        if fak > 1.0:
            mylog.error("Number of events in sample: %d, Number of events wanted: %d" % (energy.size, n_ph))
            raise ValueError("This combination of exposure time and effective area "
                             "will result in more photons being drawn than are available "
                             "in the sample!!!")
        w = earea / self.max_area
        randvec = prng.uniform(size=energy.size)
        eidxs = prng.permutation(np.where(randvec < w)[0])[:n_ph].astype("uint64")
        mylog.info("%s events detected." % n_ph)
        for key in events:
            events[key] = events[key][eidxs]
        return events

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

    def scatter_energies(self, events, prng=np.random):
        """
        Scatter photon energies with the RMF and produce the 
        corresponding channel values.

        Parameters
        ----------
        events : dict of np.ndarrays
            The energies and positions of the photons. 
        prng : :class:`~numpy.random.RandomState` object or :mod:`~numpy.random`, optional
            A pseudo-random number generator. Typically will only be specified
            if you have a reason to generate the same set of random numbers, such as for a
            test. Default is the :mod:`~numpy.random` module.
        """
        elo = self.data["ENERG_LO"]
        ehi = self.data["ENERG_HI"]
        n_de = elo.shape[0]
        mylog.info("Number of energy bins in RMF: %d" % n_de)
        mylog.info("Energy limits: %g %g" % (min(elo), max(ehi)))

        n_ch = len(self.ebounds["CHANNEL"])
        mylog.info("Number of channels in RMF: %d" % n_ch)

        eidxs = np.argsort(events["energy"])
        sorted_e = events["energy"][eidxs]

        detectedChannels = []

        # run through all photon energies and find which bin they go in
        fcurr = 0
        last = sorted_e.shape[0]

        with ProgressBar(last) as pbar:
            for (k, low), high in zip(enumerate(elo), ehi):
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
                    pbar.update(fcurr)

        for key in events:
            events[key] = events[key][eidxs]
        events[self.header["CHANTYPE"]] = np.concatenate(detectedChannels)

        return events

instrument_registry = {}
instrument_registry["xcal"] = {"name": "xcal",
                               "arf": "xrs_calorimeter.arf",
                               "rmf": "xrs_calorimeter.rmf",
                               "bkgnd": "acisi",
                               "num_pixels": 300,
                               "plate_scale": 1.0,
                               "psf": ["gaussian", 0.5]}
instrument_registry["hdxi"] = {"name": "hdxi",
                               "arf": "xrs_hdxi.arf",
                               "rmf": "xrs_hdxi.rmf",
                               "bkgnd": "acisi",
                               "num_pixels": 4096,
                               "plate_scale": 1./3.,
                               "psf": ["gaussian", 0.5]}

def add_instrument_to_registry(inst_spec):
    """
    Add an instrument specification to the registry, contained
    in either a dictionary or a JSON file.

    The *inst_spec* must have the structure as shown below. 
    The order is not important. If you use a JSON file, the
    structure is the same, but the file cannot include comments.

    >>> {
    ...     "name": "hdxi", # The short name of the instrument
    ...     "arf": "xrs_hdxi.arf", # The file containing the ARF
    ...     "rmf": "xrs_hdxi.rmf" # The file containing the RMF
    ...     "bkgnd": "acisi_particle_bkgnd.dat" # The file containing the particle background
    ...     "plate_scale": 0.33333333333, # The plate scale in arcsec
    ...     "num_pixels": 4096, # The number of pixels on a side in the FOV
    ...     "psf": ["gaussian", 0.5] # The type of PSF and its FWHM
    ... }
    """
    if isinstance(inst_spec, dict):
        inst = inst_spec
    elif os.path.exists(inst_spec):
        f = open(inst_spec, "r")
        inst = json.load(f)
        f.close()
    name = inst["name"]
    if name in instrument_registry:
        raise KeyError("The instrument with name %s is already in the registry! Assign a different name!" % name)
    instrument_registry[name] = inst
    mylog.info("The %s instrument specification has been added to the instrument registry." % name)
    return name

def get_instrument_from_registry(name):
    """
    Returns a copy of the instrument specification
    corresponding to *name*.
    """
    return instrument_registry[name].copy()

def show_instrument_registry():
    """
    Print the contents of the instrument registry.
    """
    for name, spec in instrument_registry.items():
        print("Instrument: %s" % name)
        for k, v in spec.items():
            print("    %s: %s" % (k, v))

def write_instrument_json(inst_name, filename):
    """
    Write an instrument specification to a JSON file.
    Useful if one would like to create a new specification
    by editing an existing one. 

    Parameters
    ----------
    inst_name : string
        The instrument specification to write.
    filename : string
        The filename to write to. 
    """
    inst_dict = instrument_registry[inst_name]
    fp = open(filename, 'w')
    json.dump(inst_dict, fp, indent=4)
    fp.close()

def add_background(bkgnd_name, event_params, rot_mat, bkgnd_scale, prng=np.random):

    fov = event_params["num_pixels"]*event_params["plate_scale"]*60.0

    bkgnd_spec = background_registry[bkgnd_name]

    # Generate background events

    bkg_events = {}

    # If this is an astrophysical background, detect the events with the ARF used
    # for the observation, otherwise it's a particle background and assume area = 1.
    if bkgnd_spec.bkgnd_type == "astrophysical":
        arf = AuxiliaryResponseFile(check_file_location(event_params["arf"], "files"))
        area = arf.interpolate_area(bkgnd_spec.emid)
    else:
        area = 1.0

    bkg_events["energy"] = bkgnd_spec.generate_energies(event_params["exposure_time"],
                                                        area, fov, bkgnd_scale, prng=prng)

    n_events = bkg_events["energy"].size

    bkg_events['chipx'] = np.round(prng.uniform(low=1.0, high=event_params['num_pixels'],
                                                size=n_events))
    bkg_events['chipy'] = np.round(prng.uniform(low=1.0, high=event_params['num_pixels'],
                                                size=n_events))
    bkg_events["detx"] = np.round(bkg_events["chipx"] - event_params['pix_center'][0] +
                                  prng.uniform(low=-0.5, high=0.5, size=n_events))
    bkg_events["dety"] = np.round(bkg_events["chipy"] - event_params['pix_center'][1] +
                                  prng.uniform(low=-0.5, high=0.5, size=n_events))
    pix = np.dot(rot_mat, np.array([bkg_events["detx"], bkg_events["dety"]]))
    bkg_events["xpix"] = pix[0,:] + event_params['pix_center'][0]
    bkg_events["ypix"] = pix[1,:] + event_params['pix_center'][1]

    return bkg_events

def instrument_simulator(simput_file, out_file, exp_time, instrument,
                         sky_center, clobber=False, dither_shape="square",
                         dither_size=16.0, roll_angle=0.0, astro_bkgnd="hm_cxb",
                         instr_bkgnd_scale=1.0, prng=np.random):
    """
    Take unconvolved events in a SIMPUT catalog and create an event
    file from them. This function does the following:

    1. Determines which events are observed using the ARF
    2. Pixelizes the events, applying PSF effects and dithering
    3. Adds instrumental and astrophysical background events
    4. Determines energy channels using the RMF
    5. Writes the events to a file

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
    dither_shape : string
        The shape of the dither. Currently "circle" or "square" 
        Default: "square"
    dither_size : float
        The size of the dither in arcseconds. Width of square or radius
        of circle. Default: 16.0
    roll_angle : float
        The roll angle of the observation in degrees. Default: 0.0
    astro_bkgnd : string
        The astrophysical background spectrum to assume. Must be the name of a
        background registered in the background registry. If None, no astrophysical
        background is applied. Default: "hm_cxb"
    instr_bkgnd_scale : float, optional
        The scale of the instrumental background. Default: 1.0
    prng : :class:`~numpy.random.RandomState` object or :mod:`~numpy.random`, optional
        A pseudo-random number generator. Typically will only be specified
        if you have a reason to generate the same set of random numbers, such as for a
        test. Default is the :mod:`numpy.random` module.

    Examples
    --------
    >>> instrument_simulator("sloshing_simput.fits", "sloshing_evt.fits", "hdxi",
    ...                      [30., 45.], clobber=True)
    """
    event_list, parameters = read_simput_catalog(simput_file)

    try:
        instrument_spec = instrument_registry[instrument]
    except KeyError:
        raise KeyError("Instrument %s is not in the instrument registry!" % instrument)
    arf_file = check_file_location(instrument_spec["arf"], "files")
    rmf_file = check_file_location(instrument_spec["rmf"], "files")
    arf = AuxiliaryResponseFile(arf_file)
    rmf = RedistributionMatrixFile(rmf_file)

    nx = instrument_spec["num_pixels"]
    plate_scale = instrument_spec["plate_scale"]/3600. # arcsec to deg

    event_params = {}
    event_params["exposure_time"] = exp_time
    event_params["arf"] = os.path.split(arf.filename)[-1]
    event_params["sky_center"] = sky_center
    event_params["pix_center"] = np.array([0.5*(nx+1)]*2)
    event_params["num_pixels"] = nx
    event_params["plate_scale"] = plate_scale
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
    w.wcs.cdelt = [-plate_scale, plate_scale]
    w.wcs.ctype = ["RA---TAN","DEC--TAN"]
    w.wcs.cunit = ["deg"]*2

    roll_angle = np.deg2rad(roll_angle)
    rot_mat = np.array([[np.sin(roll_angle), -np.cos(roll_angle)],
                        [-np.cos(roll_angle), -np.sin(roll_angle)]])

    all_events = {}

    first = True

    for i, evts in enumerate(event_list):

        mylog.info("Detecting events from source %d" % (i+1))

        # Step 1: Use ARF to determine which photons are observed

        mylog.info("Applying energy-dependent effective area from %s." % event_params["arf"])
        refband = [parameters["emin"][i], parameters["emax"][i]]
        events = arf.detect_events(evts, exp_time, parameters["flux"][i],
                                   refband, prng=prng)

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

            x_offset = 0.0
            y_offset = 0.0
            if dither_shape == "circle":
                r = dither_size*np.sqrt(prng.uniform(size=n_evt))
                theta = 2.*np.pi*prng.uniform(size=n_evt)
                x_offset = r*np.cos(theta)
                y_offset = r*np.sin(theta)
            elif dither_shape == "square":
                x_offset = dither_size*prng.uniform(low=-0.5, high=0.5, size=n_evt)
                y_offset = dither_size*prng.uniform(low=-0.5, high=0.5, size=n_evt)

            xpix -= x_offset
            ypix -= y_offset

            # Rotate physical coordinates to detector coordinates

            det = np.dot(rot_mat, np.array([xpix, ypix]))
            detx = det[0,:]
            dety = det[1,:]

            # PSF scattering of detector coordinates

            psf_type, psf_spec = instrument_spec["psf"]
            if psf_type == "gaussian":
                sigma = psf_spec/sigma_to_fwhm/instrument_spec["plate_scale"]
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

                events["detx"] = np.round(events["chipx"] - event_params['pix_center'][0] +
                                          prng.uniform(low=-0.5, high=0.5, size=n_evt))
                events["dety"] = np.round(events["chipy"] - event_params['pix_center'][1] +
                                          prng.uniform(low=-0.5, high=0.5, size=n_evt))

                # Convert detector coordinates back to pixel coordinates

                pix = np.dot(rot_mat, np.array([events["detx"], events["dety"]]))

                events["xpix"] = pix[0,:] + event_params['pix_center'][0] + x_offset[keep]
                events["ypix"] = pix[1,:] + event_params['pix_center'][1] + y_offset[keep]

        if n_evt > 0:
            for key in events:
                if first:
                    all_events[key] = events[key]
                else:
                    all_events[key] = np.concatenate([all_events[key], events[key]])
            first = False

    if all_events["energy"].size == 0:
        mylog.warning("No events from any of the sources in the catalog were detected!")

    # Step 3: Add astrophysical background

    if astro_bkgnd is not None:

        mylog.info("Adding in astrophysical background.")

        bkg_events = add_background(astro_bkgnd, event_params, rot_mat, 1.0, prng=prng)

        for key in all_events:
            all_events[key] = np.concatenate([all_events[key], bkg_events[key]])

    # Step 4: Add particle background

    if instr_bkgnd_scale > 0.0:

        mylog.info("Adding in instrumental background.")

        bkg_events = add_background(instrument_spec["bkgnd"], event_params, rot_mat,
                                    instr_bkgnd_scale, prng=prng)

        for key in all_events:
            all_events[key] = np.concatenate([all_events[key], bkg_events[key]])

    if all_events["energy"].size == 0:
        raise RuntimeError("No events were detected!!!")

    # Step 5: Scatter energies with RMF

    if all_events["energy"].size > 0:
        mylog.info("Scattering energies with RMF %s." % event_params['rmf'])
        all_events = rmf.scatter_energies(all_events, prng=prng)

    # Step 6: Add times to events

    all_events['time'] = np.random.uniform(size=all_events["energy"].size, low=0.0,
                                           high=event_params["exposure_time"])

    write_event_file(all_events, event_params, out_file, clobber=clobber)