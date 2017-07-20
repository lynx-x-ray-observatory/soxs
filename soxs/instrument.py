import numpy as np
import astropy.io.fits as pyfits
import astropy.wcs as pywcs
import astropy.units as u
import os
from collections import defaultdict

from soxs.constants import erg_per_keV
from soxs.simput import read_simput_catalog
from soxs.utils import mylog, check_file_location, \
    ensure_numpy_array, parse_prng, parse_value, \
    get_rot_mat
from soxs.events import write_event_file
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
        n_ph = prng.poisson(lam=rate*exp_time)
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

    def plot(self, fig=None, ax=None):
        """
        Make a quick plot of the effective area curve.

        Parameters
        ----------
        fig : :class:`~matplotlib.figure.Figure`, optional
            The figure to place the plot in. If not supplied, one will be created.
        ax : :class:`~matplotlib.axes.Axes`, optional
            The axes to place the plot in. If not supplied, one will be created.
        """
        import matplotlib.pyplot as plt
        if fig is None:
            fig = plt.figure(figsize=(10, 10))
        if ax is None:
            ax = fig.add_subplot(111)
        ax.loglog(self.emid, self.eff_area)
        ax.set_xlabel("E (keV)")
        ax.set_ylabel("$\mathrm{A\ (cm^2)}$")
        return fig

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

    def _make_channels(self, k):
        # build channel number list associated to array value,
        # there are groups of channels in rmfs with nonzero probabilities
        trueChannel = []
        f_chan = ensure_numpy_array(np.nan_to_num(self.data["F_CHAN"][k]))
        n_chan = ensure_numpy_array(np.nan_to_num(self.data["N_CHAN"][k]))
        for start, nchan in zip(f_chan, n_chan):
            if nchan == 0:
                trueChannel.append(start)
            else:
                trueChannel += list(range(start, start + nchan))
        return np.array(trueChannel)

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
            trueChannel = self._make_channels(k)
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

def perform_dither(t, dither_dict):
    if dither_dict["dither_on"]:
        a = 2.0*np.pi/dither_dict["x_period"]
        b = 2.0*np.pi/dither_dict["y_period"]
        A = dither_dict["x_amp"]/dither_dict["plate_scale"]
        B = dither_dict["y_amp"]/dither_dict["plate_scale"]
        x_offset = A*np.sin(a*t)
        y_offset = B*np.sin(b*t)
    else:
        x_offset = np.zeros(t.size)
        y_offset = np.zeros(t.size)
    return x_offset, y_offset

def generate_events(input_events, exp_time, instrument, sky_center, 
                    no_dither=False, dither_params=None, 
                    roll_angle=0.0, subpixel_res=False, prng=None):
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
    exp_time : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The exposure time to use, in seconds. 
    instrument : string
        The name of the instrument to use, which picks an instrument
        specification from the instrument registry. 
    sky_center : array, tuple, or list
        The center RA, Dec coordinates of the observation, in degrees.
    no_dither : boolean, optional
        If True, turn off dithering entirely. Default: False
    dither_params : array-like of floats, optional
        The parameters to use to control the size and period of the dither
        pattern. The first two numbers are the dither amplitude in x and y
        detector coordinates in arcseconds, and the second two numbers are
        the dither period in x and y detector coordinates in seconds. 
        Default: [8.0, 8.0, 1000.0, 707.0].
    roll_angle : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
        The roll angle of the observation in degrees. Default: 0.0
    subpixel_res: boolean, optional
        If True, event positions are not randomized within the pixels 
        within which they are detected. Default: False
    prng : :class:`~numpy.random.RandomState` object, integer, or None
        A pseudo-random number generator. Typically will only 
        be specified if you have a reason to generate the same 
        set of random numbers, such as for a test. Default is None, 
        which sets the seed based on the system time. 
    """
    import pyregion._region_filter as rfilter
    exp_time = parse_value(exp_time, "s")
    roll_angle = parse_value(roll_angle, "deg")
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

    if not instrument_spec["dither"]:
        dither_on = False
    else:
        dither_on = not no_dither
    if dither_params is None:
        dither_params = [8.0, 8.0, 1000.0, 707.0]
    dither_dict = {"x_amp": dither_params[0],
                   "y_amp": dither_params[1],
                   "x_period": dither_params[2],
                   "y_period": dither_params[3],
                   "dither_on": dither_on,
                   "plate_scale": plate_scale_arcsec}

    event_params = {}
    event_params["exposure_time"] = exp_time
    event_params["arf"] = arf.filename
    event_params["sky_center"] = sky_center
    event_params["pix_center"] = np.array([0.5*(2*nx+1)]*2)
    event_params["num_pixels"] = nx
    event_params["plate_scale"] = plate_scale
    event_params["rmf"] = rmf.filename
    event_params["channel_type"] = rmf.header["CHANTYPE"]
    event_params["telescope"] = rmf.header["TELESCOP"]
    event_params["instrument"] = instrument_spec['name']
    event_params["mission"] = rmf.header.get("MISSION", "")
    event_params["nchan"] = rmf.ebounds_header["DETCHANS"]
    event_params["roll_angle"] = roll_angle
    event_params["fov"] = instrument_spec["fov"]
    event_params["chan_lim"] = [rmf.cmin, rmf.cmax]
    event_params["chips"] = instrument_spec["chips"]
    event_params["dither_params"] = dither_dict
    event_params["aimpt_coords"] = instrument_spec["aimpt_coords"]

    w = pywcs.WCS(naxis=2)
    w.wcs.crval = event_params["sky_center"]
    w.wcs.crpix = event_params["pix_center"]
    w.wcs.cdelt = [-plate_scale, plate_scale]
    w.wcs.ctype = ["RA---TAN","DEC--TAN"]
    w.wcs.cunit = ["deg"]*2

    rot_mat = get_rot_mat(roll_angle)

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

            # Rotate physical coordinates to detector coordinates

            det = np.dot(rot_mat, np.array([xpix, ypix]))
            detx = det[0,:] + event_params["aimpt_coords"][0]
            dety = det[1,:] + event_params["aimpt_coords"][1]

            # Add times to events
            events['time'] = prng.uniform(size=n_evt, low=0.0,
                                          high=event_params["exposure_time"])

            # Apply dithering

            x_offset, y_offset = perform_dither(events["time"], dither_dict)

            detx -= x_offset
            dety -= y_offset

            # PSF scattering of detector coordinates

            if instrument_spec["psf"] is not None:
                psf_type, psf_spec = instrument_spec["psf"]
                if psf_type == "gaussian":
                    sigma = psf_spec/sigma_to_fwhm/plate_scale_arcsec
                    detx += prng.normal(loc=0.0, scale=sigma, size=n_evt)
                    dety += prng.normal(loc=0.0, scale=sigma, size=n_evt)
                else:
                    raise NotImplementedError("PSF type %s not implemented!" % psf_type)

            # Convert detector coordinates to chip coordinates.
            # Throw out events that don't fall on any chip.

            cx = np.trunc(detx)+0.5*np.sign(detx)
            cy = np.trunc(dety)+0.5*np.sign(dety)

            if event_params["chips"] is None:
                events["chip_id"] = np.zeros(n_evt, dtype='int')
                keepx = np.logical_and(cx >= -0.5*nx, cx <= 0.5*nx)
                keepy = np.logical_and(cy >= -0.5*nx, cy <= 0.5*nx)
                keep = np.logical_and(keepx, keepy)
            else:
                events["chip_id"] = -np.ones(n_evt, dtype='int')
                for i, chip in enumerate(event_params["chips"]):
                    thisc = np.ones(n_evt, dtype='bool')
                    rtype = chip[0]
                    args = chip[1:]
                    r = getattr(rfilter, rtype)(*args)
                    inside = r.inside(cx, cy)
                    thisc = np.logical_and(thisc, inside)
                    events["chip_id"][thisc] = i
                keep = events["chip_id"] > -1

            mylog.info("%d events were rejected because " % (n_evt-keep.sum()) +
                       "they do not fall on any CCD.")
            n_evt = keep.sum()

            if n_evt == 0:
                mylog.warning("No events are within the field of view for this source!!!")
            else:

                # Keep only those events which fall on a chip

                for key in events:
                    events[key] = events[key][keep]

                # Convert chip coordinates back to detector coordinates, unless the
                # user has specified that they want subpixel resolution

                if subpixel_res:
                    events["detx"] = detx[keep]
                    events["dety"] = dety[keep]
                else:
                    events["detx"] = cx[keep] + prng.uniform(low=-0.5, high=0.5, size=n_evt)
                    events["dety"] = cy[keep] + prng.uniform(low=-0.5, high=0.5, size=n_evt)

                # Convert detector coordinates back to pixel coordinates by
                # adding the dither offsets back in and applying the rotation
                # matrix again

                det = np.array([events["detx"] + x_offset[keep] - event_params["aimpt_coords"][0],
                                events["dety"] + y_offset[keep] - event_params["aimpt_coords"][1]])
                pix = np.dot(rot_mat.T, det)

                events["xpix"] = pix[0,:] + event_params['pix_center'][0]
                events["ypix"] = pix[1,:] + event_params['pix_center'][1]

        if n_evt > 0:
            for key in events:
                all_events[key] = np.concatenate([all_events[key], events[key]])

    if len(all_events["energy"]) == 0:
        mylog.warning("No events from any of the sources in the catalog were detected!")
        for key in ["xpix", "ypix", "detx", "dety", "time", "chip_id", event_params["channel_type"]]:
            all_events[key] = np.array([])
    else:
        # Step 4: Scatter energies with RMF
        mylog.info("Scattering energies with RMF %s." % os.path.split(rmf.filename)[-1])
        all_events = rmf.scatter_energies(all_events, prng=prng)

    return all_events, event_params


def make_background(exp_time, instrument, sky_center, foreground=True, 
                    ptsrc_bkgnd=True, instr_bkgnd=True, no_dither=False,
                    dither_params=None, roll_angle=0.0, subpixel_res=False, 
                    input_sources=None, absorb_model="wabs", nH=0.05, prng=None):
    """
    Make background events. 

    Parameters
    ----------
    exp_time : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
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
    no_dither : boolean, optional
        If True, turn off dithering entirely. Default: False
    dither_params : array-like of floats, optional
        The parameters to use to control the size and period of the dither
        pattern. The first two numbers are the dither amplitude in x and y
        detector coordinates in arcseconds, and the second two numbers are
        the dither period in x and y detector coordinates in seconds. 
        Default: [8.0, 8.0, 1000.0, 707.0].
    ptsrc_bkgnd : boolean, optional
        Whether or not to include the point-source background. Default: True
        Default: 0.05
    roll_angle : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
        The roll angle of the observation in degrees. Default: 0.0
    subpixel_res: boolean, optional
        If True, event positions are not randomized within the pixels 
        within which they are detected. Default: False
    input_sources : string, optional
        If set to a filename, input the point source positions, fluxes,
        and spectral indices from an ASCII table instead of generating
        them. Default: None
    absorb_model : string, optional
        The absorption model to use, "wabs" or "tbabs". Default: "wabs"
    nH : float, optional
        The hydrogen column in units of 10**22 atoms/cm**2. 
        Default: 0.05
    prng : :class:`~numpy.random.RandomState` object, integer, or None
        A pseudo-random number generator. Typically will only 
        be specified if you have a reason to generate the same 
        set of random numbers, such as for a test. Default is None, 
        which sets the seed based on the system time. 
    """
    from soxs.background import make_instrument_background, \
        make_foreground, make_ptsrc_background
    prng = parse_prng(prng)
    exp_time = parse_value(exp_time, "s")
    roll_angle = parse_value(roll_angle, "deg")
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
                                             input_sources=input_sources, 
                                             absorb_model=absorb_model,
                                             nH=nH, prng=prng)
        for key in ["ra", "dec", "energy"]:
            input_events[key].append(ptsrc_events[key])
        input_events["flux"].append(ptsrc_events["flux"])
        input_events["emin"].append(ptsrc_events["energy"].min())
        input_events["emax"].append(ptsrc_events["energy"].max())
        input_events["sources"].append("ptsrc_bkgnd")
        events, event_params = generate_events(input_events, exp_time,
                                               instrument, sky_center,
                                               no_dither=no_dither,
                                               dither_params=dither_params, 
                                               roll_angle=roll_angle,
                                               subpixel_res=subpixel_res,
                                               prng=prng)
        mylog.info("Generated %d photons from the point-source background." % len(events["energy"]))
    else:
        nx = instrument_spec["num_pixels"]
        events = defaultdict(list)
        if not instrument_spec["dither"]:
            dither_on = False
        else:
            dither_on = not no_dither
        if dither_params is None:
            dither_params = [8.0, 8.0, 1000.0, 707.0]
        dither_dict = {"x_amp": dither_params[0],
                       "y_amp": dither_params[1],
                       "x_period": dither_params[2],
                       "y_period": dither_params[3],
                       "dither_on": dither_on,
                       "plate_scale": instrument_spec["fov"]/nx*60.0}
        event_params = {"exposure_time": exp_time, 
                        "fov": instrument_spec["fov"],
                        "num_pixels": nx,
                        "pix_center": np.array([0.5*(2*nx+1)]*2),
                        "channel_type": rmf.header["CHANTYPE"],
                        "sky_center": sky_center,
                        "dither_params": dither_dict,
                        "plate_scale": instrument_spec["fov"]/nx/60.0,
                        "chan_lim": [rmf.cmin, rmf.cmax],
                        "rmf": rmf_file, "arf": arf_file,
                        "telescope": rmf.header["TELESCOP"],
                        "instrument": instrument_spec['name'],
                        "mission": rmf.header.get("MISSION", ""),
                        "nchan": rmf.ebounds_header["DETCHANS"],
                        "roll_angle": roll_angle,
                        "aimpt_coords": instrument_spec["aimpt_coords"]}

    if "chips" not in event_params:
        event_params["chips"] = instrument_spec["chips"]

    if foreground:
        mylog.info("Adding in astrophysical foreground.")
        bkg_events = make_foreground(event_params, arf, rmf, prng=prng)
        for key in bkg_events:
            events[key] = np.concatenate([events[key], bkg_events[key]])
    if instr_bkgnd and instrument_spec["bkgnd"] is not None:
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
                         ptsrc_bkgnd=True, no_dither=False, dither_params=None,
                         subpixel_res=False, input_sources=None, 
                         absorb_model="wabs", nH=0.05, prng=None):
    """
    Make an event file consisting entirely of background events. This will be 
    useful for creating backgrounds that can be added to simulations of sources.

    Parameters
    ----------
    exp_time : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
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
    instr_bkgnd : boolean, optional
        Whether or not to include the instrumental background. Default: True
    ptsrc_bkgnd : boolean, optional
        Whether or not to include the point-source background. Default: True
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
    input_sources : string, optional
        If set to a filename, input the point source positions, fluxes,
        and spectral indices from an ASCII table instead of generating
        them. Default: None
    absorb_model : string, optional
        The absorption model to use, "wabs" or "tbabs". Default: "wabs"
    nH : float, optional
        The hydrogen column in units of 10**22 atoms/cm**2. 
        Default: 0.05
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
                                           no_dither=no_dither,
                                           dither_params=dither_params, 
                                           subpixel_res=subpixel_res,
                                           input_sources=input_sources,
                                           absorb_model=absorb_model,
                                           nH=nH, prng=prng)
    write_event_file(events, event_params, out_file, overwrite=overwrite)

def instrument_simulator(input_events, out_file, exp_time, instrument,
                         sky_center, overwrite=False, instr_bkgnd=True, 
                         foreground=True, ptsrc_bkgnd=True, 
                         bkgnd_file=None, no_dither=False, 
                         dither_params=None, roll_angle=0.0, 
                         subpixel_res=False, prng=None):
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
    exp_time : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
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
    no_dither : boolean, optional
        If True, turn off dithering entirely. Default: False
    dither_params : array-like of floats, optional
        The parameters to use to control the size and period of the dither
        pattern. The first two numbers are the dither amplitude in x and y
        detector coordinates in arcseconds, and the second two numbers are
        the dither period in x and y detector coordinates in seconds. 
        Default: [8.0, 8.0, 1000.0, 707.0].
    roll_angle : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
        The roll angle of the observation in degrees. Default: 0.0
    subpixel_res: boolean, optional
        If True, event positions are not randomized within the pixels 
        within which they are detected. Default: False
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
    from soxs.background import add_background_from_file
    if not out_file.endswith(".fits"):
        out_file += ".fits"
    mylog.info("Making observation of source in %s." % out_file)
    # Make the source first
    events, event_params = generate_events(input_events, exp_time, instrument, sky_center,
                                           no_dither=no_dither, dither_params=dither_params, 
                                           roll_angle=roll_angle, subpixel_res=subpixel_res, 
                                           prng=prng)
    # If the user wants backgrounds, either make the background or add an already existing
    # background event file. It may be necessary to reproject events to a new coordinate system.
    if bkgnd_file is None:
        if not instr_bkgnd and not ptsrc_bkgnd and not foreground:
            mylog.info("No backgrounds will be added to this observation.")
        else:
            mylog.info("Adding background events.")
            bkg_events, _ = make_background(exp_time, instrument, sky_center,
                                            foreground=foreground, instr_bkgnd=instr_bkgnd, 
                                            no_dither=no_dither, dither_params=dither_params, 
                                            ptsrc_bkgnd=ptsrc_bkgnd, prng=prng, 
                                            subpixel_res=subpixel_res, roll_angle=roll_angle)
            for key in events:
                events[key] = np.concatenate([events[key], bkg_events[key]])
    else:
        mylog.info("Adding background events from the file %s." % bkgnd_file)
        if not os.path.exists(bkgnd_file):
            raise IOError("Cannot find the background event file %s!" % bkgnd_file)
        events = add_background_from_file(events, event_params, bkgnd_file)
    if len(events["energy"]) == 0:
        raise RuntimeError("No events were detected from source or background!!")
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
    exp_time : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
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
    exp_time = parse_value(exp_time, "s")
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
