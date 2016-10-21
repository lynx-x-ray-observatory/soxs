from __future__ import print_function

import json
import numpy as np
import astropy.io.fits as pyfits
from astropy.utils.console import ProgressBar
from soxs.constants import erg_per_keV
from soxs.utils import mylog, ensure_numpy_array

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
        earea = np.interp(energy, self.emid, self.eff_area, left=0.0, right=0.0)
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
                    nn = 0
                    for q in range(fcurr, last):
                        if low <= sorted_e[q] < high:
                            nn += 1
                        else:
                            break
                    channelInd = prng.choice(len(weights), size=nn, p=weights)
                    detectedChannels.append(trueChannel[channelInd])
                    fcurr += nn
                    pbar.update(fcurr)

        for key in events:
            events[key] = events[key][eidxs]
        events[self.header["CHANTYPE"]] = np.concatenate(detectedChannels)

        return events

instrument_registry = {}
instrument_registry["xcal"] = {"arf": "xrs_calorimeter.arf",
                               "rmf": "xrs_calorimeter.rmf",
                               "bkgnd": "acisi_particle_bkgnd.dat",
                               "num_pixels": 300,
                               "plate_scale": 1.0,
                               "psf": ["gaussian", 0.5]}
instrument_registry["hdxi"] = {"arf": "xrs_hdxi.arf",
                               "rmf": "xrs_hdxi.rmf",
                               "bkgnd": "acisi_particle_bkgnd.dat",
                               "num_pixels": 4096,
                               "plate_scale": 1./3.,
                               "psf": ["gaussian", 0.5]}

def add_instrument_to_registry(filename):
    """
    Add an instrument specification to the registry contained
    in a JSON file. The JSON file must have the structure as 
    shown below. The order is not important, but do not write 
    a JSON file with comments, the ones below are just for 
    clarity.

    >>> {
    ...     "name": "hdxi", # The short name of the instrument
    ...     "arf": "xrs_hdxi.arf", # The file containing the ARF
    ...     "rmf": "xrs_hdxi.rmf" # The file containing the RMF
    ...     "bkgnd": "acisi_particle_bkgnd.dat" # The file containing the particle background
    ...     "plate_scale": 0.33333333333, # The plate scale in arcsec
    ...     "num_pixels": 4096, # The number of pixels on a side in the FOV
    ...     "psf": [
    ...         "gaussian", 
    ...         0.5
    ...     ] # The type of PSF and its FWHM
    ... }

    Parameters
    ----------
    filename : string
        The JSON file containing the instrument specification.
    """
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