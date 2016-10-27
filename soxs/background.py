from __future__ import print_function
import numpy as np
import os
from soxs.spectra import Spectrum
from soxs.utils import soxs_files_path

class BackgroundSpectrum(Spectrum):
    _units = "photon/(cm**2*s*keV*arcmin**2)"
    def __init__(self, filename, bkgnd_type):
        self.bkgnd_type = bkgnd_type
        emid, flux = np.loadtxt(filename, unpack=True)
        de = np.diff(emid)[0]
        ebins = np.append(emid-0.5*de, emid[-1]+0.5*de)
        super(BackgroundSpectrum, self).__init__(ebins, flux)

    def generate_energies(self, t_exp, area, fov, bkgnd_scale, prng=None):
        """
        Generate photon energies from this background spectrum given an
        exposure time, effective area, and field of view.

        Parameters
        ----------
        t_exp : float
            The exposure time in seconds.
        area : float or NumPy array
            The effective area in cm**2. If one is creating events for a SIMPUT file,
            a constant should be used and it must be large enough so that a sufficiently
            large sample is drawn for the ARF.
        fov : float
            The width of the field of view on a side in arcminutes.
        bkgnd_scale : float
            A uniform scaling factor for the background. 
        prng : :class:`~numpy.random.RandomState` object or :mod:`~numpy.random`, optional
            A pseudo-random number generator. Typically will only be specified
            if you have a reason to generate the same set of random numbers, such as for a
            test. Default is the :mod:`numpy.random` module.
        """
        A = area*fov*fov*bkgnd_scale
        return super(BackgroundSpectrum, self).generate_energies(t_exp, A,
                                                                 prng=prng)

# ACIS-I particle background
acisi_bkgnd_file = os.path.join(soxs_files_path, "acisi_particle_bkgnd.dat")
acisi_particle_bkgnd = BackgroundSpectrum(acisi_bkgnd_file, "instrumental")

# X-ray foreground from Hickox & Markevitch 2007 (http://adsabs.harvard.edu/abs/2007ApJ...661L.117H)
hm_bkgnd_file = os.path.join(soxs_files_path, "hm_cxb_bkgnd.dat")
hm_astro_bkgnd = BackgroundSpectrum(hm_bkgnd_file, "astrophysical")

background_registry = {"acisi": acisi_particle_bkgnd,
                       "hm_cxb": hm_astro_bkgnd}

def add_background_to_registry(name, filename, bkgnd_type):
    """
    Add a background to the background registry.

    Parameters
    ----------
    name : string
        The short name of the background, which will be the key in the 
        registry.
    filename : string
        The file containing the background. It must have two columns: 
        energy in keV, and background intensity in units of
        photons/s/cm**2/arcmin**2/keV.
    bkgnd_type : string
        The type of background, either "instrumental" or "astrophysical".
    """
    background_registry[name] = BackgroundSpectrum(filename, bkgnd_type)

def show_background_registry():
    """
    Print the contents of the background registry.
    """
    for name, spec in background_registry.items():
        print("Background: %s" % name)
        print("    Type: %s" % spec.bkgnd_type)
        print("    Total Flux (%s - %s): %s" % (spec.ebins[0], spec.ebins[-1], spec.total_flux))
