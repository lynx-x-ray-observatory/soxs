import numpy as np
import os

from soxs.background.utils import BackgroundSpectrum
from soxs.simput import write_photon_list
from soxs.spatial import FillFOVModel
from soxs.utils import soxs_files_path

# X-ray foreground from Hickox & Markevitch 2007 (http://adsabs.harvard.edu/abs/2007ApJ...661L.117H)
hm_bkgnd_file = os.path.join(soxs_files_path, "hm_cxb_bkgnd.dat")
hm_astro_bkgnd = BackgroundSpectrum(hm_bkgnd_file, "astrophysical")

def make_foreground(simput_prefix, phlist_prefix, exp_time, fov, sky_center, 
                    area=40000.0, append=False, clobber=False, prng=np.random):
    r"""
    Generate a SIMPUT catalog corresponding to an astrophysical foreground. The model
    for the X-ray foreground is taken from Hickox & Markevitch 2007 
    (http://adsabs.harvard.edu/abs/2007ApJ...661L.117H).

    Parameters
    ----------
    simput_prefix : string
        The filename prefix for the SIMPUT file.
    phlist_prefix : string
        The filename prefix for the photon list file.
    exp_time : float
        The exposure time to use, in seconds. 
    fov : float
        The width of the field of view in arcminutes.
    sky_center : array, tuple, or list
        The center RA, Dec coordinates of the observation, in degrees.
    area : float, optional
        The effective area in cm**2. Used to tune the size of the sample of 
        photons that is generated.
    append : boolean, optional
        If True, append a new source an existing SIMPUT catalog. Default: False
    clobber : boolean, optional
        Set to True to overwrite previous files. Default: False
    prng : :class:`~numpy.random.RandomState` object or :mod:`~numpy.random`, optional
        A pseudo-random number generator. Typically will only be specified
        if you have a reason to generate the same set of random numbers, such as for a
        test. Default is the :mod:`numpy.random` module.
    """
    energy = hm_astro_bkgnd.generate_energies(exp_time, area, fov, prng=prng)
    fov_model = FillFOVModel(sky_center[0], sky_center[1], fov, energy.size)
    write_photon_list(simput_prefix, phlist_prefix, energy.flux, fov_model.ra, fov_model.dec,
                      energy, append=append, clobber=clobber)
