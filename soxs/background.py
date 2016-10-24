import numpy as np
import os
from soxs.spectra import Spectrum
from soxs.utils import soxs_files_path

class BackgroundSpectrum(Spectrum):
    def __init__(self, filename, bkgnd_type):
        self.bkgnd_type = bkgnd_type
        emid, flux = np.loadtxt(filename, unpack=True)
        de = np.diff(emid)[0]
        ebins = np.append(emid-0.5*de, emid[-1]+0.5*de)
        super(BackgroundSpectrum, self).__init__(ebins, flux)
        self.units = "photons/cm**2/arcmin**2/s/keV"

    def generate_energies(self, t_exp, area, fov, bkgnd_scale, prng=None):
        A = area*fov*fov*bkgnd_scale
        return super(BackgroundSpectrum, self).generate_energies(t_exp, A,
                                                                 prng=prng)

    def __repr__(self):
        s = "BackgroundSpectrum (%g - %g keV): " % (self.ebins[0], self.ebins[-1])
        s += "Total flux %g (%g) photons (erg) / cm**2 / arcmin**2 / s" % (self.tot_flux,
                                                                           self.tot_energy_flux)
        return s

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