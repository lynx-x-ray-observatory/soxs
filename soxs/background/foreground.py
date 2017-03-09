import os

from soxs.background.spectra import BackgroundSpectrum, \
    ConvolvedBackgroundSpectrum
from soxs.background.events import make_uniform_background
from soxs.utils import soxs_files_path, parse_prng

# X-ray foreground from Hickox & Markevitch 2007 (http://adsabs.harvard.edu/abs/2007ApJ...661L.117H)
hm_bkgnd_file = os.path.join(soxs_files_path, "hm_cxb_bkgnd.dat")
hm_astro_bkgnd = BackgroundSpectrum(hm_bkgnd_file)

def make_foreground(event_params, arf, rmf, prng=None):
    prng = parse_prng(prng)
    fov = event_params["fov"]
    conv_bkgnd_spec = ConvolvedBackgroundSpectrum(hm_astro_bkgnd, arf)
    energy = conv_bkgnd_spec.generate_energies(event_params["exposure_time"], 
                                               fov, prng=prng).value
    if energy.size == 0:
        raise RuntimeError("No astrophysical foreground events were detected!!!")
    return make_uniform_background(energy, event_params, rmf, prng=prng)
