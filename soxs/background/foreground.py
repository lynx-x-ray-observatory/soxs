import os

from soxs.background.spectra import BackgroundSpectrum, \
    ConvolvedBackgroundSpectrum
from soxs.background.events import make_uniform_background
from soxs.utils import soxs_files_path, parse_prng, mylog

# X-ray foreground from Hickox & Markevitch 2007 (http://adsabs.harvard.edu/abs/2007ApJ...661L.117H)
hm_bkgnd_file = os.path.join(soxs_files_path, "hm_cxb_bkgnd.h5")
hm_astro_bkgnd = BackgroundSpectrum.from_file(hm_bkgnd_file)

def make_foreground(event_params, arf, rmf, prng=None):
    prng = parse_prng(prng)

    conv_bkgnd_spec = ConvolvedBackgroundSpectrum(hm_astro_bkgnd, arf)

    energy = conv_bkgnd_spec.generate_energies(event_params["exposure_time"],
                                               event_params["fov"], prng=prng, 
                                               quiet=True).value

    if energy.size == 0:
        raise RuntimeError("No astrophysical foreground events were detected!!!")
    else:
        mylog.info("Making %d events from the astrophysical foreground." % energy.size)

    return make_uniform_background(energy, event_params, rmf, prng=prng)
