import os

from soxs.background.spectra import BackgroundSpectrum, \
    ConvolvedBackgroundSpectrum
from soxs.background.events import make_diffuse_background
from soxs.utils import soxs_files_path, parse_prng, mylog
import numpy as np

# X-ray foreground from Hickox & Markevitch 2007 
# (http://adsabs.harvard.edu/abs/2007ApJ...661L.117H)
hm_bkgnd_file = os.path.join(soxs_files_path, "hm_cxb_bkgnd.h5")
hm_astro_bkgnd = BackgroundSpectrum.from_file(hm_bkgnd_file)

def make_foreground(event_params, arf, rmf, prng=None):
    import pyregion._region_filter as rfilter

    prng = parse_prng(prng)

    conv_frgnd_spec = ConvolvedBackgroundSpectrum(hm_astro_bkgnd, arf)

    energy = conv_frgnd_spec.generate_energies(event_params["exposure_time"],
                                               event_params["fov"], prng=prng, 
                                               quiet=True).value

    prng = parse_prng(prng)

    bkg_events = {}

    n_events = energy.size

    nx = event_params["num_pixels"]
    bkg_events["detx"] = prng.uniform(low=-0.5*nx, high=0.5*nx, size=n_events)
    bkg_events["dety"] = prng.uniform(low=-0.5*nx, high=0.5*nx, size=n_events)
    bkg_events["energy"] = energy

    if event_params["chips"] is None:
        bkg_events["chip_id"] = np.zeros(n_events, dtype='int')
    else:
        bkg_events["chip_id"] = -np.ones(n_events, dtype='int')
        for i, chip in enumerate(event_params["chips"]):
            thisc = np.ones(n_events, dtype='bool')
            rtype = chip[0]
            args = chip[1:]
            r = getattr(rfilter, rtype)(*args)
            inside = r.inside(bkg_events["detx"], bkg_events["dety"])
            thisc = np.logical_and(thisc, inside)
            bkg_events["chip_id"][thisc] = i

    keep = bkg_events["chip_id"] > -1

    if keep.sum() == 0:
        raise RuntimeError("No astrophysical foreground events were detected!!!")
    else:
        mylog.info("Making %d events from the astrophysical foreground." % keep.sum())

    for key in bkg_events:
        bkg_events[key] = bkg_events[key][keep]

    return make_diffuse_background(bkg_events, event_params, rmf, prng=prng)
