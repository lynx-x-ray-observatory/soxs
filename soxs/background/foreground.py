import os

from soxs.background.spectra import BackgroundSpectrum, \
    ConvolvedBackgroundSpectrum
from soxs.background.events import make_diffuse_background
from soxs.utils import soxs_files_path, parse_prng, mylog, \
    create_region
import numpy as np
from regions import PixCoord

# X-ray foreground from Hickox & Markevitch 2007
# (http://adsabs.harvard.edu/abs/2007ApJ...661L.117H)
hm_bkgnd_file = os.path.join(soxs_files_path, "hm_cxb_bkgnd.h5")
hm_astro_bkgnd = BackgroundSpectrum.from_file(hm_bkgnd_file)


def make_foreground(event_params, arf, rmf, prng=None):

    prng = parse_prng(prng)

    conv_frgnd_spec = ConvolvedBackgroundSpectrum.convolve(hm_astro_bkgnd, arf)

    bkg_events = {"energy": [], "detx": [], "dety": [], "chip_id": []}
    pixel_area = (event_params["plate_scale"]*60.0)**2
    for i, chip in enumerate(event_params["chips"]):
        rtype = chip[0]
        args = chip[1:]
        r, bounds = create_region(rtype, args, 0.0, 0.0)
        fov = np.sqrt((bounds[1]-bounds[0])*(bounds[3]-bounds[2])*pixel_area)
        e = conv_frgnd_spec.generate_energies(event_params["exposure_time"],
                                              fov, prng=prng, quiet=True).value
        n_events = e.size
        detx = prng.uniform(low=bounds[0], high=bounds[1], size=n_events)
        dety = prng.uniform(low=bounds[2], high=bounds[3], size=n_events)
        if rtype in ["Box", "Rectangle"]:
            thisc = slice(None, None, None)
            n_det = n_events
        else:
            thisc = r.contains(PixCoord(detx, dety))
            n_det = thisc.sum()
        bkg_events["energy"].append(e[thisc])
        bkg_events["detx"].append(detx[thisc])
        bkg_events["dety"].append(dety[thisc])
        bkg_events["chip_id"].append(i*np.ones(n_det))

    for key in bkg_events:
        bkg_events[key] = np.concatenate(bkg_events[key])

    if bkg_events["energy"].size == 0:
        raise RuntimeError("No astrophysical foreground events "
                           "were detected!!!")
    else:
        mylog.info(f"Making {bkg_events['energy'].size} events from the "
                   f"astrophysical foreground.")

    bkg_events = make_diffuse_background(bkg_events, 
                                         event_params, rmf, prng=prng)
    mylog.info(f"Scattering energies with "
               f"RMF {os.path.split(rmf.filename)[-1]}.")

    return rmf.scatter_energies(bkg_events, prng=prng)

