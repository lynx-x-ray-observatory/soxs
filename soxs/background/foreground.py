import os

from soxs.spectra import ApecGenerator
from soxs.background.spectra import BackgroundSpectrum, \
    ConvolvedBackgroundSpectrum
from soxs.background.events import make_diffuse_background
from soxs.utils import parse_prng, mylog, create_region, soxs_cfg
import numpy as np
from regions import PixCoord

"""
XSPEC model used to create the "default" foreground spectrum
  model  apec + wabs*apec
            0.099
                1
                0
          1.7e-06
            0.018
            0.225
                1
                0
          7.3e-07

XSPEC model used to create the "lem" foreground spectrum
  model  apec + wabs*(apec+apec)
            0.099
                1
                0
          1.7e-06
            0.018
            0.225
                1
                0
          7.3e-07
              0.7
                1
                0
         8.76e-08 
"""


class MakeFrgndSpectrum:
    def __init__(self):
        self._make_frgnd_spectrum()

    def __call__(self, apec_vers=None, abund_table=None):
        self._make_frgnd_spectrum(apec_vers=apec_vers, abund_table=abund_table)

    def _make_frgnd_spectrum(self, apec_vers=None, abund_table=None):
        bkgnd_nH = float(soxs_cfg.get("soxs", "bkgnd_nH"))
        absorb_model = soxs_cfg.get("soxs", "bkgnd_absorb_model")
        frgnd_spec_model = soxs_cfg.get("soxs", "frgnd_spec_model")
        agen = ApecGenerator(0.1, 10.0, 10000, apec_vers=apec_vers,
                             broadening=False, abund_table=abund_table)
        spec = agen.get_spectrum(0.225, 1.0, 0.0, 7.3e-7)
        if frgnd_spec_model == "halosat":
            spec += agen.get_spectrum(0.7, 1.0, 0.0, 8.76e-8)
        spec.apply_foreground_absorption(bkgnd_nH, model=absorb_model)
        spec += agen.get_spectrum(0.099, 1.0, 0.0, 1.7e-6)
        self.spec = BackgroundSpectrum.from_spectrum(spec, 1.0)


make_frgnd_spectrum = MakeFrgndSpectrum()


def make_foreground(event_params, arf, rmf, prng=None):

    prng = parse_prng(prng)

    conv_frgnd_spec = ConvolvedBackgroundSpectrum.convolve(make_frgnd_spectrum.spec, arf)

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

