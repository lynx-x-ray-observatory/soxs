import os
import shutil
import tempfile

import numpy as np
from astropy.io import fits

from soxs.events import filter_events
from soxs.instrument import instrument_simulator
from soxs.simput import SimputCatalog, SimputSpectrum
from soxs.spectra import Spectrum


def test_filter():

    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    alpha_sim = 1.1
    nH_sim = 0.02
    norm_sim = 1.0e-4
    redshift = 0.01

    exp_time = (100.0, "ks")
    inst_name = "chandra_acisi_cy0"

    spec = Spectrum.from_powerlaw(alpha_sim, redshift, norm_sim, 0.1, 10.0, 2000)
    spec.apply_foreground_absorption(nH_sim, model="tbabs")

    pt_src = SimputSpectrum.from_spectrum("plaw_model", spec, 30.0, 45.0)
    SimputCatalog.from_source("plaw_model_simput.fits", pt_src, overwrite=True)

    instrument_simulator(
        "plaw_model_simput.fits",
        "evt.fits",
        exp_time,
        inst_name,
        [30.0, 45.0],
        instr_bkgnd=False,
        ptsrc_bkgnd=False,
        foreground=True,
        prng=24,
    )

    filter_events("evt.fits", "evt_filter_en.fits", emin=0.5, emax=2.0, overwrite=True)
    with fits.open("evt_filter_en.fits") as f:
        e = f["EVENTS"].data["ENERGY"].copy()
        assert np.logical_and(e > 500.0, e < 2000.0).all()

    reg = "# Region file format: DS9\nimage\ncircle(2430,2454,200)\n"
    filter_events("evt.fits", "evt_filter_reg.fits", region=reg, overwrite=True)
    with fits.open("evt_filter_reg.fits") as f:
        x = f["EVENTS"].data["X"].copy()
        y = f["EVENTS"].data["Y"].copy()
        r = np.sqrt((x - 2430) ** 2 + (y - 2454) ** 2)
        assert np.all(r < 201.0)

    os.chdir(curdir)
    shutil.rmtree(tmpdir)
