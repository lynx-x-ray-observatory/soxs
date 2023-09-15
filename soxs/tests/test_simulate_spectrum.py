import os
import shutil
import tempfile

from astropy.io import fits
from numpy.testing import assert_allclose

from soxs.instrument import simulate_spectrum
from soxs.instrument_registry import get_instrument_from_registry
from soxs.spectra import Spectrum


def test_simulate_spectrum():
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    alpha_sim = 1.0
    nH_sim = 0.02
    norm_sim = 1.0e-4
    redshift = 0.01

    exp_time = (50.0, "ks")
    instr = "lynx_hdxi"

    spec = Spectrum.from_powerlaw(alpha_sim, redshift, norm_sim, 0.1, 10.0, 20000)
    spec.apply_foreground_absorption(nH_sim, model="tbabs")

    out_file1 = "test1.pha"
    out_file2 = "test2.pha"

    simulate_spectrum(
        spec,
        instr,
        exp_time,
        out_file1,
        instr_bkgnd=True,
        ptsrc_bkgnd=True,
        foreground=True,
        overwrite=True,
        bkgnd_area=(1.0, "arcmin**2"),
        prng=69,
    )

    instr_spec = get_instrument_from_registry(instr)

    simulate_spectrum(
        spec,
        (instr_spec["arf"], instr_spec["rmf"], instr_spec["bkgnd"]),
        exp_time,
        out_file2,
        instr_bkgnd=True,
        ptsrc_bkgnd=True,
        foreground=True,
        overwrite=True,
        bkgnd_area=(1.0, "arcmin**2"),
        prng=69,
    )

    f1 = fits.open(out_file1)
    spec1 = f1["SPECTRUM"].data["COUNTS"].copy()
    f2 = fits.open(out_file2)
    spec2 = f2["SPECTRUM"].data["COUNTS"].copy()

    assert_allclose(spec1, spec2)

    f1.close()
    f2.close()

    os.chdir(curdir)
    shutil.rmtree(tmpdir)
