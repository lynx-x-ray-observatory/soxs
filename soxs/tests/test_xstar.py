from numpy.testing import assert_allclose
from pathlib import Path

from soxs.spectra.base import CountRateSpectrum
from soxs.utils import soxs_cfg


def test_xstar():
    answer_dir = Path(soxs_cfg.get("soxs", "soxs_answer_dir"))
    spec = CountRateSpectrum.from_xstar_model(answer_dir / "xout_spect1.fits", "emit_outward")
    rate, lum = spec.get_lum_in_band(0.5, 7.0)
    assert_allclose(rate.value, 4.04361905e54)
    assert_allclose(lum.value, 3.77093392e45)
