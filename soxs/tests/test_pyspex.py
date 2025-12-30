import numpy as np
import pytest

from soxs import Spectrum

pyspex = pytest.importorskip("pyspex", reason="pyspex not installed")


def test_pyspex():
    # For now, this is just a test to make sure it works

    s = pyspex.spex.Session()
    s.egrid(0.2, 10.0, 5000, "kev", False)
    s.abun("aspl")

    s.com("reds")
    s.com("hot")
    s.com("cie")
    s.com_rel(1, 3, np.array([1, 2]))

    s.dist(1, 0.05, "z")  # Set the distance in SPEX to z=0.05 (Luminosity).
    s.par(1, 1, "z", 0.05, thawn=True)  # Set the redshift in the model to z=0.05 (Energy shift in spectrum).
    s.par(1, 2, "nh", 2e-4)  # Set the hydrogen column density.
    s.par_fix(1, 2, "t")  # Fix the temperature of the galactic absorption to the default value
    s.par(1, 3, "norm", 1e8, thawn=True)  # Set a guess for the normalisation of the cluster component.
    s.par(1, 3, "t", 4.0, thawn=True)  # Guess that the temperature is 4.0 keV.

    s.calc()

    spec = Spectrum.from_spex_model(s)  # noqa
