import os
import shutil
import sys
import tempfile

import numpy as np
import pytest
from astropy.io import fits
from numpy.random import RandomState
from numpy.testing import assert_allclose, assert_equal


@pytest.mark.skipif(sys.platform == "win32", reason="does not run on windows")
def test_pyxsim_read():
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    pyxsim = pytest.importorskip("pyxsim")

    import yt

    import soxs

    prng = RandomState(0x4D3D3D3)

    A = 2000.0
    exp_time = 1.0e4
    redshift = 0.1

    gslr = "GasSloshingLowRes/sloshing_low_res_hdf5_plt_cnt_0300"

    ds = yt.load_sample(gslr, default_species_fields="ionized")

    sphere = ds.sphere("c", (0.5, "Mpc"))

    thermal_model = pyxsim.CIESourceModel(
        "apec", 0.1, 11.0, 10000, 0.3, thermal_broad=False, prng=prng
    )

    _, _ = pyxsim.make_photons("photons1", sphere, redshift, A, exp_time, thermal_model)

    _ = pyxsim.project_photons(
        "photons1",
        "events1",
        [1.0, -0.5, 0.2],
        [30.0, 45.0],
        absorb_model="tbabs",
        nH=0.1,
        prng=prng,
    )

    events = pyxsim.EventList("events1.h5")
    events.write_to_simput("events1", overwrite=True)

    soxs.instrument_simulator(
        "events1_simput.fits",
        "evt1.fits",
        exp_time,
        "chandra_acisi_cy0",
        [30.0, 45.0],
        overwrite=True,
        ptsrc_bkgnd=False,
        foreground=False,
        instr_bkgnd=False,
        prng=29,
    )

    soxs.instrument_simulator(
        "events1.h5",
        "evt2.fits",
        exp_time,
        "chandra_acisi_cy0",
        [30.0, 45.0],
        overwrite=True,
        ptsrc_bkgnd=False,
        foreground=False,
        instr_bkgnd=False,
        prng=29,
    )

    with fits.open("evt1.fits") as e1, fits.open("evt2.fits") as e2:
        hdu1 = e1["EVENTS"]
        hdu2 = e2["EVENTS"]
        old_cols = hdu1.data.names
        new_cols = hdu2.data.names
        assert old_cols == new_cols
        for name in old_cols:
            dtype = hdu1.data[name].dtype
            if np.issubdtype(dtype, np.integer):
                assert_equal(hdu1.data[name], hdu2.data[name])
            else:
                if np.issubdtype(dtype, np.float32):
                    rtol = 1.0e-6
                else:
                    rtol = 1.0e-8
                assert_allclose(hdu1.data[name], hdu2.data[name], rtol=rtol)

    os.chdir(curdir)
    shutil.rmtree(tmpdir)
