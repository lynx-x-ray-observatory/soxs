import os
import shutil

import numpy as np
import pytest
from astropy.io import fits
from numpy.testing import assert_allclose, assert_equal

from soxs.utils import soxs_cfg

answer_dir = soxs_cfg.get("soxs", "soxs_answer_dir")


def spectrum_answer_testing(spec, filename, answer_store, rtol=1.0e-7):
    testfile = os.path.join(answer_dir, filename)
    if answer_store:
        spec.write_hdf5_file(testfile, overwrite=True)
    else:
        answer_spec = type(spec).from_file(testfile)
        assert_allclose(answer_spec.emid.value, spec.emid.value, rtol=rtol)
        assert_allclose(answer_spec.flux.value, spec.flux.value, rtol=rtol)
        assert answer_spec.flux.unit == spec.flux.unit


def file_answer_testing(hdu, filename, answer_store):
    oldf = os.path.join(answer_dir, filename)
    if answer_store:
        shutil.copy(filename, answer_dir)
    else:
        f_old = fits.open(oldf)
        f_new = fits.open(filename)
        if hdu in ["IMAGE", "EXPMAP"]:
            for k in f_old[hdu].header:
                assert_equal(f_old[hdu].header[k], f_new[hdu].header[k])
            dtype = f_old[hdu].data.dtype
            if np.issubdtype(dtype, np.float32):
                rtol = 1.0e-6
            else:
                rtol = 1.0e-8
            assert_allclose(f_old[hdu].data, f_new[hdu].data, rtol=rtol)
        else:
            old_cols = f_old[hdu].data.names
            new_cols = f_new[hdu].data.names
            assert old_cols == new_cols
            for name in old_cols:
                dtype = f_old[hdu].data[name].dtype
                if np.issubdtype(dtype, np.integer):
                    assert_equal(f_old[hdu].data[name], f_new[hdu].data[name])
                else:
                    if np.issubdtype(dtype, np.float32):
                        rtol = 1.0e-6
                    else:
                        rtol = 1.0e-8
                    assert_allclose(
                        f_old[hdu].data[name], f_new[hdu].data[name], rtol=rtol
                    )
        f_old.close()
        f_new.close()


min_numpy_vers = pytest.mark.skipif(
    np.lib.NumpyVersion(np.__version__) < np.lib.NumpyVersion("2.0.0"),
    reason="Requires NumPy >= 2.0.0",
)
