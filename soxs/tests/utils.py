import os
from numpy.testing import assert_equal
from astropy.io import fits
from soxs.utils import soxs_cfg
import shutil


def get_test_path():
    test_path = soxs_cfg.get("soxs", "test_outputs_path")
    if not os.path.exists(test_path):
        raise IOError("The SOXS test outputs directory "
                      "%s does not exist!" % test_path)
    return test_path


def spectrum_answer_testing(spec, filename, generate):
    test_path = get_test_path()
    testfile = os.path.join(test_path, filename)
    if generate:
        spec.write_h5_file(testfile, overwrite=True)
    else:
        answer_spec = type(spec).from_file(testfile)
        assert_equal(answer_spec.emid.value, spec.emid.value)
        assert_equal(answer_spec.flux.value, spec.flux.value)
        assert answer_spec.flux.unit == spec.flux.unit


def file_answer_testing(hdu, filename, generate):
    test_path = get_test_path()
    oldf = os.path.join(test_path, filename)
    if generate:
        shutil.copy(filename, test_path)
    else:
        f_old = fits.open(oldf)
        f_new = fits.open(filename)
        old_cols = f_old[hdu].data.names
        new_cols = f_new[hdu].data.names
        assert old_cols == new_cols
        for name in old_cols:
            assert_equal(f_old[hdu].data[name], f_old[hdu].data[name])
        f_old.close()
        f_new.close()
