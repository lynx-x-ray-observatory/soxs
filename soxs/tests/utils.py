import os
from numpy.testing import assert_equal
from astropy.io import fits
from soxs.utils import soxs_cfg

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
        spec.write_file(testfile, overwrite=True)
    else:
        answer_spec = type(spec).from_file(testfile)
        assert_equal(answer_spec.emid.value, spec.emid.value)
        assert_equal(answer_spec.flux.value, spec.flux.value)
        assert answer_spec.flux.units == spec.flux.units

def file_answer_testing(file_type, old_file, new_file, generate):
    test_path = get_test_path()
    oldf = os.path.join(test_path, old_file)
    if generate:
        
    else:
        f_old = fits.open(oldf)
        f_new = fits.open(new_file)
        
        f_old.close()
        f_new.close()
