import os
from numpy.testing import assert_almost_equal
from astropy.io import fits
import shutil


def spectrum_answer_testing(spec, filename, answer_store, answer_dir):
    testfile = os.path.join(answer_dir, filename)
    if answer_store:
        spec.write_h5_file(testfile, overwrite=True)
    else:
        answer_spec = type(spec).from_file(testfile)
        assert_almost_equal(answer_spec.emid.value,
                             spec.emid.value)
        assert_almost_equal(answer_spec.flux.value,
                             spec.flux.value)
        assert answer_spec.flux.unit == spec.flux.unit


def file_answer_testing(hdu, filename, answer_store, answer_dir):
    oldf = os.path.join(answer_dir, filename)
    if answer_store:
        shutil.copy(filename, answer_dir)
    else:
        f_old = fits.open(oldf)
        f_new = fits.open(filename)
        old_cols = f_old[hdu].data.names
        new_cols = f_new[hdu].data.names
        assert old_cols == new_cols
        for name in old_cols:
            assert_almost_equal(f_old[hdu].data[name], f_new[hdu].data[name])
        f_old.close()
        f_new.close()
