import os
import logging
import astropy.io.fits as pyfits
import numpy as np
from copy import copy

soxsLogger = logging.getLogger("soxs")

ufstring = "%(name)-3s : [%(levelname)-9s] %(asctime)s %(message)s"
cfstring = "%(name)-3s : [%(levelname)-18s] %(asctime)s %(message)s"

soxs_sh = logging.StreamHandler()
# create formatter and add it to the handlers
formatter = logging.Formatter(ufstring)
soxs_sh.setFormatter(formatter)
# add the handler to the logger
soxsLogger.addHandler(soxs_sh)
soxsLogger.setLevel('INFO')
soxsLogger.propagate = False

mylog = soxsLogger

mylog.setLevel('INFO')

soxs_path = os.path.abspath(os.path.dirname(__file__))
soxs_files_path = os.path.join(soxs_path, "files")

def check_file_location(fn, subdir):
    if os.path.exists(fn):
        return os.path.abspath(fn)
    else:
        sto_fn = os.path.join(soxs_path, subdir, fn)
        if os.path.exists(sto_fn):
            return sto_fn
    raise IOError("Could not find file %s!" % fn)

def iterable(obj):
    """
    Grabbed from Python Cookbook / matploblib.cbook.  Returns true/false for
    *obj* iterable.
    """
    try: len(obj)
    except: return False
    return True

def ensure_list(obj):
    """
    This function ensures that *obj* is a list.  Typically used to convert a
    string to a list, for instance ensuring the *fields* as an argument is a
    list.
    """
    if obj is None:
        return [obj]
    if not isinstance(obj, list):
        return [obj]
    return obj

def ensure_numpy_array(obj):
    """
    This function ensures that *obj* is a numpy array. Typically used to
    convert scalar, list or tuple argument passed to functions using Cython.
    """
    if isinstance(obj, np.ndarray):
        if obj.shape == ():
            return np.array([obj])
        # We cast to ndarray to catch ndarray subclasses
        return np.array(obj)
    elif isinstance(obj, (list, tuple)):
        return np.asarray(obj)
    else:
        return np.asarray([obj])

def convert_rmf(rmffile):

    f = pyfits.open(rmffile)

    names = [ff.name for ff in f]
    idx = -1
    for i, name in enumerate(names):
        if "MATRIX" in name:
            idx = i
            break

    new_f = copy(f)

    matrix = new_f[names[idx]]
    fchan = matrix.data["F_CHAN"]
    nchan = matrix.data["N_CHAN"]
    m = matrix.data["MATRIX"]

    fchan_new = pyfits.Column(name = 'F_CHAN', format='PJ()',
                              array=np.array([fc[fc > 0] for fc in fchan], dtype=np.object))
    nchan_new = pyfits.Column(name = 'N_CHAN', format='PJ()',
                              array=np.array([nc[nc > 0] for nc in nchan], dtype=np.object))
    m_new = pyfits.Column(name = 'MATRIX', format='PE()',
                          array=np.array([mm[mm > 0] for mm in m], dtype=np.object))

    matrix_new = pyfits.BinTableHDU.from_columns([matrix.columns["ENERG_LO"],
                                                  matrix.columns["ENERG_HI"],
                                                  matrix.columns["N_GRP"],
                                                  fchan_new, nchan_new,
                                                  m_new])

    matrix_new.name = "SPECRESP MATRIX"

    new_f.pop(idx)

    new_f.append(matrix_new)

    header_keys = ["HDUCLASS", "HDUCLAS1", "HDUVERS1", "HDUCLAS2", "HDUVERS2", "HDUCLAS3",
                   "TELESCOP", "INSTRUME", "DETNAM", "FILTER", "DETCHANS", "CHANTYPE",
                   "HIERARCH LO_THRESH", "LO_THRES", "RMFVERSN", "TLMIN4", "TLMAX4"]

    for key in header_keys:
        if key in matrix.header:
            matrix_new.header[key] = matrix.header[key]

    new_f.writeto(os.path.split(rmffile)[-1], clobber=True)
