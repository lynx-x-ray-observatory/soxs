import os
import logging
import astropy.io.fits as pyfits
import numpy as np
from copy import copy
from numpy.random import RandomState
from astropy.units import Quantity
from six import string_types
from six.moves import configparser
import warnings

# Configuration

soxs_cfg_defaults = {"response_path": "/does/not/exist"}

CONFIG_DIR = os.environ.get(
    'XDG_CONFIG_HOME', os.path.join(os.path.expanduser('~'), '.config', 'soxs'))
if not os.path.exists(CONFIG_DIR):
    try:
        os.makedirs(CONFIG_DIR)
    except OSError:
        warnings.warn("unable to create soxs config directory")

CURRENT_CONFIG_FILE = os.path.join(CONFIG_DIR, 'soxs.cfg')

if not os.path.exists(CURRENT_CONFIG_FILE):
    cp = configparser.ConfigParser()
    cp.add_section("soxs")
    try:
        with open(CURRENT_CONFIG_FILE, 'w') as new_cfg:
            cp.write(new_cfg)
    except IOError:
        warnings.warn("unable to write new config file")

soxs_cfg = configparser.ConfigParser(soxs_cfg_defaults)
soxs_cfg.read([CURRENT_CONFIG_FILE, 'soxs.cfg'])
if not soxs_cfg.has_section("soxs"):
    soxs_cfg.add_section("soxs")

# Logging

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

def issue_deprecation_warning(msg):
    import warnings
    from numpy import VisibleDeprecationWarning
    warnings.warn(msg, VisibleDeprecationWarning, stacklevel=3)

soxs_path = os.path.abspath(os.path.dirname(__file__))
soxs_files_path = os.path.join(soxs_path, "files")

def parse_prng(prng):
    if isinstance(prng, RandomState):
        return prng
    else:
        return RandomState(prng)

def iterable(obj):
    """
    Grabbed from Python Cookbook / matplotlib.cbook.
    Returns true/false for *obj* iterable.
    """
    try:
        len(obj)
    except:
        return False
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
    This function ensures that *obj* is a numpy array. 
    Typically used to convert scalar, list or tuple 
    argument passed to functions using Cython.
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

    new_f.writeto(os.path.split(rmffile)[-1], overwrite=True)

def parse_value(value, default_units):
    if isinstance(value, string_types):
        v = value.split(",")
        if len(v) == 2:
            value = (float(v[0]), v[1])
        else:
            value = float(v[0])
    if hasattr(value, "to_astropy"):
        value = value.to_astropy()
    if isinstance(value, Quantity):
        q = Quantity(value.value, value.unit).to(default_units)
    elif iterable(value):
        q = Quantity(value[0], value[1]).to(default_units)
    else:
        q = Quantity(value, default_units)
    return q.value

def get_rot_mat(roll_angle):
    roll_angle = np.deg2rad(roll_angle)
    rot_mat = np.array([[np.cos(roll_angle), -np.sin(roll_angle)],
                        [np.sin(roll_angle), np.cos(roll_angle)]])
    return rot_mat

def downsample(myarr,factor,estimator=np.mean):
    """
    Downsample a 2D array by averaging over *factor* pixels in each axis.
    Crops upper edge if the shape is not a multiple of factor.

    This code is pure numpy and should be fast.

    keywords:
        estimator - default to mean.  You can downsample by summing or
            something else if you want a different estimator
            (e.g., downsampling error: you want to sum & divide by sqrt(n))
    """
    ys,xs = myarr.shape
    crarr = myarr[:ys-(ys % int(factor)),:xs-(xs % int(factor))]
    dsarr = estimator(np.concatenate([[crarr[i::factor,j::factor]
                                       for i in range(factor)]
                                      for j in range(factor)]), axis=0)
    return dsarr
