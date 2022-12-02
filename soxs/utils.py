import os
import logging
import numpy as np
from numpy.random import RandomState
import astropy.units as u
from astropy.units import Quantity
import warnings
from configparser import ConfigParser
import regions
import appdirs
from scipy.interpolate import interp1d


# Configuration

soxs_cfg_defaults = {"soxs_data_dir": "/does/not/exist",
                     "abund_table": "angr",
                     "apec_vers": "3.0.9",
                     "spex_vers": "3.06.01",
                     "bkgnd_nH": 0.018,
                     "bkgnd_absorb_model": "tbabs",
                     "frgnd_spec_model": "default",
                     "frgnd_velocity": 0.0}

CONFIG_DIR = os.environ.get('XDG_CONFIG_HOME',
                            os.path.join(os.path.expanduser('~'),
                                         '.config', 'soxs'))
if not os.path.exists(CONFIG_DIR):
    try:
        os.makedirs(CONFIG_DIR)
    except OSError:
        warnings.warn("unable to create soxs config directory")

CURRENT_CONFIG_FILE = os.path.join(CONFIG_DIR, 'soxs.cfg')

if not os.path.exists(CURRENT_CONFIG_FILE):
    cp = ConfigParser()
    cp.add_section("soxs")
    try:
        with open(CURRENT_CONFIG_FILE, 'w') as new_cfg:
            cp.write(new_cfg)
    except IOError:
        warnings.warn("unable to write new config file")


soxs_cfg = ConfigParser(soxs_cfg_defaults)
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

if soxs_cfg.get("soxs", "soxs_data_dir") == "/does/not/exist":
    soxs_data_dir = appdirs.user_cache_dir("soxs")
    mylog.warning(f"Setting 'soxs_data_dir' to {soxs_data_dir} for this session. "
                  f"Please update your configuration if you want it somewhere else.")
    soxs_cfg.set("soxs", "soxs_data_dir", appdirs.user_cache_dir("soxs"))


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


def parse_value(value, default_units, equivalence=None):
    if isinstance(value, str):
        v = value.split(",")
        if len(v) == 2:
            value = (float(v[0]), v[1])
        else:
            value = float(v[0])
    if hasattr(value, "to_astropy"):
        value = value.to_astropy()
    if isinstance(value, Quantity):
        q = Quantity(value.value, value.unit)
    elif iterable(value):
        q = Quantity(value[0], value[1])
    else:
        q = Quantity(value, default_units)
    return q.to(default_units, equivalencies=equivalence).value


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


def line_width_equiv(rest):
    from astropy.constants import c
    ckms = c.to_value('km/s')
    forward = lambda x: rest*x/ckms
    backward = lambda x: x/rest*ckms
    return [(u.km/u.s, u.keV, forward, backward)]


class DummyPbar(object):
    def __init__(self):
        pass

    def update(self, *args, **kwargs):
        pass

    def close(self):
        pass


def create_region(rtype, args, dx, dy):
    if rtype in ["Rectangle", "Box"]:
        xctr, yctr, xw, yw = args
        x = xctr+dx
        y = yctr+dy
        center = regions.PixCoord(x=x, y=y)
        reg = regions.RectanglePixelRegion(center=center, width=xw, height=yw)
        bounds = [x-0.5*xw, x+0.5*xw, y-0.5*yw, y+0.5*yw]
    elif rtype == "Circle":
        xctr, yctr, radius = args
        x = xctr+dx
        y = yctr+dy
        center = regions.PixCoord(x=x, y=y)
        reg = regions.CirclePixelRegion(center=center, radius=radius)
        bounds = [x-radius, x+radius, y-radius, y+radius]
    elif rtype == "Polygon":
        x = np.array(args[0])+dx
        y = np.array(args[1])+dy
        vertices = regions.PixCoord(x=x, y=y)
        reg = regions.PolygonPixelRegion(vertices=vertices)
        bounds = [x.min(), x.max(), y.min(), y.max()]
    else:
        raise NotImplementedError
    return reg, bounds


def process_fits_string(fitsstr):
    import re
    from astropy.io import fits
    fn = fitsstr.split("[")[0]
    brackets = re.findall(r"[^[]*\[([^]]*)\]", fitsstr)
    with fits.open(fn) as f:
        if len(brackets) == 0:
            imgs = np.array([hdu.is_image for hdu in f])
            if imgs.sum() > 1:
                raise IOError("Multiple HDUs in this file, "
                              "please specify one to read!")
            ext = np.where(imgs)[0][0]
        else:
            ext = brackets[0]
            if ext.isdigit():
                ext = int(ext)
        imhdu = f[ext]
    return imhdu


class PoochHandle:
    r"""
    Container for a pooch object used to fetch remote response that isn't
    already stored locally.
    """
    def __init__(self, cache_dir=None):
        import json
        import pooch
        import pkg_resources
        if cache_dir is None:
            if os.path.isdir(soxs_cfg.get("soxs", "soxs_data_dir")):
                cache_dir = soxs_cfg.get("soxs", "soxs_data_dir")
            else:
                cache_dir = pooch.os_cache("soxs")
        self._registry = json.load(
            pkg_resources.resource_stream("soxs", "file_hash_registry.json"))
        self.pooch_obj = pooch.create(
            path=cache_dir,
            registry=self._registry,
            env="SOXS_DATA_DIR",
            base_url="https://hea-www.cfa.harvard.edu/soxs/soxs_responses/"
        )
        self.dl = pooch.HTTPDownloader(progressbar=True)

    def fetch(self, fname):
        return self.pooch_obj.fetch(fname, downloader=self.dl)


finley = PoochHandle()


def get_data_file(fn):
    soxs_data_dir = soxs_cfg.get("soxs", "soxs_data_dir")
    rel_fn = os.path.split(fn)[-1]
    data_fn = os.path.join(soxs_data_dir, rel_fn)
    if os.path.exists(fn):
        return fn
    elif os.path.exists(data_fn):
        return data_fn
    else:
        return finley.fetch(rel_fn)


def image_pos(im, nph, prng):
    im[im < 0.0] = 0.0
    im = im/im.sum()
    idxs = prng.choice(im.size, size=nph, p=im.flatten())
    x, y = np.unravel_index(idxs, im.shape)
    dx = prng.uniform(low=0.5, high=1.5, size=x.size)
    dy = prng.uniform(low=0.5, high=1.5, size=y.size)
    return x+dx, y+dy


def set_soxs_config(option, value):
    """
    Set SOXS configuration values.

    Parameters
    ----------
    option : string
        The option to change.
    value : number or string
        The value to set the option to.
    """
    soxs_cfg.set("soxs", option, value=str(value))


def set_mission_config(mission):
    """
    Set configuration options most appropriate for a specific
    *mission*. Currently only takes "lem".
    """
    if mission == "lem":
        frgnd_spec_model = "halosat"
        bkgnd_absorb_model = "tbabs"
        frgnd_velocity = 100.0 # km/s
        set_soxs_config("frgnd_spec_model", frgnd_spec_model)
        set_soxs_config("bkgnd_absorb_model", bkgnd_absorb_model)
        set_soxs_config("frgnd_velocity", frgnd_velocity)
    else:
        raise RuntimeError(f"Mission '{mission}' is not implemented!")


def regrid_spectrum(ebins_new, ebins, spec):
    cspec = np.insert(np.cumsum(spec, axis=-1), 0, 0.0, axis=-1)
    f = interp1d(ebins, cspec, axis=-1, fill_value=0.0,
                 assume_sorted=True, copy=False)
    return np.diff(f(ebins_new), axis=-1)
