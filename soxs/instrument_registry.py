import json
from soxs.utils import mylog, parse_value, finley
import os
from copy import deepcopy

# The Instrument Registry


class InstrumentRegistry:
    def __init__(self):
        self.registry = {}

    def __getitem__(self, key):
        return self.registry[key]

    def __setitem__(self, key, value):
        self.registry[key] = value

    def keys(self):
        return self.registry.keys()

    def items(self):
        return self.registry.items()

    def __contains__(self, item):
        return item in self.registry

    def get(self, key, default=None):
        return self.registry.get(key, default)

    def set(self, key, value, default=None):
        self.registry.set(key, value, default)


instrument_registry = InstrumentRegistry()

## Lynx

# High-Definition X-ray Imager (HDXI)

instrument_registry["lynx_hdxi"] = {"name": "lynx_hdxi",
                                    "arf": "xrs_hdxi_3x10.arf",
                                    "rmf": "xrs_hdxi.rmf",
                                    "bkgnd": "acisi",
                                    "fov": 22.0,
                                    "num_pixels": 4096,
                                    "aimpt_coords": [0.0, 0.0],
                                    "chips": [["Box", 0, 0, 4096, 4096]],
                                    "focal_length": 10.0,
                                    "dither": True,
                                    "psf": ["gaussian", 0.5],
                                    "imaging": True,
                                    "grating": False}

# Micro-calorimeter

instrument_registry["lynx_lxm"] = {"name": "lynx_lxm",
                                   "arf": "xrs_mucal_3x10_3.0eV.arf",
                                   "rmf": "xrs_mucal_3.0eV.rmf",
                                   "bkgnd": "mucal",
                                   "fov": 5.0,
                                   "num_pixels": 300,
                                   "aimpt_coords": [0.0, 0.0],
                                   "chips": [["Box", 0, 0, 300, 300]],
                                   "focal_length": 10.0,
                                   "dither": True,
                                   "psf": ["gaussian", 0.5],
                                   "imaging": True,
                                   "grating": False}

instrument_registry["lynx_lxm_enh"] = {"name": "lynx_lxm_enh",
                                       "arf": "xrs_mucal_3x10_1.5eV.arf",
                                       "rmf": "xrs_mucal_1.5eV.rmf",
                                       "bkgnd": "mucal",
                                       "fov": 1.0,
                                       "num_pixels": 120,
                                       "aimpt_coords": [0.0, 0.0], 
                                       "chips": [["Box", 0, 0, 120, 120]],
                                       "focal_length": 10.0,
                                       "dither": True,
                                       "psf": ["gaussian", 0.5],
                                       "imaging": True,
                                       "grating": False}

instrument_registry["lynx_lxm_ultra"] = {"name": "lynx_lxm_ultra",
                                         "arf": "xrs_mucal_3x10_0.3eV.arf",
                                         "rmf": "xrs_mucal_0.3eV.rmf",
                                         "bkgnd": "mucal",
                                         "fov": 1.0,
                                         "num_pixels": 60,
                                         "aimpt_coords": [0.0, 0.0],
                                         "chips": [["Box", 0, 0, 60, 60]],
                                         "focal_length": 10.0,
                                         "dither": True,
                                         "psf": ["gaussian", 0.5],
                                         "imaging": True,
                                         "grating": False}


# Gratings (for spectra only)

instrument_registry["lynx_xgs"] = {"name": "lynx_xgs",
                                   "arf": "xrs_cat.arf",
                                   "rmf": "xrs_cat.rmf",
                                   "bkgnd": None,
                                   "focal_length": 10.0,
                                   "imaging": False,
                                   "grating": True}

## Athena

# WFI

instrument_registry["athena_wfi"] = {"name": "athena_wfi",
                                     "arf": "athena_wfi_rib2.3_B4C_20190122_wo_filter_OnAxis.arf",
                                     "rmf": "athena_wfi_rib2.3_B4C_20190122_wo_filter_OnAxis.rmf",
                                     "bkgnd": "athena_wfi",
                                     "fov": 40.147153,
                                     "num_pixels": 1078,
                                     "aimpt_coords": [53.69, -53.69],
                                     "chips": [["Box", -283, -283, 512, 512],
                                               ["Box", 283, -283, 512, 512],
                                               ["Box", -283, 283, 512, 512],
                                               ["Box", 283, 283, 512, 512]],
                                     "focal_length": 12.0,
                                     "dither": False,
                                     "psf": ["gaussian", 5.0],
                                     "imaging": True,
                                     "grating": False}

# XIFU

instrument_registry["athena_xifu"] = {"name": "athena_xifu",
                                      "arf": "XIFU_CC_BASELINECONF_2018_10_10.arf",
                                      "rmf": "XIFU_CC_BASELINECONF_2018_10_10.rmf",
                                      "bkgnd": "athena_xifu",
                                      "fov": 5.991992621478149,
                                      "num_pixels": 84,
                                      "aimpt_coords": [0.0, 0.0],
                                      "chips": [["Polygon", 
                                                 [-33, 0, 33, 33, 0, -33],
                                                 [20, 38, 20, -20, -38, -20]]],
                                      "focal_length": 12.0,
                                      "dither": False,
                                      "psf": ["gaussian", 5.0],
                                      "imaging": True, 
                                      "grating": False}

## Chandra

# ACIS-I, Cycle 0 and 20

for cycle in [0, 22]:
    name = f"chandra_acisi_cy{cycle}"
    instrument_registry[name] = {"name": name, 
                                 "arf": f"acisi_aimpt_cy{cycle}.arf",
                                 "rmf": f"acisi_aimpt_cy{cycle}.rmf",
                                 "bkgnd": "acisi",
                                 "fov": 20.008,
                                 "num_pixels": 2440,
                                 "aimpt_coords": [86.0, 57.0],
                                 "chips": [["Box", -523, -523, 1024, 1024],
                                           ["Box", 523, -523, 1024, 1024],
                                           ["Box", -523, 523, 1024, 1024],
                                           ["Box", 523, 523, 1024, 1024]],
                                 "psf": ["gaussian", 0.5],
                                 "focal_length": 10.0,
                                 "dither": True,
                                 "imaging": True,
                                 "grating": False}

# ACIS-S, Cycle 0 and 22

for cycle in [0, 22]:
    name = f"chandra_aciss_cy{cycle}"
    instrument_registry[name] = {"name": name,
                                 "arf": f"aciss_aimpt_cy{cycle}.arf",
                                 "rmf": f"aciss_aimpt_cy{cycle}.rmf",
                                 "bkgnd": ["acisi", "aciss",
                                           "acisi", "aciss",
                                           "acisi", "acisi"],
                                 "fov": 50.02,
                                 "num_pixels": 6100,
                                 "aimpt_coords": [206.0, 0.0],
                                 "chips": [["Box", -2605, 0, 1024, 1024],
                                           ["Box", -1563, 0, 1024, 1024],
                                           ["Box", -521, 0, 1024, 1024],
                                           ["Box", 521, 0, 1024, 1024],
                                           ["Box", 1563, 0, 1024, 1024],
                                           ["Box", 2605, 0, 1024, 1024]],
                                 "psf": ["gaussian", 0.5],
                                 "focal_length": 10.0,
                                 "dither": True,
                                 "imaging": True, 
                                 "grating": False}


# ACIS-S, Cycle 0 and 19 HETG

orders = {"p1": 1, "m1": -1}

for energy in ["meg", "heg"]:
    for order in ["p1", "m1"]:
        for cycle in [0, 22]:
            name = f"chandra_aciss_{energy}_{order}_cy{cycle}"
            resp_name = f"chandra_aciss_{energy}{orders[order]}_cy{cycle}"
            instrument_registry[name] = {"name": name,
                                         "arf": "%s.garf" % resp_name,
                                         "rmf": "%s.grmf" % resp_name,
                                         "bkgnd": None,
                                         "focal_length": 10.0,
                                         "imaging": False,
                                         "grating": True}

## Hitomi

# SXS

instrument_registry["xrism_resolve"] = {"name": "xrism_resolve",
                                        "arf": "xarm_res_flt_pa_20170818.arf",
                                        "rmf": "xarm_res_h5ev_20170818.rmf",
                                        "bkgnd": "hitomi_sxs",
                                        "num_pixels": 6,
                                        "fov": 3.06450576,
                                        "aimpt_coords": [0.0, 0.0],
                                        "chips": [["Box", 0, 0, 6, 6]],
                                        "focal_length": 5.6,
                                        "dither": False,
                                        "psf": ["gaussian", 72.0],
                                        "imaging": True,
                                        "grating": False}

## AXIS

instrument_registry["axis"] = {"name": "axis",
                               "arf": "axis-31jan18.arf",
                               "rmf": "axis-31jan18.rmf",
                               "bkgnd": "axis",
                               "num_pixels": 5200,
                               "fov": 15.0,
                               "aimpt_coords": [0.0, 0.0],
                               "chips": [["Box", 0, 0, 5200, 5200]],
                               "focal_length": 9.5,
                               "dither": False,
                               "psf": ["gaussian", 0.3],
                               "imaging": True,
                               "grating": False}


def add_instrument_to_registry(inst_spec):
    """
    Add an instrument specification to the registry, contained
    in either a dictionary or a JSON file.

    The *inst_spec* must have the structure as shown below. 
    The order is not important. If you use a JSON file, the
    structure is the same, but the file cannot include comments,
    and use "null" instead of "None", and "true" or "false"
    instead of "True" or "False".

    For the "chips" entry, "None" means no chips and the detector
    field of view is a single square. If you want to have multiple
    chips, they must be specified in a format described in the 
    online documentation.

    >>> {
    ...     "name": "lynx_hdxi", # The short name of the instrument
    ...     "arf": "xrs_hdxi_3x10.arf", # The file containing the ARF
    ...     "rmf": "xrs_hdxi.rmf", # The file containing the RMF
    ...     "bkgnd": "acisi", # The name of the particle background
    ...     "fov": 20.0, # The field of view in arcminutes
    ...     "focal_length": 10.0, # The focal length in meters
    ...     "num_pixels": 4096, # The number of pixels on a side in the FOV
    ...     "dither": True, # Whether or not to dither the instrument
    ...     "psf": ["gaussian", 0.5], # The type of PSF and its HPD
    ...     "chips": [["Box", 0, 0, 4096, 4096]], # The specification for the chips
    ...     "aimpt_coords": [0.0, 0.0], # The detector coordinates of the aimpoint
    ...     "imaging": True # Whether or not this is a imaging instrument
    ...     "grating": False # Whether or not this is a grating instrument
    ... }
    """
    if isinstance(inst_spec, dict):
        inst = inst_spec
    elif os.path.exists(inst_spec):
        with open(inst_spec, "r") as f:
            inst = json.load(f)
    name = inst["name"]
    if name in instrument_registry:
        raise KeyError("The instrument with name %s is already in the registry! Assign a different name!" % name)
    # Catch older JSON files which don't distinguish between imagings and non-imagings
    if "imaging" not in inst:
        mylog.warning("Instrument specifications must now include an 'imaging' item, which "
                      "determines whether or not this instrument specification supports "
                      "imaging. Default is True.")
        inst["imaging"] = True
    if "grating" not in inst:
        mylog.warning("Instrument specifications must now include an 'grating' item, which "
                      "determines whether or not this instrument specification corresponds "
                      "to a gratings instrument. Default is False.")
        inst["grating"] = False
    if inst["grating"] and inst["imaging"]:
        raise RuntimeError("Currently, gratings instrument specifications cannot have "
                           "'imaging' == True!")
    if inst['imaging']:
        default_set = {"name", "arf", "rmf", "bkgnd", "fov", "chips",
                       "aimpt_coords", "focal_length", "num_pixels",
                       "dither", "psf", "imaging", "grating"}
    else:
        default_set = {"name", "arf", "rmf", "bkgnd", "focal_length", 
                       "imaging", "grating"}
    my_keys = set(inst.keys())
    if my_keys != default_set:
        missing = default_set.difference(my_keys)
        raise RuntimeError("One or more items is missing from the instrument specification!\n"
                           "Items needed: %s" % missing)
    instrument_registry[name] = inst
    mylog.debug("The %s instrument specification has been added to the instrument registry." % name)
    return name


def get_instrument_from_registry(name):
    """
    Returns a copy of the instrument specification
    corresponding to *name*.
    """
    if name not in instrument_registry:
        raise KeyError("Instrument '%s' not in registry!" % name)
    return deepcopy(instrument_registry[name])


def show_instrument_registry():
    """
    Print the contents of the instrument registry.
    """
    for name, spec in instrument_registry.items():
        print("Instrument: %s" % name)
        for k, v in spec.items():
            print("    %s: %s" % (k, v))


def write_instrument_json(inst_name, filename):
    """
    Write an instrument specification to a JSON file.
    Useful if one would like to create a new specification
    by editing an existing one. 

    Parameters
    ----------
    inst_name : string
        The instrument specification to write.
    filename : string
        The filename to write to.
    """
    inst_dict = instrument_registry[inst_name]
    with open(filename, 'w') as f:
        json.dump(inst_dict, f, indent=4)


def make_simple_instrument(base_inst, new_inst, fov, num_pixels,
                           no_bkgnd=False, no_psf=False, no_dither=False):
    """
    Using an existing imaging instrument specification, 
    make a simple square instrument given a field of view 
    and a resolution.

    Parameters
    ----------
    base_inst : string
        The name for the instrument specification to base the 
        new one on.
    new_inst : string
        The name for the new instrument specification.
    fov : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The field of view in arcminutes.
    num_pixels : integer
        The number of pixels on a side.
    no_bkgnd : boolean, optional
        Set this new instrument to have no particle background. 
        Default: False
    no_psf : boolean, optional
        Set this new instrument to have no spatial PSF. 
        Default: False
    no_dither : boolean, optional
        Set this new instrument to have no dithering. 
        Default: False
    """
    sq_inst = get_instrument_from_registry(base_inst)
    if sq_inst["imaging"] is False:
        raise RuntimeError("make_simple_instrument only works with "
                           "imaging instruments!")
    sq_inst["name"] = new_inst
    sq_inst["chips"] = [["Box", 0, 0, num_pixels, num_pixels]]
    sq_inst["fov"] = parse_value(fov, "arcmin")
    sq_inst["num_pixels"] = num_pixels
    if no_bkgnd:
        sq_inst["bkgnd"] = None
    elif base_inst.startswith("aciss"):
        # Special-case ACIS-S to use the BI background on S3
        sq_inst["bkgnd"] = "aciss"
    if no_psf:
        sq_inst["psf"] = None
    if sq_inst["dither"]:
        sq_inst["dither"] = not no_dither
    add_instrument_to_registry(sq_inst)