from __future__ import print_function
import json
from soxs.utils import mylog
import os

# The Instrument Registry

instrument_registry = {}

# Lynx

# High-Definition X-ray Imager
instrument_registry["hdxi"] = {"name": "hdxi_3x10",
                               "arf": "xrs_hdxi_3x10.arf",
                               "rmf": "xrs_hdxi.rmf",
                               "bkgnd": "acisi",
                               "fov": 20.0,
                               "num_pixels": 4096,
                               "focal_length": 10.0,
                               "dither": True,
                               "psf": ["gaussian", 0.5]}
instrument_registry["hdxi_3x10"] = instrument_registry["hdxi"]

# Micro-calorimeter
instrument_registry["mucal"] = {"name": "mucal_3x10",
                                "arf": "xrs_mucal_3x10.arf",
                                "rmf": "xrs_mucal.rmf",
                                "bkgnd": "mucal",
                                "fov": 5.0,
                                "num_pixels": 300,
                                "focal_length": 10.0,
                                "dither": True,
                                "psf": ["gaussian", 0.5]}
instrument_registry["mucal_3x10"] = instrument_registry["mucal"]
# The next line is for backwards-compatibility
instrument_registry["xcal"] = instrument_registry["mucal"]

# Account for different ARFs
for det in ["hdxi", "mucal"]:
    for mirror in ["3x15", "3x20", "6x20"]:
        instrument_registry["%s_%s" % (det, mirror)] = instrument_registry[det].copy()
        instrument_registry["%s_%s" % (det, mirror)]["name"] = "%s_%s" % (det, mirror)
        instrument_registry["%s_%s" % (det, mirror)]["arf"] = "xrs_%s_%s.arf" % (det, mirror)
        instrument_registry["%s_%s" % (det, mirror)]["focal_length"] = float(mirror.split("x")[-1])

# Athena

instrument_registry["athena_wfi"] = {"name": "athena_wfi",
                                     "arf": "athena_wfi_1469_onaxis_w_filter_v20150326.arf",
                                     "rmf": "athena_wfi_rmf_v20150326.rmf",
                                     "bkgnd": "athena_wfi",
                                     "fov": 60.0,
                                     "num_pixels": 1024,
                                     "focal_length": 12.0,
                                     "dither": False,
                                     "psf": ["gaussian", 5.0]}

instrument_registry["athena_xifu"] = {"name": "athena_xifu",
                                      "arf": "athena_xifu_1469_onaxis_pitch249um_v20160401.arf",
                                      "rmf": "athena_xifu_rmf_v20160401.rmf",
                                      "bkgnd": "athena_xifu",
                                      "fov": 10.0,
                                      "num_pixels": 66,
                                      "focal_length": 12.0,
                                      "dither": False,
                                      "psf": ["gaussian", 5.0]}

def add_instrument_to_registry(inst_spec):
    """
    Add an instrument specification to the registry, contained
    in either a dictionary or a JSON file.

    The *inst_spec* must have the structure as shown below. 
    The order is not important. If you use a JSON file, the
    structure is the same, but the file cannot include comments.

    >>> {
    ...     "name": "hdxi_3x10", # The short name of the instrument
    ...     "arf": "xrs_hdxi_3x10.arf", # The file containing the ARF
    ...     "rmf": "xrs_hdxi.rmf" # The file containing the RMF
    ...     "bkgnd": "acisi" # The name of the particle background
    ...     "fov": 20.0, # The field of view in arcminutes
    ...     "focal_length": 10.0, # The focal length in meters
    ...     "num_pixels": 4096, # The number of pixels on a side in the FOV
    ...     "dither": True, # Whether or not to dither the instrument
    ...     "psf": ["gaussian", 0.5] # The type of PSF and its HPD
    ... }
    """
    if isinstance(inst_spec, dict):
        inst = inst_spec
    elif os.path.exists(inst_spec):
        f = open(inst_spec, "r")
        inst = json.load(f)
        f.close()
    name = inst["name"]
    if name in instrument_registry:
        raise KeyError("The instrument with name %s is already in the registry! Assign a different name!" % name)
    # Catch old JSON files with plate scale
    if "plate_scale" in inst:
        mylog.warning("Instrument specifications with the 'plate_scale' item are deprecated, and will "
                      "not work in a future release. Please specify the field of view in arcminutes "
                      "with 'fov' instead.")
        inst["fov"] = inst["num_pixels"]*inst["plate_scale"]/60.0
        inst.pop("plate_scale")
    default_set = {"name", "arf", "rmf", "bkgnd", "fov", 
                   "focal_length", "num_pixels", "dither", "psf"}
    if set(inst.keys()) != default_set:
        raise RuntimeError("One or more items is missing from the instrument specification!\n"
                           "Items present: %s\nItems needed: %s" % (set(inst.keys()), default_set))
    instrument_registry[name] = inst
    mylog.debug("The %s instrument specification has been added to the instrument registry." % name)
    return name

def get_instrument_from_registry(name):
    """
    Returns a copy of the instrument specification
    corresponding to *name*.
    """
    return instrument_registry[name].copy()

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
    fp = open(filename, 'w')
    json.dump(inst_dict, fp, indent=4)
    fp.close()
