from __future__ import print_function
import json
from soxs.utils import mylog
import os

# The Instrument Registry

instrument_registry = {}

# Lynx

# High-Definition X-ray Imager (HDXI)

instrument_registry["hdxi"] = {"name": "hdxi_3x10",
                               "arf": "xrs_hdxi_3x10.arf",
                               "rmf": "xrs_hdxi.rmf",
                               "bkgnd": "acisi",
                               "fov": 20.0,
                               "num_pixels": 4096,
                               "aimpt_coords": [0.0, 0.0],
                               "chips": None,
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
                                "aimpt_coords": [0.0, 0.0],
                                "chips": None,
                                "focal_length": 10.0,
                                "dither": True,
                                "psf": ["gaussian", 0.5]}

instrument_registry["mucal_3x10"] = instrument_registry["mucal"]

# Account for different ARFs
for det in ["hdxi", "mucal"]:
    for mirror in ["3x15", "3x20", "6x20"]:
        instrument_registry["%s_%s" % (det, mirror)] = instrument_registry[det].copy()
        instrument_registry["%s_%s" % (det, mirror)]["name"] = "%s_%s" % (det, mirror)
        instrument_registry["%s_%s" % (det, mirror)]["arf"] = "xrs_%s_%s.arf" % (det, mirror)
        instrument_registry["%s_%s" % (det, mirror)]["focal_length"] = float(mirror.split("x")[-1])

# Athena

# XIFU

instrument_registry["athena_xifu"] = {"name": "athena_xifu",
                                      "arf": "athena_xifu_1469_onaxis_pitch249um_v20160401.arf",
                                      "rmf": "athena_xifu_rmf_v20160401.rmf",
                                      "bkgnd": "athena_xifu",
                                      "fov": 5.991992621478149,
                                      "num_pixels": 84,
                                      "aimpt_coords": [0.0, 0.0],
                                      "chips": [["Polygon", 
                                                 [-33, 0, 33, 33, 0, -33],
                                                 [20, 38, 20, -20, -38, -20]]],
                                      "focal_length": 12.0,
                                      "dither": False,
                                      "psf": ["gaussian", 5.0]}

# WFI

instrument_registry["athena_wfi"] = {"name": "athena_wfi",
                                     "arf": "athena_wfi_1469_onaxis_w_filter_v20150326.arf",
                                     "rmf": "athena_wfi_rmf_v20150326.rmf",
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
                                     "psf": ["gaussian", 5.0]}

# Old specs

instrument_registry["athena_wfi_old"] = {"name": "athena_wfi_old",
                                         "arf": "athena_wfi_1469_onaxis_w_filter_v20150326.arf",
                                         "rmf": "athena_wfi_rmf_v20150326.rmf",
                                         "bkgnd": "athena_wfi",
                                         "fov": 40.0,
                                         "num_pixels": 1092,
                                         "aimpt_coords": [0.0, 0.0],
                                         "chips": None,
                                         "focal_length": 12.0,
                                         "dither": False,
                                         "psf": ["gaussian", 5.0]}

instrument_registry["athena_xifu_old"] = {"name": "athena_xifu_old",
                                          "arf": "athena_xifu_1469_onaxis_pitch249um_v20160401.arf",
                                          "rmf": "athena_xifu_rmf_v20160401.rmf",
                                          "bkgnd": "athena_xifu",
                                          "fov": 5.0,
                                          "num_pixels": 70,
                                          "aimpt_coords": [0.0, 0.0],
                                          "chips": None,
                                          "focal_length": 12.0,
                                          "dither": False,
                                          "psf": ["gaussian", 5.0]}

# Chandra

# ACIS-I, Cycle 0 

instrument_registry["acisi_cy0"] = {"name": "acisi_cy0", 
                                    "arf": "acisi_aimpt_cy0.arf",
                                    "rmf": "acisi_aimpt_cy0.rmf",
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
                                    "dither": True}

# ACIS-I, Cycle 18

instrument_registry["acisi_cy18"] = {"name": "acisi_cy18",
                                     "arf": "acisi_aimpt_cy18.arf",
                                     "rmf": "acisi_aimpt_cy18.rmf",
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
                                     "dither": True}

# Old specs

instrument_registry["acisi_cy0_old"] = {"name": "acisi_cy0_old",
                                        "arf": "acisi_aimpt_cy0.arf",
                                        "rmf": "acisi_aimpt_cy0.rmf",
                                        "bkgnd": "acisi",
                                        "fov": 16.892,
                                        "num_pixels": 2060,
                                        "aimpt_coords": [0.0, 0.0],
                                        "chips": None,
                                        "psf": ["gaussian", 0.5],
                                        "focal_length": 10.0,
                                        "dither": True}

instrument_registry["acisi_cy18_old"] = {"name": "acisi_cy18",
                                         "arf": "acisi_aimpt_cy18.arf",
                                         "rmf": "acisi_aimpt_cy18.rmf",
                                         "bkgnd": "acisi",
                                         "fov": 16.892,
                                         "num_pixels": 2060,
                                         "aimpt_coords": [0.0, 0.0],
                                         "chips": None,
                                         "psf": ["gaussian", 0.5],
                                         "focal_length": 10.0,
                                         "dither": True}

# Hitomi

# SXS 

instrument_registry["hitomi_sxs"] = {"name": "hitomi_sxs",
                                     "arf": "hitomi_sxs_ptsrc.arf",
                                     "rmf": "hitomi_sxs.rmf",
                                     "bkgnd": "hitomi_sxs",
                                     "num_pixels": 6,
                                     "fov": 3.06450576,
                                     "aimpt_coords": [0.0, 0.0],
                                     "chips": None,
                                     "focal_length": 5.6,
                                     "dither": False,
                                     "psf": ["gaussian", 72.0]}


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
    ...     "name": "hdxi_3x10", # The short name of the instrument
    ...     "arf": "xrs_hdxi_3x10.arf", # The file containing the ARF
    ...     "rmf": "xrs_hdxi.rmf", # The file containing the RMF
    ...     "bkgnd": "acisi", # The name of the particle background
    ...     "fov": 20.0, # The field of view in arcminutes
    ...     "focal_length": 10.0, # The focal length in meters
    ...     "num_pixels": 4096, # The number of pixels on a side in the FOV
    ...     "dither": True, # Whether or not to dither the instrument
    ...     "psf": ["gaussian", 0.5] # The type of PSF and its HPD
    ...     "chips": None # The specification for the chips
    ...     "aimpt_coords": [0.0, 0.0] # The detector coordinates of the aimpoint
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
    # Catch older JSON files without chip definitions
    if "chips" not in inst:
        mylog.warning("Instrument specifications must now include a 'chips' item, which details "
                      "the layout of the chips if there are more that one. Assuming None for "
                      "one chip that covers the entire field of view.")
        inst["chips"] = None
    # Catch older JSON files without aimpoint coordinates
    if "aimpt_coords" not in inst:
        mylog.warning("Instrument specifications must now include a 'aimpt_coords' item, which "
                      "details the position in detector coordinates of the nominal aimpoint. "
                      "Assuming [0.0, 0.0].")
        inst["aimpt_coords"] = [0.0, 0.0]
    default_set = {"name", "arf", "rmf", "bkgnd", "fov", "chips",
                   "aimpt_coords", "focal_length", "num_pixels", "dither", "psf"}
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
