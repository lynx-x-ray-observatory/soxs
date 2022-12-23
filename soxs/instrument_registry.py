import json
import os
from copy import deepcopy

from soxs.utils import PoochHandle, mylog, parse_value

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

    def fetch_files(self, key, loc=None):
        """
        A handy method to fetch ARF, RMF, background,
        and PSF files to a location of one's choice.
        Files are only actually downloaded if they are
        not present already.

        Parameters
        ----------
        key : string
            The instrument specification to download
            the files for.
        loc : string, optional
            The path to download the files to. If not
            specified, it will download them to the
            current working directory.
        """
        inst_spec = self[key]
        if loc is None:
            loc = os.getcwd()
        dog = PoochHandle(cache_dir=loc)
        log_msg = f'Downloading %s "%s" for instrument "{key}".'
        fns = [inst_spec["arf"], inst_spec["rmf"]]
        logs = ["ARF", "RMF"]
        if inst_spec["bkgnd"] is not None:
            bkgnd = inst_spec["bkgnd"][0]
            if isinstance(bkgnd, list):
                for b in inst_spec["bkgnd"]:
                    fns.append(b[0])
            else:
                fns.append(bkgnd)
            logs.append("instrumental background model")
        if inst_spec["psf"] is not None:
            if "image" in inst_spec["psf"][0]:
                fns.append(inst_spec["psf"][1])
                logs.append("PSF model")
        for fn, log in zip(fns, logs):
            mylog.info(log_msg, log, fn)
            dog.fetch(fn)


instrument_registry = InstrumentRegistry()

# Lynx High-Definition X-ray Imager (HDXI)

instrument_registry["lynx_hdxi"] = {
    "name": "lynx_hdxi",
    "arf": "xrs_hdxi_3x10.arf",
    "rmf": "xrs_hdxi.rmf",
    "bkgnd": ["lynx_hdxi_particle_bkgnd.pha", 1.0],
    "fov": 22.0,
    "num_pixels": 4096,
    "aimpt_coords": [0.0, 0.0],
    "chips": [["Box", 0, 0, 4096, 4096]],
    "focal_length": 10.0,
    "dither": True,
    "psf": ["image", "chandra_psf.fits", 6],
    "imaging": True,
    "grating": False,
}

# Lynx Micro-calorimeter

instrument_registry["lynx_lxm"] = {
    "name": "lynx_lxm",
    "arf": "xrs_mucal_3x10_3.0eV.arf",
    "rmf": "xrs_mucal_3.0eV.rmf",
    "bkgnd": ["lynx_lxm_particle_bkgnd.pha", 1.0],
    "fov": 5.0,
    "num_pixels": 300,
    "aimpt_coords": [0.0, 0.0],
    "chips": [["Box", 0, 0, 300, 300]],
    "focal_length": 10.0,
    "dither": True,
    "psf": ["image", "chandra_psf.fits", 6],
    "imaging": True,
    "grating": False,
}

instrument_registry["lynx_lxm_enh"] = {
    "name": "lynx_lxm_enh",
    "arf": "xrs_mucal_3x10_1.5eV.arf",
    "rmf": "xrs_mucal_1.5eV.rmf",
    "bkgnd": ["lynx_lxm_enh_particle_bkgnd.pha", 1.0],
    "fov": 1.0,
    "num_pixels": 120,
    "aimpt_coords": [0.0, 0.0],
    "chips": [["Box", 0, 0, 120, 120]],
    "focal_length": 10.0,
    "dither": True,
    "psf": ["image", "chandra_psf.fits", 6],
    "imaging": True,
    "grating": False,
}

instrument_registry["lynx_lxm_ultra"] = {
    "name": "lynx_lxm_ultra",
    "arf": "xrs_mucal_3x10_0.3eV.arf",
    "rmf": "xrs_mucal_0.3eV.rmf",
    "bkgnd": ["lynx_lxm_ultra_particle_bkgnd.pha", 1.0],
    "fov": 1.0,
    "num_pixels": 60,
    "aimpt_coords": [0.0, 0.0],
    "chips": [["Box", 0, 0, 60, 60]],
    "focal_length": 10.0,
    "dither": True,
    "psf": ["image", "chandra_psf.fits", 6],
    "imaging": True,
    "grating": False,
}

# Lynx Gratings (for spectra only)

instrument_registry["lynx_xgs"] = {
    "name": "lynx_xgs",
    "arf": "xrs_cat.arf",
    "rmf": "xrs_cat.rmf",
    "bkgnd": None,
    "focal_length": 10.0,
    "imaging": False,
    "grating": True,
}

# Athena WFI

instrument_registry["athena_wfi"] = {
    "name": "athena_wfi",
    "arf": "athena_sixte_wfi_wo_filter_v20190122.arf",
    "rmf": "athena_wfi_sixte_v20150504.rmf",
    "bkgnd": ["sixte_wfi_particle_bkg_20190829.pha", 79552.92570677],
    "fov": 40.147153,
    "num_pixels": 1078,
    "aimpt_coords": [53.69, -53.69],
    "chips": [
        ["Box", -283, -283, 512, 512],
        ["Box", 283, -283, 512, 512],
        ["Box", -283, 283, 512, 512],
        ["Box", 283, 283, 512, 512],
    ],
    "focal_length": 12.0,
    "dither": True,
    "psf": ["multi_image", "athena_psf_15row.fits"],
    "imaging": True,
    "grating": False,
}

# Athena XIFU

instrument_registry["athena_xifu"] = {
    "name": "athena_xifu",
    "arf": "sixte_xifu_cc_baselineconf_20180821.arf",
    "rmf": "XIFU_CC_BASELINECONF_2018_10_10.rmf",
    "bkgnd": ["xifu_nxb_20181209.pha", 79552.92570677],
    "fov": 5.991992621478149,
    "num_pixels": 84,
    "aimpt_coords": [0.0, 0.0],
    "chips": [["Polygon", [-33, 0, 33, 33, 0, -33], [20, 38, 20, -20, -38, -20]]],
    "focal_length": 12.0,
    "dither": True,
    "psf": ["multi_image", "athena_psf_15row.fits"],
    "imaging": True,
    "grating": False,
}

# Chandra ACIS-I, Cycle 0 and 20

for cycle in [0, 22]:
    name = f"chandra_acisi_cy{cycle}"
    instrument_registry[name] = {
        "name": name,
        "arf": f"acisi_aimpt_cy{cycle}.arf",
        "rmf": f"acisi_aimpt_cy{cycle}.rmf",
        "bkgnd": [f"chandra_acisi_cy{cycle}_particle_bkgnd.pha", 1.0],
        "fov": 20.008,
        "num_pixels": 2440,
        "aimpt_coords": [86.0, 57.0],
        "chips": [
            ["Box", -523, -523, 1024, 1024],
            ["Box", 523, -523, 1024, 1024],
            ["Box", -523, 523, 1024, 1024],
            ["Box", 523, 523, 1024, 1024],
        ],
        "psf": ["multi_image", "chandra_psf.fits"],
        "focal_length": 10.0,
        "dither": True,
        "imaging": True,
        "grating": False,
    }

# Chandra ACIS-S, Cycle 0 and 22

for cycle in [0, 22]:
    name = f"chandra_aciss_cy{cycle}"
    instrument_registry[name] = {
        "name": name,
        "arf": f"aciss_aimpt_cy{cycle}.arf",
        "rmf": f"aciss_aimpt_cy{cycle}.rmf",
        "bkgnd": [
            [f"chandra_aciss_cy{cycle}_fi_particle_bkgnd.pha", 1.0],
            [f"chandra_aciss_cy{cycle}_bi_particle_bkgnd.pha", 1.0],
            [f"chandra_aciss_cy{cycle}_fi_particle_bkgnd.pha", 1.0],
            [f"chandra_aciss_cy{cycle}_bi_particle_bkgnd.pha", 1.0],
            [f"chandra_aciss_cy{cycle}_fi_particle_bkgnd.pha", 1.0],
            [f"chandra_aciss_cy{cycle}_fi_particle_bkgnd.pha", 1.0],
        ],
        "fov": 50.02,
        "num_pixels": 6100,
        "aimpt_coords": [206.0, 0.0],
        "chips": [
            ["Box", -2605, 0, 1024, 1024],
            ["Box", -1563, 0, 1024, 1024],
            ["Box", -521, 0, 1024, 1024],
            ["Box", 521, 0, 1024, 1024],
            ["Box", 1563, 0, 1024, 1024],
            ["Box", 2605, 0, 1024, 1024],
        ],
        "psf": ["multi_image", "chandra_psf.fits"],
        "focal_length": 10.0,
        "dither": True,
        "imaging": True,
        "grating": False,
    }


# Chandra ACIS-S, Cycle 0 and 19 HETG (for spectra only)

orders = {"p1": 1, "m1": -1}

for energy in ["meg", "heg"]:
    for order in ["p1", "m1"]:
        for cycle in [0, 22]:
            name = f"chandra_aciss_{energy}_{order}_cy{cycle}"
            resp_name = f"chandra_aciss_{energy}{orders[order]}_cy{cycle}"
            instrument_registry[name] = {
                "name": name,
                "arf": f"{resp_name}.garf",
                "rmf": f"{resp_name}.grmf",
                "bkgnd": None,
                "focal_length": 10.0,
                "imaging": False,
                "grating": True,
            }

# XRISM Resolve

instrument_registry["xrism_resolve"] = {
    "name": "xrism_resolve",
    "arf": "resolve_pnt_heasim_noGV_20190701.arf",
    "rmf": "resolve_h5ev_2019a.rmf",
    "bkgnd": ["resolve_h5ev_2019a_rslnxb.pha", 9.130329009932256],
    "num_pixels": 6,
    "fov": 3.06450576,
    "aimpt_coords": [0.0, 0.0],
    "chips": [["Box", 0, 0, 6, 6]],
    "focal_length": 5.6,
    "dither": False,
    "psf": ["multi_image", "sxs_psfimage_20140618.fits"],
    "imaging": True,
    "grating": False,
}

instrument_registry["xrism_resolve_withGV"] = deepcopy(
    instrument_registry["xrism_resolve"]
)
instrument_registry["xrism_resolve_withGV"][
    "arf"
] = "resolve_pnt_heasim_withGV_20190701.arf"

# XRISM Xtend

instrument_registry["xrism_xtend"] = {
    "name": "xrism_xtend",
    "arf": "sxt-i_140505_ts02um_int01.8r_intall_140618psf.arf",
    "rmf": "ah_sxi_20120702.rmf",
    "bkgnd": ["ah_sxi_pch_nxb_full_20110530.pi", 1422.6292229683816],
    "num_pixels": 1296,
    "fov": 38.18845555660526,
    "aimpt_coords": [-244.0, -244.0],
    "chips": [
        ["Box", -327, 327, 640, 640],
        ["Box", -327, -327, 640, 640],
        ["Box", 327, 327, 640, 640],
        ["Box", 327, -327, 640, 640],
    ],
    "focal_length": 5.6,
    "dither": False,
    "psf": ["eef", "eef_from_sxi_psfimage_20140618.fits", 1],
    "imaging": True,
    "grating": False,
}

# AXIS

instrument_registry["axis"] = {
    "name": "axis",
    "arf": "axis_onaxis_20221116.arf",
    "rmf": "axis_ccd_20221101.rmf",
    "bkgnd": ["axis_nxb_FOV_10Msec_20221215.pha", 697.06],
    "num_pixels": 2952,
    "fov": 27.06194257961904,
    "aimpt_coords": [-109, 109],
    "chips": [
        ["Box", -756, -756, 1440, 1440],
        ["Box", -756, 756, 1440, 1440],
        ["Box", 756, -756, 1440, 1440],
        ["Box", 756, 756, 1440, 1440],
    ],
    "focal_length": 9.0,
    "dither": False,
    "psf": ["multi_eef", "AXIS_EEF_2022-02-16.fits", 2],
    "imaging": True,
    "grating": False,
}

# STAR-X

instrument_registry["star-x"] = {
    "name": "star-x",
    "arf": "starx_2020-11-26_fov_avg.arf",
    "rmf": "starx.rmf",
    "bkgnd": None,
    "num_pixels": 3600,
    "fov": 60.0,
    "aimpt_coords": [0.0, 0.0],
    "chips": [["Box", 0, 0, 3600, 3600]],
    "focal_length": 4.5,
    "dither": True,
    "psf": ["gaussian", 3.0],
    "imaging": True,
    "grating": False,
}

# LEM

instrument_registry["lem_2eV"] = {
    "name": "lem_2eV",
    "arf": "lem_300522.arf",
    "rmf": "lem_2ev_110422.rmf",
    "bkgnd": ["lem_2eV_171222_fov_bkg.pi", 900.0],
    "num_pixels": 128,
    "fov": 32.0,
    "aimpt_coords": [0.0, 0.0],
    "chips": [["Box", 0, 0, 128, 128]],
    "focal_length": 4.0,
    "dither": True,
    "psf": ["gaussian", 10.0],
    "imaging": True,
    "grating": False,
}

instrument_registry["lem_0.9eV"] = {
    "name": "lem_0.9eV",
    "arf": "lem_300522.arf",
    "rmf": "lem_09ev_110422.rmf",
    "bkgnd": ["lem_09eV_171222_fov_bkg.pi", 900.0],
    "num_pixels": 128,
    "fov": 32.0,
    "aimpt_coords": [0.0, 0.0],
    "chips": [["Box", 0, 0, 128, 128]],
    "focal_length": 4.0,
    "dither": True,
    "psf": ["gaussian", 10.0],
    "imaging": True,
    "grating": False,
}

instrument_registry["lem_2eV_0422"] = deepcopy(instrument_registry["lem_2eV"])
instrument_registry["lem_2eV_0422"]["arf"] = "lem_110422.arf"
instrument_registry["lem_0.9eV_0422"] = deepcopy(instrument_registry["lem_0.9eV"])
instrument_registry["lem_0.9eV_0422"]["arf"] = "lem_110422.arf"
instrument_registry["lem_2eV_0322"] = deepcopy(instrument_registry["lem_2eV"])
instrument_registry["lem_2eV_0322"]["arf"] = "lem_030322a.arf"
instrument_registry["lem_2eV_0322"]["rmf"] = "lem_2ev_030322.rmf"
instrument_registry["lem_0.9eV_0322"] = deepcopy(instrument_registry["lem_0.9eV"])
instrument_registry["lem_0.9eV_0322"]["arf"] = "lem_030322a.arf"
instrument_registry["lem_0.9eV_0322"]["rmf"] = "lem_09ev_030322.rmf"


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
    ...     "bkgnd": ["lynx_hdxi_particle_bkgnd.pha", 1.0], # The name of the particle background file and the area of extraction
    ...     "fov": 20.0, # The field of view in arcminutes
    ...     "focal_length": 10.0, # The focal length in meters
    ...     "num_pixels": 4096, # The number of pixels on a side in the FOV
    ...     "dither": True, # Whether to dither the instrument
    ...     "psf": ["image", "chandra_psf.fits", 6], # The type of PSF and associated parameters
    ...     "chips": [["Box", 0, 0, 4096, 4096]], # The specification for the chips
    ...     "aimpt_coords": [0.0, 0.0], # The detector coordinates of the aimpoint
    ...     "imaging": True # Whether this is an imaging instrument
    ...     "grating": False # Whether this is a grating instrument
    ... }
    """
    if isinstance(inst_spec, dict):
        inst = inst_spec
    elif os.path.exists(inst_spec):
        with open(inst_spec, "r") as f:
            inst = json.load(f)
    name = inst["name"]
    if name in instrument_registry:
        raise KeyError(
            f"The instrument with name {name} is already in the "
            f"registry! Assign a different name!"
        )
    # Catch older JSON files which don't distinguish between imagings
    # and non-imagings
    if "imaging" not in inst:
        mylog.warning(
            "Instrument specifications must now include an 'imaging' "
            "item, which determines whether or not this instrument "
            "specification supports imaging. Default is True."
        )
        inst["imaging"] = True
    if "grating" not in inst:
        mylog.warning(
            "Instrument specifications must now include an 'grating' "
            "item, which determines whether or not this instrument "
            "specification corresponds to a gratings instrument. "
            "Default is False."
        )
        inst["grating"] = False
    if inst["grating"] and inst["imaging"]:
        raise RuntimeError(
            "Currently, gratings instrument specifications cannot "
            "have 'imaging' == True!"
        )
    if inst["imaging"]:
        default_set = {
            "name",
            "arf",
            "rmf",
            "bkgnd",
            "fov",
            "chips",
            "aimpt_coords",
            "focal_length",
            "num_pixels",
            "dither",
            "psf",
            "imaging",
            "grating",
        }
    else:
        default_set = {
            "name",
            "arf",
            "rmf",
            "bkgnd",
            "focal_length",
            "imaging",
            "grating",
        }
    if inst["imaging"]:
        if inst["psf"] is not None and not isinstance(inst["psf"], list):
            raise RuntimeError(
                "The 'psf' option in the instrument needs to be "
                "a two-element list specifying the PSF type and "
                "additional arguments. See the SOXS documentation "
                "for details."
            )
    if inst["bkgnd"] is not None and not isinstance(inst["bkgnd"], list):
        raise RuntimeError(
            "The 'bkgnd' option in the instrument needs to be "
            "either a two-element list specifying the background "
            "file and its area in square arcminutes, or a list of "
            "such two-element lists if the background is different "
            "on different chips. See the SOXS documentation for "
            "details."
        )
    if "dep_name" in inst:
        mylog.warning(
            "The 'dep_name' option is no longer supported. Dropping it "
            "from the instrument specification."
        )
        inst.pop("dep_name")
    my_keys = set(inst.keys())
    if my_keys != default_set:
        missing = default_set.difference(my_keys)
        raise RuntimeError(
            f"One or more items is missing from the instrument "
            f"specification!\nItems needed: {missing}"
        )
    instrument_registry[name] = inst
    mylog.debug(
        "The %s instrument specification has been added " "to the instrument registry.",
        name,
    )
    return name


def get_instrument_from_registry(name):
    """
    Returns a copy of the instrument specification
    corresponding to *name*.
    """
    if name not in instrument_registry:
        raise KeyError(f"Instrument '{name}' not in registry!")
    return deepcopy(instrument_registry[name])


def show_instrument_registry():
    """
    Print the contents of the instrument registry.
    """
    for name, spec in instrument_registry.items():
        print(f"Instrument: {name}")
        for k, v in spec.items():
            print(f"    {k}: {v}")


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
    with open(filename, "w") as f:
        json.dump(inst_dict, f, indent=4)


def make_simple_instrument(
    base_inst, new_inst, fov, num_pixels, no_bkgnd=False, no_psf=False, no_dither=False
):
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
        raise RuntimeError(
            "make_simple_instrument only works with " "imaging instruments!"
        )
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
