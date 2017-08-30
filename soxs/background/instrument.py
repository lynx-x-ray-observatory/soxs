import os
from soxs.utils import soxs_files_path, parse_prng, \
    parse_value, mylog
from soxs.background.spectra import \
    InstrumentalBackgroundSpectrum
from soxs.background.events import make_diffuse_background
import numpy as np
from six import string_types

# ACIS-I particle background
acisi_bkgnd_file = os.path.join(soxs_files_path, "acisi_particle_bkgnd.h5")
acisi_particle_bkgnd = InstrumentalBackgroundSpectrum.from_file(acisi_bkgnd_file, 10.0)

# ACIS-S BI particle background
aciss_bkgnd_file = os.path.join(soxs_files_path, "aciss_particle_bkgnd.h5")
aciss_particle_bkgnd = InstrumentalBackgroundSpectrum.from_file(aciss_bkgnd_file, 10.0)

# Athena-like microcalorimeter background 
# (http://adsabs.harvard.edu/abs/2014A%26A...569A..54L)
mucal_bkgnd_file = os.path.join(soxs_files_path, "mucal_particle_bkgnd.h5")
mucal_particle_bkgnd = InstrumentalBackgroundSpectrum.from_file(mucal_bkgnd_file, 10.0)

# Athena microcalorimeter background 
# (http://adsabs.harvard.edu/abs/2014A%26A...569A..54L)
xifu_bkgnd_file = os.path.join(soxs_files_path, "xifu_bkgnd.h5")
athena_xifu_bkgnd = InstrumentalBackgroundSpectrum.from_file(xifu_bkgnd_file, 12.0)

# Athena imager background 
wfi_bkgnd_file = os.path.join(soxs_files_path, "wfi_bkgnd.h5")
athena_wfi_bkgnd = InstrumentalBackgroundSpectrum.from_file(wfi_bkgnd_file, 12.0)

# Hitomi SXS background
sxs_bkgnd_file = os.path.join(soxs_files_path, "hitomi_sxs_bkgnd.h5")
hitomi_sxs_bkgnd = InstrumentalBackgroundSpectrum.from_file(sxs_bkgnd_file, 5.6)

# AXIS wide-field imager background
axis_bkgnd_file = os.path.join(soxs_files_path, "axis_leo_bkgnd.h5")
axis_bkgnd = InstrumentalBackgroundSpectrum.from_file(axis_bkgnd_file, 9.5)

instrument_backgrounds = {"acisi": acisi_particle_bkgnd,
                          "aciss": aciss_particle_bkgnd,
                          "mucal": mucal_particle_bkgnd,
                          "athena_xifu": athena_xifu_bkgnd,
                          "athena_wfi": athena_wfi_bkgnd,
                          "hitomi_sxs": hitomi_sxs_bkgnd,
                          "axis": axis_bkgnd}

def add_instrumental_background(name, filename, default_focal_length):
    """
    Add a particle/instrument background to the list 
    of known backgrounds.

    Parameters
    ----------
    name : string
        The short name of the background, which will 
        be the key in the registry.
    filename : string
        The file containing the background. It must 
        have two columns: energy in keV, and background 
        intensity in units of photons/s/cm**2/arcmin**2/keV.
    default_focal_length : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The focal length of the telescope that this background
        is scaled to. Used for rescaling the background if an
        alternative focal length is provided in an instrument
        specification.
    """
    default_focal_length = parse_value(default_focal_length, "m")
    spec = InstrumentalBackgroundSpectrum.from_file(filename, default_focal_length)
    instrument_backgrounds[name] = spec

def make_instrument_background(bkgnd_name, event_params, focal_length, rmf, 
                               prng=None):
    import pyregion._region_filter as rfilter

    prng = parse_prng(prng)

    if event_params["chips"] is None:
        bkgnd_spec = [instrument_backgrounds[bkgnd_name]]
    else:
        nchips = len(event_params["chips"])
        if isinstance(bkgnd_name, string_types):
            bkgnd_spec = [instrument_backgrounds[bkgnd_name]]*nchips
        else:
            bkgnd_spec = [instrument_backgrounds[name] for name in bkgnd_name]

    bkg_events = {}

    nx = event_params["num_pixels"]

    if event_params["chips"] is None:
        bkg_events["energy"] = bkgnd_spec[0].generate_energies(event_params["exposure_time"],
                                                               event_params["fov"], prng=prng,
                                                               quiet=True).value
        n_events = bkg_events["energy"].size
        bkg_events["chip_id"] = np.zeros(n_events, dtype='int')
        bkg_events["detx"] = prng.uniform(low=-0.5*nx, high=0.5*nx, size=n_events)
        bkg_events["dety"] = prng.uniform(low=-0.5*nx, high=0.5*nx, size=n_events)
    else:
        bkg_events["energy"] = []
        bkg_events["detx"] = []
        bkg_events["dety"] = []
        bkg_events["chip_id"] = []
        for i, chip in enumerate(event_params["chips"]):
            e = bkgnd_spec[i].generate_energies(event_params["exposure_time"],
                                                event_params["fov"], prng=prng,
                                                quiet=True).value
            n_events = e.size
            detx = prng.uniform(low=-0.5*nx, high=0.5*nx, size=n_events)
            dety = prng.uniform(low=-0.5*nx, high=0.5*nx, size=n_events)
            thisc = np.ones(n_events, dtype='bool')
            rtype = chip[0]
            args = chip[1:]
            r = getattr(rfilter, rtype)(*args)
            inside = r.inside(detx, dety)
            thisc = np.logical_and(thisc, inside)
            bkg_events["energy"].append(e[thisc])
            bkg_events["detx"].append(detx[thisc])
            bkg_events["dety"].append(dety[thisc])
            bkg_events["chip_id"].append(i*np.ones(thisc.sum()))
        for key in bkg_events:
            bkg_events[key] = np.concatenate(bkg_events[key])

    if bkg_events["energy"].size == 0:
        raise RuntimeError("No instrumental background events were detected!!!")
    else:
        mylog.info("Making %d events from the instrumental background." % bkg_events["energy"].size)

    return make_diffuse_background(bkg_events, event_params, rmf, prng=prng)
