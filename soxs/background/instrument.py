import os
from soxs.utils import soxs_files_path, parse_prng, \
    parse_value, mylog
from soxs.background.spectra import \
    InstrumentalBackgroundSpectrum
from soxs.background.events import make_uniform_background

# ACIS-I particle background
acisi_bkgnd_file = os.path.join(soxs_files_path, "acisi_particle_bkgnd.h5")
acisi_particle_bkgnd = InstrumentalBackgroundSpectrum.from_file(acisi_bkgnd_file, 10.0)

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

instrument_backgrounds = {"acisi": acisi_particle_bkgnd,
                          "mucal": mucal_particle_bkgnd,
                          "athena_xifu": athena_xifu_bkgnd,
                          "athena_wfi": athena_wfi_bkgnd,
                          "hitomi_sxs": hitomi_sxs_bkgnd}

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
    prng = parse_prng(prng)

    bkgnd_spec = instrument_backgrounds[bkgnd_name]

    # Generate background events

    energy = bkgnd_spec.generate_energies(event_params["exposure_time"],
                                          event_params["fov"], 
                                          focal_length=focal_length,
                                          prng=prng, quiet=True).value

    if energy.size == 0:
        raise RuntimeError("No instrumental background events were detected!!!")
    else:
        mylog.info("Making %d events from the instrumental background." % energy.size)

    return make_uniform_background(energy, event_params, rmf, prng=prng)
