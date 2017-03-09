import os
from soxs.utils import soxs_files_path, parse_prng
from soxs.background.spectra import BackgroundSpectrum
from soxs.background.events import make_uniform_background

class InstrumentalBackgroundSpectrum(BackgroundSpectrum):
    def __init__(self, filename, default_focal_length):
        super(InstrumentalBackgroundSpectrum, self).__init__(filename)
        self.default_focal_length = default_focal_length

# ACIS-I particle background
acisi_bkgnd_file = os.path.join(soxs_files_path, "acisi_particle_bkgnd.dat")
acisi_particle_bkgnd = InstrumentalBackgroundSpectrum(acisi_bkgnd_file, 10.0)

# Athena-like microcalorimeter background 
# (http://adsabs.harvard.edu/abs/2014A%26A...569A..54L)
mucal_bkgnd_file = os.path.join(soxs_files_path, "mucal_particle_bkgnd.dat")
mucal_particle_bkgnd = InstrumentalBackgroundSpectrum(mucal_bkgnd_file, 10.0)

# Athena microcalorimeter background 
# (http://adsabs.harvard.edu/abs/2014A%26A...569A..54L)
xifu_bkgnd_file = os.path.join(soxs_files_path, "xifu_bkgnd.dat")
athena_xifu_bkgnd = InstrumentalBackgroundSpectrum(xifu_bkgnd_file, 12.0)

# Athena imager background 
wfi_bkgnd_file = os.path.join(soxs_files_path, "wfi_bkgnd.dat")
athena_wfi_bkgnd = InstrumentalBackgroundSpectrum(wfi_bkgnd_file, 12.0)

instrument_backgrounds = {"acisi": acisi_particle_bkgnd,
                          "mucal": mucal_particle_bkgnd,
                          "athena_xifu": athena_xifu_bkgnd,
                          "athena_wfi": athena_wfi_bkgnd}

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
    default_focal_length : float
        The focal length of the telescope that this background
        is scaled to. Used for rescaling the background if an
        alternative focal length is provided in an instrument
        specification.
    """
    spec = InstrumentalBackgroundSpectrum(filename, default_focal_length)
    instrument_backgrounds[name] = spec

def make_instrument_background(bkgnd_name, event_params, focal_length, rmf, 
                               prng=None):
    prng = parse_prng(prng)
    fov = event_params["fov"]

    bkgnd_spec = instrument_backgrounds[bkgnd_name]

    # Generate background events

    area = (focal_length / bkgnd_spec.default_focal_length) ** 2
    energy = bkgnd_spec.generate_energies(event_params["exposure_time"], area, 
                                          fov, prng=prng).value

    if energy.size == 0:
        raise RuntimeError("No instrumental background events were detected!!!")

    return make_uniform_background(energy, event_params, rmf, prng=prng)
