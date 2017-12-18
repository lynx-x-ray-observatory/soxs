from soxs.background.spectra import \
    GratingsBackgroundSpectrum
from soxs.utils import soxs_files_path
import os

# ACIS-S HETG background
aciss_hetg_bkgnd_file = os.path.join(soxs_files_path, "aciss_hetg_bkgnd.h5")
aciss_hetg_bkgnd = GratingsBackgroundSpectrum.from_file(aciss_hetg_bkgnd_file)

# Lynx gratings background
lynx_gratings_bkgnd_file = os.path.join(soxs_files_path, "lynx_gratings_bkgnd.h5")
lynx_gratings_bkgnd = GratingsBackgroundSpectrum.from_file(lynx_gratings_bkgnd_file)

# Arcus background
arcus_bkgnd_file = os.path.join(soxs_files_path, "arcus_bkgnd.h5")
arcus_bkgnd = GratingsBackgroundSpectrum.from_file(arcus_bkgnd_file)

gratings_backgrounds = {"aciss_hetg": aciss_hetg_bkgnd,
                        "lynx_gratings": lynx_gratings_bkgnd,
                        "arcus": arcus_bkgnd}

def add_instrumental_background(name, filename):
    """
    Add a gratings background to the list of known backgrounds.

    Parameters
    ----------
    name : string
        The short name of the background, which will 
        be the key in the registry.
    filename : string
        The file containing the background. It must 
        have two columns: energy in keV, and background 
        count rate in units of photons/s/keV.
    """
    spec = GratingsBackgroundSpectrum.from_file(filename)
    gratings_backgrounds[name] = spec
