from ._version import get_versions

from soxs.utils import soxs_cfg

from soxs.simput import \
    write_photon_list, \
    read_simput_catalog, \
    PhotonList, SimputCatalog

from soxs.spectra import \
    Spectrum, \
    ApecGenerator, \
    ConvolvedSpectrum

from soxs.instrument import \
    instrument_simulator, \
    AuxiliaryResponseFile, \
    RedistributionMatrixFile, \
    make_background_file, \
    FlatResponse, simulate_spectrum

from soxs.instrument_registry import \
    add_instrument_to_registry, \
    show_instrument_registry, \
    write_instrument_json, \
    get_instrument_from_registry, \
    instrument_registry, \
    make_simple_instrument

from soxs.spatial import \
    PointSourceModel, \
    RadialFunctionModel, \
    RadialArrayModel, \
    RadialFileModel, \
    AnnulusModel, \
    BetaModel, \
    FillFOVModel, \
    RectangleModel, \
    SpatialModel

from soxs.events import \
    write_spectrum, \
    write_image, \
    write_radial_profile, \
    plot_spectrum, \
    make_exposure_map

from soxs.background import \
    add_instrumental_background, \
    BackgroundSpectrum, \
    ConvolvedBackgroundSpectrum, \
    InstrumentalBackgroundSpectrum

from soxs.cosmology import \
    make_cosmological_sources_file

from soxs.background import \
    make_point_sources_file, \
    make_point_source_list

__version__ = get_versions()['version']
del get_versions
