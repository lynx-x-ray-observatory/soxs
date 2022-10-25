from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("soxs")
except PackageNotFoundError:
    # package is not installed
    pass


from soxs.background import \
    BackgroundSpectrum, \
    ConvolvedBackgroundSpectrum, \
    make_point_sources_file, \
    make_point_source_list

from soxs.cosmology import \
    make_cosmological_sources_file

from soxs.events import \
    write_spectrum, \
    write_image, \
    write_radial_profile, \
    plot_spectrum, \
    make_exposure_map, \
    plot_image, \
    filter_events

from soxs.instrument import \
    instrument_simulator, \
    make_background_file, \
    simulate_spectrum, \
    simple_event_list

from soxs.instrument_registry import \
    add_instrument_to_registry, \
    show_instrument_registry, \
    write_instrument_json, \
    get_instrument_from_registry, \
    instrument_registry, \
    make_simple_instrument

from soxs.mosaic import \
    make_mosaic_events, \
    make_mosaic_image

from soxs.response import \
    AuxiliaryResponseFile, \
    RedistributionMatrixFile, \
    FlatResponse

from soxs.simput import \
    read_simput_catalog, \
    SimputPhotonList, \
    SimputCatalog, \
    SimputSpectrum, \
    write_photon_list

from soxs.spatial import \
    PointSourceModel, \
    RadialFunctionModel, \
    RadialArrayModel, \
    RadialFileModel, \
    AnnulusModel, \
    BetaModel, \
    DoubleBetaModel, \
    FillFOVModel, \
    RectangleModel, \
    SpatialModel

from soxs.spectra import \
    Spectrum, \
    CountRateSpectrum, \
    ConvolvedSpectrum

from soxs.thermal_spectra import \
    ApecGenerator, \
    SpexGenerator, \
    IGMGenerator, \
    CloudyCIEGenerator, \
    MekalGenerator

from soxs.utils import \
    set_soxs_config, \
    set_mission_config
