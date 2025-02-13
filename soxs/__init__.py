from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("soxs")
except PackageNotFoundError:
    # package is not installed
    pass


from soxs.background import (
    BackgroundSpectrum,
    ConvolvedBackgroundSpectrum,
    make_point_source_list,
    make_point_sources_file,
)
from soxs.cosmology import make_cosmological_sources_file
from soxs.events import (
    fill_regions,
    filter_events,
    make_exposure_map,
    make_image,
    merge_event_files,
    plot_image,
    plot_spectrum,
    write_image,
    write_radial_profile,
    write_spectrum,
)
from soxs.instrument import (
    instrument_simulator,
    make_background_file,
    simple_event_list,
    simulate_spectrum,
)
from soxs.instrument_registry import (
    add_instrument_to_registry,
    get_instrument_from_registry,
    instrument_registry,
    make_simple_instrument,
    show_instrument_registry,
    write_instrument_json,
)
from soxs.mosaic import make_mosaic_events, make_mosaic_image
from soxs.response import AuxiliaryResponseFile, FlatResponse, RedistributionMatrixFile
from soxs.simput import (
    SimputCatalog,
    SimputPhotonList,
    SimputSpectrum,
    make_bkgnd_simput,
    write_photon_list,
)
from soxs.spatial import (
    AnnulusModel,
    BetaModel,
    DoubleBetaModel,
    FillFOVModel,
    PointSourceModel,
    RadialArrayModel,
    RadialFileModel,
    RadialFunctionModel,
    RectangleModel,
    SpatialModel,
)
from soxs.spectra import (
    ACX2Generator,
    ApecGenerator,
    CloudyCIEGenerator,
    CloudyPionGenerator,
    ConvolvedSpectrum,
    CountRateSpectrum,
    IGMGenerator,
    MekalGenerator,
    OneACX2Generator,
    Spectrum,
    SpexGenerator,
    download_spectrum_tables,
)
from soxs.utils import set_mission_config, set_soxs_config
