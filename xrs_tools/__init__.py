from setuptools_scm import get_version

__version__ = get_version()

from xrs_tools.simput import \
    write_simput_phlist, \
    read_simput_phlist

from xrs_tools.spectra import \
    Spectrum, \
    ApecGenerator

from xrs_tools.instrument import \
    add_instrument_to_registry, \
    show_instrument_registry

from xrs_tools.events import \
    make_event_file, \
    add_background_events