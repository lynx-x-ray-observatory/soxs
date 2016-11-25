from soxs.instrument_registry import show_instrument_registry
from soxs.background import show_background_registry

def test_show():
    # Just make sure these don't throw any errors
    show_background_registry()
    show_instrument_registry()