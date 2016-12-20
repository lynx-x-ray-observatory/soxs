from soxs.instrument_registry import show_instrument_registry
from soxs.background.instrument import show_instrument_backgrounds

def test_show():
    # Just make sure these don't throw any errors
    show_instrument_backgrounds()
    show_instrument_registry()