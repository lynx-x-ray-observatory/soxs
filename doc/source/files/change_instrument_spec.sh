#!/bin/bash

# First, we make the SIMPUT file using the Python script which calls pyXSIM
python make_sloshing.py

# Next, we make three event files, using a different instrument specification for each

# Normal HDXI with 0.5 arcsec PSF
instrument_simulator sloshing_simput.fits evt.fits 50000.0 hdxi 30.0,45.0 --clobber

# HDXI with 2 arcsec PSF
instrument_simulator sloshing_simput.fits evt_2.fits 50000.0 hdxi_2.json 30.0,45.0 --clobber

# HDXI with 5 arcsec PSF
instrument_simulator sloshing_simput.fits evt_5.fits 50000.0 hdxi_5.json 30.0,45.0 --clobber
