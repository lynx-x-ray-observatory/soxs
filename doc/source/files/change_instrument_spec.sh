#!/bin/bash

# First, we make a SIMPUT file of gas sloshing in a galaxy cluster core using the Python 
# script which calls pyXSIM
python make_sloshing.py

# Next, we make three event files, using a different instrument specification for each

# Normal HDXI with 0.5 arcsec PSF
instrument_simulator sloshing_simput.fits evt.fits 50.0,ks lynx_hdxi 30.0,45.0 --overwrite

# HDXI with 2 arcsec PSF
instrument_simulator sloshing_simput.fits evt_2.fits 50.0,ks hdxi_2.json 30.0,45.0 --overwrite

# HDXI with 5 arcsec PSF
instrument_simulator sloshing_simput.fits evt_5.fits 50.0,ks hdxi_5.json 30.0,45.0 --overwrite
