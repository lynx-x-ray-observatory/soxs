#!/bin/bash

# We will use this script to illustrate how to change the dithering options. 

# First, make an absorbed power-law spectrum, which will be for a point source.
make_powerlaw_spectrum 1.1 0.05 1.0e-4 powerlaw_spec.dat --absorb --nh 0.04 --clobber

# Take this spectrum and make a SIMPUT catalog with a point source photon list
make_point_source my_cat point_source 30.0 45.0 powerlaw_spec.dat 100000.0 --clobber

# Next, we make three event files, using a different dither for each

# Normal HDXI with the default 16.0 arcsec square dither
instrument_simulator my_cat_simput.fits evt_square.fits 50000.0 hdxi 30.0,45.0 --clobber

# HDXI with a 16.0 arcsec radius circle dither
instrument_simulator my_cat_simput.fits evt_circle.fits 50000.0 hdxi 30.0,45.0 --dither_shape=circle --dither_size=16.0 --clobber

# HDXI with no dither
instrument_simulator my_cat_simput.fits evt_none.fits 50000.0 hdxi 30.0,45.0 --dither_shape=None --clobber
