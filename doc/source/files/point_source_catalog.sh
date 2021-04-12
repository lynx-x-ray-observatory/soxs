#!/bin/sh
    
# Use this script to create a distribution of photons from point sources 
# and simulate an observation.

# First, make the photons from the point sources into a SIMPUT catalog,
# saving the point source properties to a table file and choosing "24"
# as a random seed to insure we get the same point sources every time.
make_point_sources my_cat.simput ptsrc 300.0,ks 20. 22.,-27.0 --random_seed=24 --output_sources=point_source_table.dat
 
# Take the SIMPUT catalog and make an event file. Since we already made a
# distribution of point sources, turn the point-source background off. 
instrument_simulator my_cat.simput ptsrc_cat_evt.fits 300.0,ks lynx_hdxi 22.,-27.0 --overwrite --no_ptsrc_bkgnd
