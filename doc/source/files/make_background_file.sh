#!/bin/sh

# Use this script to make a background event file using a pre-made 
# set of point sources.

# First, make the point source properties using the make_point_source_list
# command which saves the positions, fluxes, and spectral indices of the
# sources to an ASCII table which can be read in by make_background_file.
make_point_source_list point_source_table.dat 20. 22.,-27.0

# Take the SIMPUT catalog and make an event file. Since we already made a
# distribution of point sources, turn the point-source background off. 
make_background_file bkgnd_evt.fits 300.0,ks lynx_hdxi 22.,-27.0 --overwrite --input_pt_sources=point_source_table.dat
