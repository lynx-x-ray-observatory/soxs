import pyxsim
import yt

# This example simulates photons for a sloshing galaxy cluster core.

# Load the file using yt
# This file can be obtained at http://yt-project.org/data/GasSloshing.tar.gz
ds = yt.load("GasSloshing/sloshing_nomag2_hdf5_plt_cnt_0150")

# We create a sphere at the center of the simulation domain with a radius of 500 kpc
sp = ds.sphere("c", (500.0, "kpc"))

# We define our emission model to be a thermal model using the APEC tables.
# Metallicity is the same everywhere with Z = 0.3 solar
source_model = pyxsim.CIESourceModel("apec", 0.05, 11.0, 1000, 0.3, thermal_broad=True)

# We set up some basic parameters to determine the sample
exp_time = (50.0, "ks")  # exposure time
area = (3.0, "m**2")  # collecting area
redshift = 0.1

# This line generates the photons
n_photons, n_cells = pyxsim.make_photons(
    "my_photons", sp, redshift, area, exp_time, source_model
)

# Project the photons along the z-axis and absorb them
n_events = pyxsim.project_photons(
    "my_photons", "my_events", "z", (30.0, 45.0), absorb_model="tbabs", nH=0.04
)

# Write the events to a SIMPUT catalog
events = pyxsim.EventList("my_events.h5")
events.write_to_simput("sloshing", overwrite=True)
