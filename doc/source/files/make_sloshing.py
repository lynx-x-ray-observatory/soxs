import yt
import pyxsim

# Load the file using yt
# This file can be obtained at http://yt-project.org/data/GasSloshing.tar.gz
ds = yt.load("GasSloshing/sloshing_nomag2_hdf5_plt_cnt_0150")

# We create a sphere at the center of the simulation domain with a radius of 500 kpc
sp = ds.sphere("c", (500.,"kpc"))

# We define our emission model to be a thermal model using the APEC tables.
# Metallicity is the same everywhere with Z = 0.3 solar
spec_model = pyxsim.TableApecModel(0.05, 11.0, 10000, thermal_broad=True)
source_model = pyxsim.ThermalSourceModel(spec_model, Zmet=0.3)

# We set up some basic parameters to determine the sample
exp_time = (50., "ks") # exposure time
area = (30000.0, "cm**2") # collecting area
redshift = 0.2

# This line generates the photons
photons = pyxsim.PhotonList.from_data_source(sp, redshift, area, exp_time, source_model)

# Here we define the foreground absorption using the tbabs model
tbabs_model = pyxsim.TBabsModel(0.04)

# Project the photons along the z-axis and absorb them
events = photons.project_photons("z", absorb_model=tbabs_model, sky_center=(30.,45.))

# Write the events to a SIMPUT catalog
events.write_simput_file("sloshing", clobber=True)
