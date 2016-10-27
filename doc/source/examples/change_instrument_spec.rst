.. _change-instrument-spec:

Changing the Instrument Specification
=====================================

This script shows how to simulate a single source with different instrument specifications, 
including JSON files. We will make a small tweak to HDXI, by increasing the PSF size. This 
example requires pyXSIM.

.. code-block:: bash

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

The ``make_sloshing.py`` script that is called:

.. code-block:: python

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
    redshift = 0.1
    
    # This line generates the photons
    photons = pyxsim.PhotonList.from_data_source(sp, redshift, area, exp_time, source_model)
    
    # Here we define the foreground absorption using the tbabs model
    tbabs_model = pyxsim.TBabsModel(0.04)
    
    # Project the photons along the z-axis and absorb them
    events = photons.project_photons("z", absorb_model=tbabs_model, sky_center=(30.,45.))
    
    # Write the events to a SIMPUT catalog
    events.write_simput_file("sloshing", clobber=True)
    
Download these scripts and JSON files here: 

* `change_instrument_spec.sh <../change_instrument_spec.sh>`_
* `make_sloshing.py <../make_sloshing.py>`_
* `hdxi_2.json <../hdxi_2.json>`_
* `hdxi_5.json <../hdxi_5.json>`_

The result of this script are three event files with an increasingly broadened PSF, as shown
here: