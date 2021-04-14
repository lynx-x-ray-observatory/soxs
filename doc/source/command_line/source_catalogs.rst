.. _cmd-source-catalogs:

Command Line Scripts for Source Catalogs
========================================

This section documents command-line scripts for generating a SIMPUT catalog
of photons from a halo catalog drawn from a cosmological simulation. 

For more information about how the source catalog generation is implemented
in SOXS, see :ref:`source-catalogs`. 

.. _cmd-make-cosmo-sources:

``make_cosmological_sources``
-----------------------------

.. code-block:: text

    usage: make_cosmological_sources [-h] [--cat_center CAT_CENTER] [--absorb_model ABSORB_MODEL] [--nh NH] [--area AREA]
                                     [--src_filename SRC_FILENAME] [--append] [--overwrite] [--output_sources OUTPUT_SOURCES]
                                     [--write_regions WRITE_REGIONS] [--random_seed RANDOM_SEED]
                                     filename name exp_time fov sky_center
    
    Create a SIMPUT photon list catalog of a cosmological background.
    
    positional arguments:
      filename              The filename of the SIMPUT catalog to be used as the root of the catalog. If it does not exist, it will be
                            created.
      name                  The name of the source in the SIMPUT catalog.
      exp_time              The exposure time to use, in seconds.
      fov                   The field of view on a side in arcminutes.
      sky_center            The center RA, Dec coordinates of the observation, in degrees, comma-separated
    
    optional arguments:
      -h, --help            show this help message and exit
      --cat_center CAT_CENTER
                            The center of the field in the coordinates of the halo catalog, which range from -5.0 to 5.0 degrees in both
                            directions. If not set, a center will be randomly chosen.
      --absorb_model ABSORB_MODEL
                            The absorption model to use for foreground galactic absorption. Default: 'wabs'
      --nh NH               The hydrogen column in units of 10**22 atoms/cm**2. Default: 0.05
      --area AREA           The collecting area to use, in cm^2. Default: 30000.0
      --src_filename SRC_FILENAME
                            An optional filename to store the source instead of the SIMPUT catalog file.
      --append              If set, append a new source an existing SIMPUT catalog.
      --overwrite           Overwrite an existing file with the same name.
      --output_sources OUTPUT_SOURCES
                            Output the source properties to the specified file.
      --write_regions WRITE_REGIONS
                            Write ds9 circle region files corresponding to the positions and r500 of the halos.
      --random_seed RANDOM_SEED
                            A constant integer random seed to produce a consistent set of random numbers.

Examples
++++++++

Generate photons from halos with a field of view of 10.0 arcminutes, to a new SIMPUT
catalog, with an exposure time of 100 ks. Let a random location in the halo catalog
be chosen:

.. code-block:: bash

    [~]$ make_cosmological_sources halos.simput halos 100.0,ks 10.0 22.0,-12.0 --overwrite

The same as before, but choose a particular location in the halo catalog:

.. code-block:: bash

    [~]$ make_cosmological_sources halos.simput halos 100.0,ks 10.0 22.0,-12.0 --cat_center=-0.1,2.0 --overwrite

Append the halo photons to an existing SIMPUT catalog, "my_cat.simput":

.. code-block:: bash

    [~]$ make_cosmological_sources my_cat.simput halos 100.0,ks 10.0 22.0,-12.0 --append

Append the halo photons to an existing SIMPUT catalog, "my_cat.simput", and 
write the source to a different file:

.. code-block:: bash

    [~]$ make_cosmological_sources my_cat.simput halos 100.0,ks 10.0 22.0,-12.0 --append --src_filename=halos.fits

Change the Galactic hydrogen column to :math:`2 \times 10^{20}~cm^{-2}`, and 
use the "tbabs" model:

.. code-block:: bash

    [~]$ make_cosmological_sources halos.simput halos 100.0,ks 10.0 22.0,-12.0 --absorb_model="tbabs" --nh=0.02 --overwrite

Write the source properties to an ASCII text file:

.. code-block:: bash

    [~]$ make_cosmological_sources halos.simput halos 100.0,ks 10.0 22.0,-12.0 --output_sources=my_halos.txt --overwrite

Write out ds9 regions corresponding to the positions and the :math:`r_{500}` of
the sources:

.. code-block:: bash

    [~]$ make_cosmological_sources halos.simput halos 100.0,ks 10.0 22.0,-12.0 --write_regions=halos.reg --overwrite

``make_point_sources``
----------------------

.. code-block:: text

    usage: make_point_sources [-h] [--absorb_model ABSORB_MODEL] [--nh NH] [--area AREA] [--src_filename SRC_FILENAME] [--append]
                              [--overwrite] [--input_sources INPUT_SOURCES] [--output_sources OUTPUT_SOURCES] [--random_seed RANDOM_SEED]
                              filename name exp_time fov sky_center
    
    Create a SIMPUT photon list catalog of a point-source background.
    
    positional arguments:
      filename              The filename of the SIMPUT catalog to be used as the root of the catalog. If it does not exist, it will be
                            created.
      name                  The name of the source in the SIMPUT catalog.
      exp_time              The exposure time to use, in seconds.
      fov                   The field of view on a side in arcminutes.
      sky_center            The center RA, Dec coordinates of the observation, in degrees, comma-separated.
    
    optional arguments:
      -h, --help            show this help message and exit
      --absorb_model ABSORB_MODEL
                            The absorption model to use for foreground galactic absorption. Default: 'wabs'
      --nh NH               The galactic hydrogen column in units of 10**22 atoms/cm**2. Default: 0.05
      --area AREA           The collecting area to use, in cm^2. Default: 30000.0
      --src_filename SRC_FILENAME
                            An optional filename to store the source instead of the SIMPUT catalog file.
      --append              If set, append a new source an existing SIMPUT catalog.
      --overwrite           Overwrite an existing file with the same name.
      --input_sources INPUT_SOURCES
                            Use a previously written table of sources as input instead of generating them.
      --output_sources OUTPUT_SOURCES
                            Output the source properties to the specified file.
      --random_seed RANDOM_SEED
                            A constant integer random seed to produce a consistent set of random numbers.
    
Examples
++++++++

Generate photons from point sources with a field of view of 5.0 arcminutes, to a new SIMPUT
catalog, with an exposure time of 75 ks:

.. code-block:: bash

    [~]$ make_point_sources pt_src.simput pt_src 75.0,ks 5.0 90.0,-10.0 --overwrite

Append the point source photons to an existing SIMPUT catalog, "my_cat":

.. code-block:: bash

    [~]$ make_point_sources my_cat.simput pt_src 75.0,ks 5.0 90.0,-10.0 --append

Append the point source photons to an existing SIMPUT catalog, "my_cat", and 
write the source to a different file:

.. code-block:: bash

    [~]$ make_point_sources my_cat.simput pt_src 75.0,ks 5.0 90.0,-10.0 --append --src_filename=sources.fits

Change the Galactic hydrogen column to :math:`3.5 \times 10^{20}~cm^{-2}`, 
and use the "tbabs" model:

.. code-block:: bash

    [~]$ make_point_sources pt_src.simput pt_src 75.0,ks 5.0 90.0,-10.0 --absorb_model="tbabs" --nh=0.035 --overwrite

Write the source properties to an ASCII text file:

.. code-block:: bash

    [~]$ make_point_sources pt_src.simput pt_src 75.0,ks 5.0 90.0,-10.0 --output_sources=my_ptsrc.txt --overwrite

Use a previously written ASCII text file of point source properties as input:

.. code-block:: bash
        
    [~]$ make_point_sources pt_src.simput pt_src 75.0,ks 5.0 90.0,-10.0 --input_sources=my_ptsrc.txt --overwrite

.. _cmd-make-point-source-list:

``make_point_source_list``
--------------------------

.. code-block:: text

    usage: make_point_source_list [-h] [--random_seed RANDOM_SEED] output_file fov sky_center
    
    Make a list of point source properties and write it to an ASCII table file.
    
    positional arguments:
      output_file           The ASCII table file to write the source properties to.
      fov                   The field of view on a side in arcminutes.
      sky_center            The center RA, Dec coordinates of the observation, in degrees, comma-separated.
    
    optional arguments:
      -h, --help            show this help message and exit
      --random_seed RANDOM_SEED
                            A constant integer random seed to produce a consistent set of random numbers.

Examples
++++++++

Generate point source properties and write them to an ASCII table, assuming a field
of view of 30 arcminutes:

.. code-block:: bash

    [~]$ make_point_source_list my_ptsrc_list.dat 30.0 90.0,-10.0
