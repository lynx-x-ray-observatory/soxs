.. _cmd-cosmology:

Command Line Scripts for Cosmological Sources
=============================================

This section documents command-line scripts for generating a SIMPUT catalog
of photons from a halo catalog drawn from a cosmological simulation. 

For more information about how the cosmological source generation is implemented
in SOXS, see :ref:`cosmology`. 

``make_cosmological_sources``
-----------------------------

usage: make_cosmological_sources [-h] [--cat_center CAT_CENTER] [--nh NH]
                                 [--area AREA] [--append] [--clobber]
                                 [--output_sources OUTPUT_SOURCES]
                                 [--random_seed RANDOM_SEED]
                                 simput_prefix phlist_prefix exp_time fov
                                 sky_center

Create a SIMPUT photon list of a cosmological background.

positional arguments:
  simput_prefix         The prefix of the SIMPUT file to be used as the root
                        of the catalog. If it does not exist, it will be
                        created.
  phlist_prefix         The prefix of the photon list file to be written.
  exp_time              The exposure time to use, in seconds.
  fov                   The field of view on a side in arcminutes.
  sky_center            The center RA, Dec coordinates of the observation, in
                        degrees, comma-separated

optional arguments:
  -h, --help            show this help message and exit
  --cat_center CAT_CENTER
                        The center of the field in the coordinates of the halo
                        catalog, which range from -5.0 to 5.0 degrees in both
                        directions. If not set, a center will be randomly
                        chosen.
  --nh NH               The hydrogen column in units of 10**22 atoms/cm**2.
                        Default: 0.05
  --area AREA           The collecting area to use, in cm^2. Default: 30000.0
  --append              If set, append a new source an existing SIMPUT
                        catalog.
  --clobber             Overwrite an existing file with the same name.
  --output_sources OUTPUT_SOURCES
                        Output the source properties to the specified file.
  --random_seed RANDOM_SEED
                        A constant integer random seed to produce a consistent
                        set of random numbers.

Examples
++++++++

Generate photons from halos with a field of view of 10.0 arcminutes, to a new SIMPUT
catalog, with an exposure time of 100 ks. Let a random location in the halo catalog
be chosen:

.. code-block:: bash

    [~]$ make_cosmological_sources halos halos 10.0 100000.0 20.0 22.0,-12.0 --clobber

The same as before, but choose a particular location in the halo catalog:

.. code-block:: bash

    [~]$ make_cosmological_sources halos halos 10.0 100000.0 20.0 22.0,-12.0 --cat_center=-0.1,2.0 --clobber

Append the halos to an existing SIMPUT catalog, "my_cat":

.. code-block:: bash

    [~]$ make_cosmological_sources my_cat halos 10.0 100000.0 20.0 22.0,-12.0 --append

Change the Galactic hydrogen column to :math:`2 \times 10^{20}~cm^{-2}`:

.. code-block:: bash

    [~]$ make_cosmological_sources halos halos 10.0 100000.0 20.0 22.0,-12.0 --nh=0.02 --clobber

Write the source properties to an ASCII text file:

.. code-block:: bash

    [~]$ make_cosmological_sources halos halos 10.0 100000.0 20.0 22.0,-12.0 --output_sources=my_halos.txt --clobber
