.. _cmd-spatial:

Command Line Scripts for Spatial Models
=======================================

These are scripts that create photon lists for SIMPUT catalogs that can serve
as inputs to the instrument simulator.

``make_point_source``
---------------------

.. code-block:: text

    usage: make_point_source [-h] [--area AREA] [--append] [--clobber]
                             simput_prefix phlist_prefix ra0 dec0 specfile
                             exp_time
    
    Create a SIMPUT photon list of a point source from a spectrum supplied in a
    file.
    
    positional arguments:
      simput_prefix  The prefix of the SIMPUT file to be used as the root of the
                     catalog. If it does not exist, it will be created.
      phlist_prefix  The prefix of the photon list file to be written.
      ra0            The right ascension of the source in degrees.
      dec0           The declination of the source in degrees.
      specfile       The file containing the spectrum to be used.
      exp_time       The exposure time to use, in seconds.
    
    optional arguments:
      -h, --help     show this help message and exit
      --area AREA    The collecting area to use, in cm^2. Default: 30000.0
      --append       If set, append a new source an existing SIMPUT catalog.
      --clobber      Whether or not to clobber an existing file with the same
                     name.
                     
Examples
++++++++

``make_beta_model``
-------------------

.. code-block:: text

    usage: make_beta_model [-h] [--area AREA] [--append] [--clobber]
                           simput_prefix phlist_prefix ra0 dec0 r_c beta specfile
                           exp_time
    
    Create a SIMPUT photon list of a beta-model source.
    
    positional arguments:
      simput_prefix  The prefix of the SIMPUT file to be used as the root of the
                     catalog. If it does not exist, it will be created.
      phlist_prefix  The prefix of the photon list file to be written.
      ra0            The right ascension of the source center in degrees.
      dec0           The declination of the source center in degrees.
      r_c            The core radius in arcseconds.
      beta           The beta parameter.
      specfile       The file containing the spectrum to be used.
      exp_time       The exposure time to use, in seconds.
    
    optional arguments:
      -h, --help     show this help message and exit
      --area AREA    The collecting area to use, in cm^2. Default: 30000.0
      --append       If set, append a new source an existing SIMPUT catalog.
      --clobber      Whether or not to clobber an existing file with the same
                     name.
                     
Examples
++++++++

``make_annulus_source``
-----------------------

.. code-block:: text

    usage: make_annulus_source [-h] [--area AREA] [--append] [--clobber]
                               simput_prefix phlist_prefix ra0 dec0 r_in r_out
                               specfile exp_time
    
    Create a SIMPUT photon list of an annulus source with uniform surface
    brightness from a spectrum supplied in a file.
    
    positional arguments:
      simput_prefix  The prefix of the SIMPUT file to be used as the root of the
                     catalog. If it does not exist, it will be created.
      phlist_prefix  The prefix of the photon list file to be written.
      ra0            The right ascension of the source center in degrees.
      dec0           The declination of the source center in degrees.
      r_in           The inner annulus of the source center in arcseconds.
      r_out          The outer annulus of the source center in arcseconds.
      specfile       The file containing the spectrum to be used.
      exp_time       The exposure time to use, in seconds.
    
    optional arguments:
      -h, --help     show this help message and exit
      --area AREA    The collecting area to use, in cm^2. Default: 30000.0
      --append       If set, append a new source an existing SIMPUT catalog.
      --clobber      Whether or not to clobber an existing file with the same
                     name.

Examples
++++++++

``make_fov_source``
-------------------

.. code-block:: text

    usage: make_fov_source [-h] [--area AREA] [--append] [--clobber]
                           simput_prefix phlist_prefix ra0 dec0 fov specfile
                           exp_time
    
    Create a SIMPUT photon list of a uniformly filled field of view source from a
    spectrum supplied in a file.
    
    positional arguments:
      simput_prefix  The prefix of the SIMPUT file to be used as the root of the
                     catalog. If it does not exist, it will be created.
      phlist_prefix  The prefix of the photon list file to be written.
      ra0            The right ascension of the source center in degrees.
      dec0           The declination of the source center in degrees.
      fov            The field of view on a side in arcminutes.
      specfile       The file containing the spectrum to be used.
      exp_time       The exposure time to use, in seconds.
    
    optional arguments:
      -h, --help     show this help message and exit
      --area AREA    The collecting area to use, in cm^2. Default: 30000.0
      --append       If set, append a new source an existing SIMPUT catalog.
      --clobber      Whether or not to clobber an existing file with the same
                     name.

Examples
++++++++
