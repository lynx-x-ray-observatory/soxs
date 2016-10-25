.. _cmd-spatial:

Command Line Scripts for Spatial Models
=======================================

These are scripts that create photon lists for SIMPUT catalogs that can serve
as inputs to the instrument simulator.

``make_point_source``
---------------------

.. code-block:: text

    usage: make_point_source [-h] [--area AREA] [--clobber]
                             simput_prefix phlist_prefix ra dec specfile exp_time
    
    Create a SIMPUT photon list of a point source froma spectrum supplied in a
    file.
    
    positional arguments:
      simput_prefix  The prefix of the SIMPUT file to be used as the root of the
                     catalog. If it does not exist, it will be created.
      phlist_prefix  The prefix of the photon list file to be written.
      ra             The right ascension of the source in degrees.
      dec            The declination of the source in degrees.
      specfile       The file containing the spectrum to be used.
      exp_time       The exposure time to use, in seconds.
    
    optional arguments:
      -h, --help     show this help message and exit
      --area AREA    The collecting area to use, in cm^2.
      --clobber      Whether or not to clobber an existing file with the same
                     name.

Examples
++++++++

``make_beta_model``
-------------------

.. code-block:: text

    usage: make_beta_model [-h] [--velocity VELOCITY] [--absorb] [--nh NH]
                           [--area AREA] [--clobber]
                           simput_prefix phlist_prefix ra0 dec0 r_c beta kT abund
                           redshift flux emin emax exp_time
    
    Create a SIMPUT photon list of an isothermal beta-model source.
    
    positional arguments:
      simput_prefix        The prefix of the SIMPUT file to be used as the root of
                           the catalog. If it does not exist, it will be created.
      phlist_prefix        The prefix of the photon list file to be written.
      ra0                  The right ascension of the source center in degrees.
      dec0                 The declination of the source center in degrees.
      r_c                  The core radius in arcseconds.
      beta                 The beta parameter.
      kT                   The temperature in keV.
      abund                The metallicity in solar units.
      redshift             The temperature in keV.
      flux                 The total flux in units of erg/cm**2/s.
      emin                 The lower reference energy in keV.
      emax                 The upper reference energy in keV.
      exp_time             The exposure time to use, in seconds.
    
    optional arguments:
      -h, --help           show this help message and exit
      --velocity VELOCITY  The velocity broadening parameter, in units of km/s.
      --absorb             Whether or not to apply foreground galactic absorption.
      --nh NH              The hydrogen column in units of 10**22 atoms/cm**2
      --area AREA          The collecting area to use, in cm^2.
      --clobber            Whether or not to clobber an existing file with the
                           same name.

Examples
++++++++

``make_astrophysical_background``
---------------------------------

.. code-block:: text

    usage: make_astrophysical_background [-h] [--bkgnd_file BKGND_FILE]
                                         [--area AREA] [--clobber]
                                         simput_prefix phlist_prefix ra_pnt
                                         dec_pnt fov exp_time
    
    Create a SIMPUT photon list of a point source froma spectrum supplied in a
    file.
    
    positional arguments:
      simput_prefix         The prefix of the SIMPUT file to be used as the root
                            of the catalog. If it does not exist, it will be
                            created.
      phlist_prefix         The prefix of the photon list file to be written.
      ra_pnt                The right ascension of the source in degrees.
      dec_pnt               The declination of the source in degrees.
      fov                   The width of the field of view in arcminutes.
      exp_time              The exposure time to use, in seconds.
    
    optional arguments:
      -h, --help            show this help message and exit
      --bkgnd_file BKGND_FILE
                            The file containing the spectrum to be used to make
                            the background.
      --area AREA           The collecting area to use, in cm^2.
      --clobber             Whether or not to clobber an existing file with the
                            same name.

Examples
++++++++
