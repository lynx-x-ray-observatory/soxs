.. _cmd-events:

Command Line Scripts for Processing Events
==========================================

These command-line scripts allow one to make derived products from 
event files. For details on what's going on under the hood, see 
:ref:`event-tools`.

``make_exposure_map``
---------------------

This script takes an event file made by SOXS and makes a SOXS exposure map for it. 

.. code-block:: text

    usage: make_exposure_map [-h] [--energy ENERGY] [--weightsfile WEIGHTSFILE]
                             [--asol_file ASOL_FILE] [--overwrite]
                             [--nhistx NHISTX] [--nhisty NHISTY]
                             [--normalize | --no_normalize] [--reblock REBLOCK]
                             event_file expmap_file
    
    Make a SOXS exposure map from an event file.
    
    positional arguments:
      event_file            The event file to use to make the exposure map.
      expmap_file           The file to write the exposure map to.
    
    optional arguments:
      -h, --help            show this help message and exit
      --energy ENERGY       The reference energy to use when making the exposure
                            map. This parameter will be ignored if a 'weightsfile'
                            is set.
      --weightsfile WEIGHTSFILE
                            A file containing two columns: energy in keV and
                            spectral weights, to create an exposure map weighted
                            over an energy band.
      --asol_file ASOL_FILE
                            If set, write the aspect solution to this file.
      --overwrite           Overwrite an existing file with the same name.
      --nhistx NHISTX       The number of bins in the aspect histogram in the DETX
                            direction. Default: 16
      --nhisty NHISTY       The number of bins in the aspect histogram in the DETY
                            direction. Default: 16
      --normalize           Normalize the exposure map by the exposure time. This
                            is the default.
      --no_normalize        Don't normalize the exposure map by the exposure time.
      --reblock REBLOCK     Supply an integer power of two to set the binning of
                            the exposure map. Default: 1

Examples
++++++++

Make an exposure map from an event file at a single energy of 3.0 keV.

.. code-block:: bash

    [~]$ make_exposure_map evt.fits expmap.fits --energy=3.0 --overwrite

Also write an aspect solution file.

.. code-block:: bash

    [~]$ make_exposure_map evt.fits expmap.fits --energy=3.0 --asol_file=asol.fits --overwrite

Make an exposure map, using an ASCII text file with two columns, energy and flux, 
to weight the exposure. 

.. code-block:: bash

    [~]$ make_exposure_map evt.fits expmap.fits --weightsfile=spec.dat --overwrite

Make an exposure map and change the binning of the aspect histogram. 

.. code-block:: bash

    [~]$ make_exposure_map evt.fits expmap.fits --energy=3.0 --overwrite --nhistx=32 --nhisty=32

Make an exposure map, but don't normalize by the exposure time. 

.. code-block:: bash

    [~]$ make_exposure_map evt.fits expmap.fits --energy=3.0 --overwrite --no_normalize

Reblock the exposure map by 4.

.. code-block:: bash

    [~]$ make_exposure_map evt.fits expmap.fits --energy=3.0 --overwrite --reblock=4

``make_image``
--------------

.. code-block:: text

    usage: make_image [-h] [--coord_type COORD_TYPE] [--emin EMIN] [--emax EMAX]
                      [--overwrite] [--expmap_file EXPMAP_FILE]
                      [--reblock REBLOCK]
                      event_file out_file
    
    Make a FITS image from a SOXS event file.
    
    positional arguments:
      event_file            The event file to use to make the image.
      out_file              The file to write the image to.
    
    optional arguments:
      -h, --help            show this help message and exit
      --coord_type COORD_TYPE
                            The type of coordinate to bin into the image. Can be
                            'sky' or 'det'. Default: 'sky'
      --emin EMIN           The minimum energy of the photons to put in the image,
                            in keV.
      --emax EMAX           The maximum energy of the photons to put in the image,
                            in keV.
      --overwrite           Overwrite an existing file with the same name.
      --expmap_file EXPMAP_FILE
                            Supply an exposure map file to divide this image by to
                            get a flux map.
      --reblock REBLOCK     Change this value to reblock the image to larger pixel
                            sizes (reblock >= 1). Only supported for sky
                            coordinates. Default: 1

Examples
++++++++

Make an image in celestial coordinates from an event file.

.. code-block:: bash

    [~]$ make_image evt.fits img.fits --overwrite

The same image, but with a restricted energy band.

.. code-block:: bash

    [~]$ make_image evt.fits img.fits --emin=0.5 --emax=7.0 --overwrite

Make an image in detector coordinates.

.. code-block:: bash

    [~]$ make_image evt.fits det_img.fits --overwrite --coord_type=det

Make an image and divide it by an exposure map.

.. code-block:: bash

    [~]$ make_image evt.fits flux_img.fits --overwrite --expmap_file=expmap.fits

Reblock the image by 4.

.. code-block:: bash

    [~]$ make_image evt.fits img.fits --emin=0.5 --emax=7.0 --overwrite --reblock=4

``make_radial_profile``
-----------------------

.. code-block:: text

    usage: make_radial_profile [-h] [--ctr_type CTR_TYPE] [--emin EMIN]
                               [--emax EMAX] [--overwrite]
                               [--expmap_file EXPMAP_FILE]
                               event_file out_file ctr rmin rmax nbins
    
    Make a FITS radial profile from a SOXS event file.
    
    positional arguments:
      event_file            The event file to use to make the profile.
      out_file              The file to write the profile to.
      ctr                   The central coordinate of the profile. Can either be
                            in celestial coordinates (the default) or "physical"
                            pixel coordinates. If the former, the ``ctr_type``
                            keyword argument must be explicity set to "physical".
      rmin                  The minimum radius of the profile, in arcseconds.
      rmax                  The maximum radius of the profile, in arcseconds.
      nbins                 The number of bins in the profile.
    
    optional arguments:
      -h, --help            show this help message and exit
      --ctr_type CTR_TYPE   The type of center coordinate. Either 'celestial' for
                            (RA, Dec) coordinates (the default), or 'physical' for
                            pixel coordinates.
      --emin EMIN           The minimum energy of the photons to put in the profile,
                            in keV.
      --emax EMAX           The maximum energy of the photons to put in the profile,
                            in keV.
      --overwrite           Overwrite an existing file with the same name.
      --expmap_file EXPMAP_FILE
                            Supply an exposure map file to divide the profile by to
                            obtain flux-based quantities.

Examples
++++++++

Make a radial profile from an event file, using (RA, Dec) = (30.0, 45.0) as the
central coordinates of the profile. The profile runs from 0.0 arcseconds to 100.0
arcseconds, with 50 linearly spaced bins. 

.. code-block:: bash

    [~]$ make_radial_profile evt.fits profile.fits 30.0,45.0 0.0 100.0 50 --overwrite

The same profile, but with a restricted energy band.

.. code-block:: bash

    [~]$ make_radial_profile evt.fits profile.fits 30.0,45.0 0.0 100.0 50 --emin=0.2 --emax=3.0 --overwrite

The same profile, but specifying the center in physical coordinates instead.

.. code-block:: bash

    [~]$ make_radial_profile evt.fits profile.fits 1024.0,300.0 0.0 100.0 50 --ctr_type=physical --overwrite

Include an exposure map, allowing flux-based quantities to also be computed. 

.. code-block:: bash

    [~]$ make_radial_profile evt.fits profile.fits 30.0,45.0 0.0 100.0 50 --overwrite --expmap_file=expmap.fits
