.. _cmd-events:

Command Line Scripts for Processing Events
==========================================

These command-line scripts allow one to make derived products from 
event files. For details on what's going on under the hood, see 
:ref:`event-tools`.

``make_image``
--------------

.. code-block:: text

    usage: make_image [-h] [--coord_type COORD_TYPE] [--emin EMIN] [--emax EMAX]
                      [--overwrite] [--expmap_file EXPMAP_FILE]
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
