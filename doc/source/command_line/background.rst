.. _cmd-background:

Command Line Scripts for Generating Backgrounds
===============================================

These command line scripts allow one to generate background event files.

``make_background_file``
------------------------

The ``make_background_file`` generates a simulated observation of background
in a standard event file format which can then be used as the background for
and observation or processed by standard tools such as CIAO, HEATOOLS, XSPEC,
etc.

.. code-block:: text

    usage: make_background_file [-h] [--overwrite] [--input_pt_sources INPUT_PT_SOURCES]
                                [--no_diffuse_unresolved]
                                [--drop_brightest DROP_BRIGHTEST] [--subpixel_res]
                                [--no_dither] [--dither_params DITHER_PARAMS]
                                [--random_seed RANDOM_SEED] [--ptsrc_bkgnd |
                                --no_ptsrc_bkgnd] [--instr_bkgnd | --no_instr_bkgnd]
                                [--foreground | --no_foreground]
                                [--instr_bkgnd_scale INSTR_BKGND_SCALE]
                                out_file exp_time instrument sky_center

    Run the instrument simulator and produce a simulated background event file.

    positional arguments:
      out_file              The name of the event file to be written.
      exp_time              The exposure time to use, in seconds.
      instrument            The name of the instrument to use, or alternatively the name
                            of a JSON file which contains an instrument specification.
      sky_center            The center RA, Dec coordinates of the observation, in degrees,
                            comma-separated

    options:
      -h, --help            show this help message and exit
      --overwrite           Overwrite an existing file with the same name.
      --input_pt_sources INPUT_PT_SOURCES
                            Use a previously written table of point sources as input
                            instead of generating them.
      --no_diffuse_unresolved
                            If set, the diffuse component across the entire field of view
                            to represent the unresolved flux from sources at very small
                            fluxes will be turned off.
      --drop_brightest DROP_BRIGHTEST
                            This many brightest sources will be dropped from the point
                            sources.
      --subpixel_res        Don't uniformly distribute event positions within pixels.
      --no_dither           Turn dithering off entirely.
      --dither_params DITHER_PARAMS
                            The parameters controlling the size and period of dither. Four
                            floats joined by commas, in the form of
                            x_amp,y_amp,x_period,y_period. The first two numbers are in
                            arcseconds and the second are in seconds. Default:
                            8.0,8.0,1000.0,707.0
      --random_seed RANDOM_SEED
                            A constant integer random seed to produce a consistent set of
                            random numbers.
      --ptsrc_bkgnd         Turn the point-source background on.
      --no_ptsrc_bkgnd      Turn the point-source background off.
      --instr_bkgnd         Turn the instrumental background on.
      --no_instr_bkgnd      Turn the instrumental background off.
      --foreground          Turn the galactic foreground on.
      --no_foreground       Turn the galactic foreground off.
      --instr_bkgnd_scale INSTR_BKGND_SCALE
                            A scaling factor for the instrumental background. Default: 1.0.

Examples
++++++++

Changing Instrument Specification
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This example uses the pre-built HDXI instrument specification, assuming a 50 ks observation
with the pointing (RA, Dec) = (30, 45) degrees.

.. code-block:: bash

    [~]$ make_background_file bkg_evt.fits 50.0,ks hdxi 30.,45. --overwrite

The same, but use the HDXI specification with mirror diameter of :math:`d` = 3 m and focal length of
:math:`f` = 20 m:

.. code-block:: bash

    [~]$ make_background_file bkg_evt.fits 50.0,ks hdxi_3x20 30.,45. --overwrite

See :ref:`instrument-arg` for details on the options for the ``instrument`` argument.

This example uses a JSON file created by the user, which contains a custom instrument specification. See
:ref:`instrument-registry` for details on how to do this.

.. code-block:: bash

    [~]$ make_background_file bkg_evt.fits 50.0,ks my_inst.json 30.,45. --overwrite

Changing Dither
~~~~~~~~~~~~~~~

Change the dither amplitudes to 32 arcseconds and the periods to 707 and 1200 seconds:

.. code-block:: bash

    [~]$ make_background_file bkg_evt.fits 50.0,ks hdxi 30.,45. --dither_params=32.,32.,707.,1200. --overwrite

Turn dither off entirely:

.. code-block:: bash

    [~]$ make_background_file bkg_evt.fits 50.0,ks hdxi 30.,45. --no_dither --overwrite

Customizing the Background
~~~~~~~~~~~~~~~~~~~~~~~~~~

Turn off the instrumental background:

.. code-block:: bash

    [~]$ make_background_file bkg_evt.fits 50.0,ks hdxi 30.,45. --no_instr_bkgnd --overwrite

Turn off the Galactic foreground:

.. code-block:: bash

    [~]$ make_background_file bkg_evt.fits 50.0,ks hdxi 30.,45. --no_foreground --overwrite

Turn off the point-source background:

.. code-block:: bash

    [~]$ make_background_file bkg_evt.fits 50.0,ks hdxi 30.,45. --no_ptsrc_bkgnd --overwrite

Any combination of these may be used to turn multiple components off or all
of them.

Use a pre-made ASCII table of point-source properties to generate the point-source background:

.. code-block:: bash

    [~]$ make_background_file bkg_evt.fits 50.0,ks hdxi 30.,45. --input_sources=my_ptsrc.dat --overwrite

Turn off the diffuse unresolved CXB:

.. code-block:: bash

    [~]$ make_background_file bkg_evt.fits 50.0,ks hdxi 30.,45. --no_diffuse_unresolved --overwrite

Drop the 50 brightest point sources:

.. code-block:: bash

    [~]$ make_background_file bkg_evt.fits 50.0,ks hdxi 30.,45. --drop_brightest=50 --overwrite
