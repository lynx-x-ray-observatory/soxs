.. _cmd-instrument:

Command Line Scripts for the Instrument Simulator
=================================================

These command-line scripts allow one to run and modify the instrument simulator.
For details on what's going on under the hood, see :ref:`instrument`.

``instrument_simulator``
------------------------

The ``instrument_simulator`` script takes a SIMPUT catalog and generates a 
simulated observation in a standard event file format which can then be 
processed by standard tools such as CIAO, HEATOOLS, XSPEC, etc. 

.. code-block:: text

    usage: instrument_simulator [-h] [--overwrite] [--roll_angle ROLL_ANGLE]
                                [--bkgnd_file BKGND_FILE] [--subpixel_res]
                                [--no_dither] [--dither_params DITHER_PARAMS]
                                [--aimpt_shift AIMPT_SHIFT]
                                [--random_seed RANDOM_SEED]
                                [--ptsrc_bkgnd | --no_ptsrc_bkgnd]
                                [--instr_bkgnd | --no_instr_bkgnd]
                                [--foreground | --no_foreground]
                                simput_file out_file exp_time instrument
                                sky_center
    
    Run the instrument simulator and produce a simulated event file.
    
    positional arguments:
      simput_file           The SIMPUT file to be used as input, or "None" if you
                            only want to simulate backgrounds.
      out_file              The name of the event file to be written.
      exp_time              The exposure time to use, in seconds.
      instrument            The name of the instrument to use, or alternatively
                            the name of a JSON file which contains an instrument
                            specification.
      sky_center            The center RA, Dec coordinates of the observation, in
                            degrees, comma-separated
    
    optional arguments:
      -h, --help            show this help message and exit
      --overwrite           Overwrite an existing file with the same name.
      --roll_angle ROLL_ANGLE
                            The roll angle in degrees. Default: 0.0
      --bkgnd_file BKGND_FILE
                            Use background stored in a file instead of generating
                            one.
      --subpixel_res        Don't uniformly distribute event positions within
                            pixels.
      --no_dither           Turn dithering off entirely.
      --dither_params DITHER_PARAMS
                            The parameters controlling the size and period of
                            dither. Four floats joined by commas, in the form of
                            x_amp,y_amp,x_period,y_period. The first two numbers
                            are in arcseconds and the second are in seconds.
                            Default: 8.0,8.0,1000.0,707.0
      --aimpt_shift AIMPT_SHIFT
                            The shift of the aimpoint on the detector in both
                            directions from the nominal aimpoint in arcseconds.
                            Default: [0.0, 0.0]
      --random_seed RANDOM_SEED
                            A constant integer random seed to produce a consistent
                            set of random numbers.
      --ptsrc_bkgnd         Turn the point-source background on.
      --no_ptsrc_bkgnd      Turn the point-source background off.
      --instr_bkgnd         Turn the instrumental background on.
      --no_instr_bkgnd      Turn the instrumental background off.
      --foreground          Turn the galactic foreground on.
      --no_foreground       Turn the galactic foreground off.

Examples
++++++++

Changing Instrument Specification
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This example uses the pre-built HDXI instrument specification, assuming a 50 ks observation
with the pointing (RA, Dec) = (30, 45) degrees.

.. code-block:: bash

    [~]$ instrument_simulator sloshing_simput.fits evt.fits 50.0,ks hdxi 30.,45. --overwrite

The same, but use the HDXI specification with mirror diameter of :math:`d` = 3 m and focal length of
:math:`f` = 20 m:

.. code-block:: bash

    [~]$ instrument_simulator sloshing_simput.fits evt.fits 50.0,ks hdxi_3x20 30.,45. --overwrite

See :ref:`instrument-arg` for details on the options for the ``instrument`` argument.

This example uses a JSON file created by the user, which contains a custom instrument specification. See
:ref:`instrument-registry` for details on how to do this.

.. code-block:: bash

    [~]$ instrument_simulator sloshing_simput.fits evt.fits 50.0,ks my_inst.json 30.,45. --overwrite

The following details how to change the other options, for more info see :ref:`other-mods`.

Changing Roll Angle, Dither, and Aimpoint Shift
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Change the roll angle to 45 degrees:

.. code-block:: bash

    [~]$ instrument_simulator sloshing_simput.fits evt.fits 50.0,ks hdxi 30.,45. --roll_angle=45.0 --overwrite

Change the dither amplitudes to 32 arcseconds and the periods to 707 and 1200 seconds:

.. code-block:: bash

    [~]$ instrument_simulator sloshing_simput.fits evt.fits 50.0,ks hdxi 30.,45. --dither_params=32.,32.,707.,1200. --overwrite

Turn dither off entirely:

.. code-block:: bash

    [~]$ instrument_simulator sloshing_simput.fits evt.fits 50.0,ks hdxi 30.,45. --no_dither --overwrite

Shift the aimpoint (units are in arcseconds):

.. code-block:: bash

    [~]$ instrument_simulator sloshing_simput.fits evt.fits 50.0,ks hdxi 30.,45. --aimpt_shift=10.0,-15.0 --overwrite

Customizing Backgrounds
~~~~~~~~~~~~~~~~~~~~~~~

Turn off the instrumental background:

.. code-block:: bash

    [~]$ instrument_simulator sloshing_simput.fits evt.fits 50.0,ks hdxi 30.,45. --no_instr_bkgnd --overwrite

Turn off the Galactic foreground:

.. code-block:: bash

    [~]$ instrument_simulator sloshing_simput.fits evt.fits 50.0,ks hdxi 30.,45. --no_foreground --overwrite

Turn off the point-source background:

.. code-block:: bash

    [~]$ instrument_simulator sloshing_simput.fits evt.fits 50.0,ks hdxi 30.,45. --no_ptsrc_bkgnd --overwrite

Any combination of these may be used to turn multiple components off or all 
of them. 

To use a background stored in an event file:

.. code-block:: bash

    [~]$ instrument_simulator sloshing_simput.fits evt.fits 50.0,ks hdxi 30.,45. --bkgnd_file="bkg_evt.fits" --overwrite

.. note::

    If you use a background stored in an event file, the background will be 
    entirely determined from the contents of this file and any of the above
    background flags will be ignored.

.. _cmd-simulate-spectrum:

``simulate_spectrum``
---------------------

Generate a PI or PHA spectrum from a spectrum in an ASCII table (such as 
one made by one of the commands detailed in :ref:`cmd-spectra`) by convolving
it with responses. To be used if one wants to create a spectrum without 
worrying about spatial response, or if the underlying instrument supports
only simulating spectra. Similar to XSPEC's "fakeit". 

.. code-block:: bash

    usage: simulate_spectrum [-h] [--overwrite] [--bkgnd_area BKGND_AREA]
                             [--random_seed RANDOM_SEED]
                             [--ptsrc_bkgnd | --no_ptsrc_bkgnd]
                             [--instr_bkgnd | --no_instr_bkgnd]
                             [--foreground | --no_foreground]
                             spec_file instrument exp_time out_file
    
    Convolve a spectrum with an ARF and RMF and produce a PHA or PI spectrum.
    
    positional arguments:
      spec_file             The file containing the spectrum to be used. If None,
                            then only a simulated background may be generated if
                            they are turned on.
      instrument            The name of the instrument to use, or alternatively
                            the name of a JSON file which contains an instrument
                            specification.
      exp_time              The exposure time to use, in seconds.
      out_file              The file to write the convolved spectrum to.
    
    optional arguments:
      -h, --help            show this help message and exit
      --overwrite           Overwrite an existing file with the same name.
      --bkgnd_area BKGND_AREA
                            The area on the sky for the background components, in
                            square arcminutes. Default: None. Must be specified if
                            any of the background components are turned on.
      --random_seed RANDOM_SEED
                            A constant integer random seed to produce a consistent
                            set of random numbers.
      --ptsrc_bkgnd         Turn the unresolved point-source background on.
      --no_ptsrc_bkgnd      Turn the unresolved point-source background off.
      --instr_bkgnd         Turn the instrumental background on.
      --no_instr_bkgnd      Turn the instrumental background off.
      --foreground          Turn the galactic foreground on.
      --no_foreground       Turn the galactic foreground off.

Examples
++++++++

Simulate a Lynx microcalorimeter spectrum.

.. code-block:: bash

    [~]$ simulate_spectrum power_law_spec.dat lynx_lxm 300.0,ks plaw_spec.pha

The same spectrum, but with point-source, foreground, and instrumental backgrounds
added. Two square arcminutes of background assumed. 

.. code-block:: bash

    [~]$ simulate_spectrum power_law_spec.dat lynx_lxm 300.0,ks plaw_spec.pha --bkgnd_area 2.0 --ptsrc_bkgnd --foreground --instr_bkgnd

Simulate backgrounds only.

.. code-block:: bash

    [~]$ simulate_spectrum None lynx_lxm 300.0,ks plaw_spec.pha --bkgnd_area 2.0 --ptsrc_bkgnd --foreground --instr_bkgnd

.. _cmd-get-instrument-data:

``get_instrument_data``
-----------------------

Download the files needed for a particular instrument to a location of one's 
choosing.

.. code-block:: bash

    usage: get_instrument_data [-h] [--loc LOC] instrument
    
    Download files associated with a particular instrument model.
    
    positional arguments:
      instrument  The name of the instrument to download files for.
    
    optional arguments:
      -h, --help  show this help message and exit
      --loc LOC   The path to download the files to. Defaults to the current
                  working directory.

Examples
++++++++

Get the instrument files for the ``"lynx_hdxi"`` instrument model, saved to
the current working directory. 

.. code-block:: bash

    [~]$ get_instrument_data lynx_hdxi

Get the instrument files for the ``"xrism_resolve"`` instrument model, saved
to a particular directory.

.. code-block:: bash

    [~]$ get_instrument_data lynx_hdxi --loc /Users/jzuhone/Data/soxs
