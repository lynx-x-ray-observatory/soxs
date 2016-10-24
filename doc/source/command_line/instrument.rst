.. _cmd-events:

Instrument Simulator Scripts
============================

These command-line scripts allow one to run and modify the instrument simulator. 

``instrument_simulator``
------------------------

The ``instrument_simulator`` script takes a SIMPUT file and generates a simulated observation
in a standard event file format which can then be processed by standard tools such as CIAO, 
HEATOOLS, XSPEC, etc. 

.. code-block:: text

    usage: instrument_simulator [-h] [--clobber] [--dither_shape DITHER_SHAPE]
                                [--dither_size DITHER_SIZE] [--roll_angle ROLL_ANGLE]
                                [--bkgnd_scale BKGND_SCALE]
                                simput_file out_file exp_time instrument sky_center
    
    Create a simulated event file.
    
    positional arguments:
      simput_file           The SIMPUT file to be used as input.
      out_file              The name of the event file to be written.
      exp_time              The exposure time to use, in seconds.
      instrument            The name of the instrument to use, or alternatively
                            the name of a JSON file which contains an instrument
                            specification.
      sky_center            The center RA, Dec coordinates of the observation, in
                            degrees, comma-separated
    
    optional arguments:
      -h, --help            show this help message and exit
      --clobber             Whether or not to clobber an existing file with the
                            same name.
      --dither_shape DITHER_SHAPE
                            The shape of the dither pattern: "square", "circle",
                            or "None". Default: "square"
      --dither_size DITHER_SIZE
                            The size of the dither pattern in arcseconds. For a
                            circle, thesize is the radius; for a square, the size
                            is the width. Default: 16.0
      --roll_angle ROLL_ANGLE
                            The roll angle in degrees. Default: 0.0
      --bkgnd_scale BKGND_SCALE
                            The scale of the background.

Examples
++++++++

This example uses the pre-built HDXI instrument specification. 

.. code-block:: bash

    [~]$ instrument_simulator sloshing_simput.fits evt.fits 50000.0 hdxi 30.,45. --clobber

This example uses a JSON file created by the user. 

.. code-block:: bash

    [~]$ instrument_simulator sloshing_simput.fits evt.fits 50000.0 my_inst.json 30.,45. --clobber
