.. _cmd-instrument:

Command Line Scripts for the Instrument Simulator
=================================================

These command-line scripts allow one to run and modify the instrument simulator. For details on
what's going on under the hood, see :ref:`instrument`.

``instrument_simulator``
------------------------

The ``instrument_simulator`` script takes a SIMPUT catalog and generates a simulated observation
in a standard event file format which can then be processed by standard tools such as CIAO, 
HEATOOLS, XSPEC, etc. 

.. code-block:: text

    usage: instrument_simulator [-h] [--clobber] [--dither_shape DITHER_SHAPE]
                                [--dither_size DITHER_SIZE]
                                [--roll_angle ROLL_ANGLE]
                                [--astro_bkgnd ASTRO_BKGND]
                                [--instr_bkgnd_scale INSTR_BKGND_SCALE]
                                simput_file out_file exp_time instrument
                                sky_center
    
    Run the instrument simulator and produce a simulated event file.
    
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
                            The shape of the dither pattern: square, circle, or
                            None. Default: square
      --dither_size DITHER_SIZE
                            The size of the dither pattern in arcseconds. For a
                            circle, thesize is the radius; for a square, the size
                            is the width. Default: 16.0
      --roll_angle ROLL_ANGLE
                            The roll angle in degrees. Default: 0.0
      --astro_bkgnd ASTRO_BKGND
                            The astrophysical background to use. Default: hm_cxb
      --instr_bkgnd_scale INSTR_BKGND_SCALE
                            The scale of the instrumental background. Default: 1.0

Examples
++++++++

Changing Instrument Specification
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This example uses the pre-built HDXI instrument specification, assuming a 50 ks observation
with the pointing (RA, Dec) = (30, 45) degrees.

.. code-block:: bash

    [~]$ instrument_simulator sloshing_simput.fits evt.fits 50000.0 hdxi 30.,45. --clobber

The same, but use the HDXI specification with mirror diameter of :math:`d` = 3 m and focal length of
:math:`f` = 20 m:

.. code-block:: bash

    [~]$ instrument_simulator sloshing_simput.fits evt.fits 50000.0 hdxi_3x20 30.,45. --clobber

This example uses a JSON file created by the user, which contains a custom instrument specification. See
:ref:`instrument-registry` for details on how to do this.

.. code-block:: bash

    [~]$ instrument_simulator sloshing_simput.fits evt.fits 50000.0 my_inst.json 30.,45. --clobber

Changing Roll Angle and Dither
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Change the roll angle to 45 degrees:

.. code-block:: bash

    [~]$ instrument_simulator sloshing_simput.fits evt.fits 50000.0 hdxi 30.,45. --roll_angle=45.0 --clobber

Change the dither shape to a circle and make the dither radius 32 arcsec:

.. code-block:: bash

    [~]$ instrument_simulator sloshing_simput.fits evt.fits 50000.0 hdxi 30.,45. --dither_shape=circle --dither_size=32.0 --clobber

Turn dithering off entirely:

.. code-block:: bash

    [~]$ instrument_simulator sloshing_simput.fits evt.fits 50000.0 hdxi 30.,45. --dither_shape=None --clobber

Changing Backgrounds
~~~~~~~~~~~~~~~~~~~~

Rescale the instrumental background by 1/3:

.. code-block:: bash

    [~]$ instrument_simulator sloshing_simput.fits evt.fits 50000.0 hdxi 30.,45. --instr_bkgnd_scale=0.33333333 --clobber

Set the astrophysical background to the "my_bkg" model, which must be in the 
`background registry <../python/background.html>`_:

.. code-block:: bash

    [~]$ instrument_simulator sloshing_simput.fits evt.fits 50000.0 hdxi 30.,45. --astro_bkgnd=my_bkg --clobber

Turn the astrophysical background off entirely:

.. code-block:: bash

    [~]$ instrument_simulator sloshing_simput.fits evt.fits 50000.0 hdxi 30.,45. --astro_bkgnd=None --clobber
