.. _cmd-background:

Command Line Scripts for Generating Backgrounds
===============================================

``make_background_file``
------------------------

.. code-block:: text

    usage: make_background_file [-h] [--clobber]
                                [--dither_shape DITHER_SHAPE]
                                [--dither_size DITHER_SIZE]
                                [--random_seed RANDOM_SEED]
                                [--ptsrc_bkgnd | --no_ptsrc_bkgnd]
                                [--instr_bkgnd | --no_instr_bkgnd]
                                [--foreground | --no_foreground]
                                simput_file out_file exp_time instrument
                                sky_center
    
    Run the instrument simulator and produce a simulated background event file.
    
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
      --clobber             Overwrite an existing file with the same name.
      --dither_shape DITHER_SHAPE
                            The shape of the dither pattern: square, circle, or
                            None. Default: square
      --dither_size DITHER_SIZE
                            The size of the dither pattern in arcseconds. For a
                            circle, thesize is the radius; for a square, the size
                            is the width. Default: 16.0
      --random_seed RANDOM_SEED
                            A constant integer random seed to produce a consistent
                            set of random numbers.
      --ptsrc_bkgnd         Turn the point-source background on.
      --no_ptsrc_bkgnd      Turn the point-source background off.
      --instr_bkgnd         Turn the instrumental background on.
      --no_instr_bkgnd      Turn the instrumental background off.
      --foreground          Turn the galactic foreground on.
      --no_foreground       Turn the galactic foreground off.
