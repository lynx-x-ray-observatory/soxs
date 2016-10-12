.. _make-event-file:

``make_event_file``
===================

The ``make_event_file`` script takes a SIMPUT file and generates a simulated observation
in a standard event file format which can then be processed by standard tools such as 
CIAO, HEATOOLS, XSPEC, etc. 

.. code-block:: text

    usage: make_event_file [-h] [--clobber]
                           simput_file out_file exp_time instrument sky_center
    
    Create a simulated event file.
    
    positional arguments:
      simput_file  The SIMPUT file to be used as input.
      out_file     The name of the event file to be written.
      exp_time     The exposure time to use, in seconds.
      instrument   The name of the instrument to use, or alternatively the name of
                   a JSON file which contains an instrument specification.
      sky_center   The center RA, Dec coordinates of the observation, in degrees,
                   comma-separated
    
    optional arguments:
      -h, --help   show this help message and exit
      --clobber    Whether or not to clobber an existing file with the same name.