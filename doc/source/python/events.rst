.. _events:

Creating Event Files
====================

Creating Events from SIMPUT Sources
-----------------------------------

The end product of a mock observation is a "standard" event file which has been 
convolved with a model for the telescope. In SOX, this is handled by the
instrument simulator. 

:func:`~sox.events.make_event_file` reads in a SIMPUT file and creates a
standard event file using the instrument simulator. :func:`~sox.events.make_event_file`
performs the following actions:

1. Uses the effective area curve to determine which events will actually be detected.
2. Projects these events onto the detector plane and perform dithering of their positions.
3. Convolves the event energies with the response matrix to produce channels.

A typical invocation of :func:`~sox.events.make_event_file` looks like the following:

.. code-block:: python

    from sox import make_event_file
    simput_file = "snr_simput.fits" # SIMPUT file to be read
    out_file = "evt_xcal.fits" # event file to be written
    exp_time = 30000. # The exposure time in seconds
    instrument = "xcal" # short name for instrument to be used
    sky_center = [30., 45.] # RA, Dec of pointing in degrees
    make_event_file(simput_file, out_file, exp_time, instrument, sky_center, clobber=True)
    
If you have your own JSON-based instrument specification, you can supply it to the instrument
argument instead:

.. code-block:: python

    instrument = "my_imager.json"
    make_event_file(simput_file, out_file, instrument, sky_center, clobber=True)

Available instruments are:

* ``"hdxi"``: 
* ``"xcal"``:

Adding Background Events
------------------------


