.. _instrument:

The Instrument Simulator
========================

Running the Instrument Simulator
--------------------------------

The end product of a mock observation is a "standard" event file which has been 
convolved with a model for the telescope. In SOXS, this is handled by the
instrument simulator. 

:func:`~soxs.instrument.instrument_simulator` reads in a SIMPUT catalog and creates a
standard event file using the instrument simulator. :func:`~soxs.instrument.instrument_simulator`
performs the following actions:

1. Uses the effective area curve to determine which events will actually be detected.
2. Projects these events onto the detector plane and perform PSF blurring and dithering 
   of their positions.
3. Add particle/instrumental background events. 
4. Convolves the event energies with the response matrix to produce channels.
5. Writes everything to an event file.

All of the photon lists in the SIMPUT catalog will be processed. A typical invocation of 
:func:`~soxs.instrument.instrument_simulator` looks like the following:

.. code-block:: python

    from soxs import instrument_simulator
    simput_file = "snr_simput.fits" # SIMPUT file to be read
    out_file = "evt_xcal.fits" # event file to be written
    exp_time = 30000. # The exposure time in seconds
    instrument = "xcal" # short name for instrument to be used
    sky_center = [30., 45.] # RA, Dec of pointing in degrees
    instrument_simulator(simput_file, out_file, exp_time, instrument, sky_center, clobber=True)
 
The ``clobber`` argument allows an existing file to be overwritten. Instruments must exist
in the instrument registry, unless you have your own JSON-based instrument specification, 
which you can then supply as the instrument argument instead:

.. code-block:: python

    instrument = "my_imager.json"
    instrument_simulator(simput_file, out_file, instrument, sky_center, clobber=True)

Available instruments are:

* ``"hdxi"``: 
* ``"xcal"``:

You can also change other aspects of the observation with :func:`~soxs.instrument.instrument_simulator`. 
For example, you can change the shape and size of the dither pattern:


You can also specify a non-zero roll angle:

The particle background scale can be set using the ``bkgnd_scale`` argument:

Customizing the Instrument Simulator
------------------------------------

The Instrument Registry
+++++++++++++++++++++++

Making Custom Instruments with JSON Files
+++++++++++++++++++++++++++++++++++++++++

