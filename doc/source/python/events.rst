.. _events:

Creating Event Files
====================

The end product of a mock observation is a "standard" event file which has been 
convolved with a model for the telescope. In XRStools, this is handled by the 
instrument simulator. 

:func:`~xrs_tools.events.make_event_file` reads in a SIMPUT file and creates a 
standard event file using the instrument simulator. :func:`~xrs_tools.events.make_event_file`
performs the following actions:

1. Uses the effective area curve to determine which events will actually be detected.
2. Projects these events onto the detector plane and perform dithering of their positions.
3. Convolves the event energies with the response matrix to produce channels.

Typical invocations of :func:`~xrs_tools.events.make_event_file` look like this:

.. code-block:: python

    make_event_file()