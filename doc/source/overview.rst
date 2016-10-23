.. _overview:

SOXS Overview
=============

Generating Spectra and Events
-----------------------------

SIMPUT Files
------------

The Instrument Simulator
------------------------

The instrument simulator in SOXS takes unconvolved events in the form of a
SIMPUT file and performs the following operations:
 
1. Uses the effective area curve from an ARF to determine which events will 
   actually be detected.
2. Projects these events onto the detector plane and perform PSF blurring and 
   dithering of their positions.
3. Add particle/instrumental background events. 
4. Convolves the event energies with the response matrix from an RMF to produce
   channels.
5. Writes everything to an event file.
