.. _overview:

SOXS Overview
=============

This page provides a general overview of the capabilities of SOXS. 

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

Working with Other Tools
------------------------

SOXS is designed so that the SIMPUT files it produces may be used in other
X-ray simulation tools, and other tools which produce SIMPUT files may then be read
in by SOXS. Here we list a few of the tools and software that you may have interest
in using:

MARX
++++

Website: 

SIMX
++++

Website: http://hea-www.cfa.harvard.edu/simx

simx simulates a photon-counting detector's response to an input source, including 
a simplified model of any telescope. The code is not a full ray-trace, but a convolution 
tool that uses standard descriptions of telescope PSF (via either a simple Gaussian 
parameter, an energy-dependent encircled-energy function, or an image of the PSF) and
the detector response (using the OGIP response function) to model how sources will appear.

The SIMPUT files produced by SOXS can be used as inputs to SIMX, and may be useful for 
simulating observations using other instruments, such as Athena, Hitomi, etc.

pyXSIM
++++++

Website: http://hea-www.cfa.harvard.edu/~jzuhone/pyxsim

pyXSIM is a Python package for simulating X-ray observations from 3-D models of
astrophysical sources. pyXSIM makes it possible to generate synthetic X-ray 
observations of these sources from a wide variety of models, whether from grid-based 
simulation codes such as FLASH, Enzo, and Athena, to particle-based codes such as 
Gadget and AREPO, and even from datasets that have been created “by hand”, such as from
NumPy arrays. pyXSIM also provides facilities for manipulating the synthetic observations 
it produces in various ways, as well as ways to export the simulated X-ray events to other
software packages to simulate the end products of specific X-ray observatories.

pyXSIM can be used to produce SIMPUT files which can be ingested by SOXS for making
simulated X-ray Surveyor observations.
