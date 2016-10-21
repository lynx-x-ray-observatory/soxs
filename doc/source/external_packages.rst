.. _external-packages:

External Packages
=================

The following are external packages that may be used with SOXS or in place of some of its
features.

SIMX
----

Website: http://hea-www.cfa.harvard.edu/simx

simx simulates a photon-counting detector's response to an input source, including 
a simplified model of any telescope. The code is not a full ray-trace, but a convolution 
tool that uses standard descriptions of telescope PSF (via either a simple Gaussian 
parameter, an energy-dependent encircled-energy function, or an image of the PSF) and
the detector response (using the OGIP response function) to model how sources will appear.

simx is NOT a full ray-trace code, but rather uses a predefined set of PSFs, vignetting 
information, and instrumental responses and outputs to make the simulation. It is designed
to be a "approximation" tool to estimate issues such as source confusion, background effects,
pileup, and other similar issues.

The SIMPUT files produced by SOXS can be used as inputs to SIMX.

pyXSIM
------

Website: http://hea-www.cfa.harvard.edu/~jzuhone/pyxsim

pyXSIM is a Python package for simulating X-ray observations from 3-D models of
astrophysical sources. pyXSIM makes it possible to generate synthetic X-ray 
observations of these sources from a wide variety of models, whether from grid-based 
simulation codes such as FLASH, Enzo, and Athena, to particle-based codes such as 
Gadget and AREPO, and even from datasets that have been created “by hand”, such as from
NumPy arrays. pyXSIM also provides facilities for manipulating the synthetic observations 
it produces in various ways, as well as ways to export the simulated X-ray events to other
software packages to simulate the end products of specific X-ray observatories.

pyXSIM can be used to produce SIMPUT files which can be ingested by SOXS.
