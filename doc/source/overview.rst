.. _overview:

SOXS Overview
=============

This section provides a general overview of the capabilities of SOXS. SOXS has 
three main goals

Generating Simplified Models of Sources
---------------------------------------

SIMPUT Files
------------

For storage and representation of source models, SOXS uses the SIMPUT file format. SIMPUT
stands for "SIMulated inPUT", and was developed at the University of Erlangen-Nuremberg
for use with the `SIXTE <http://www.sternwarte.uni-erlangen.de/research/sixte/index.php>`_
mock X-ray observation package. However, other similar packages, such as 
`SIMX <http://hea-www.cfa.harvard.edu/simx/>`_ and `MARX <http://space.mit.edu/CXC/MARX/>`_
use SIMPUT as a possible input format for source models, making it possible to used models
produced in or for SOXS with these packages to compare the expected performance of X-ray
Surveyor with other instruments. 

SOXS currently provides facilities for both input and output of SIMPUT catalogs with
photon lists. The SIMPUT outputs can be used to compare simulations made with SOXS for 
X-ray Surveyor to simulations of other X-ray instruments past, current, and future, and 
likewise SIMPUT catalogs with photon lists not created with SOXS can be used as inputs 
to SOXS's instrument simulator. 

The Instrument Simulator
------------------------

.. |simulate_events_cmd| replace:: ``simulate_events`` command-line script
.. _simulate_events_cmd: command_line/instrument.html#simulate-events

.. |simulate_events_py| replace:: ``simulate_events()`` Python function
.. _simulate_events_py: python/instrument.html#creating-event-files

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

The instrument simulator is called using either the |simulate_events_py|_ or the
|simulate_events_cmd|_. 

The instrument simulator can simulate "default" instrument configurations for 
X-ray Surveyor, but one can also supply a modified instrument configuration for use
with the instrument simulator, which is laid out in more detail in :ref:`instrument`. 

Working with Other Tools
------------------------

SOXS is designed so that the SIMPUT files it produces may be used in other
X-ray simulation tools, and other tools which produce SIMPUT files may then be read
in by SOXS. Here we list a few of the tools and software that you may have interest
in using:

MARX
++++

Website: http://space.mit.edu/CXC/MARX/

MARX is a set of programs developed to provide a detailed ray-tracing simulation of the
on-orbit performance of *Chandra*. The SIMPUT files produced by SOXS can be used as inputs
to MARX to simulate *Chandra* observations to compare with X-ray Surveyor. 

SIMX
++++

Website: http://hea-www.cfa.harvard.edu/simx/

SIMX simulates a photon-counting detector's response to an input source, including 
a simplified model of any telescope. SIMX is a "convolution tool" that uses standard 
descriptions of telescope PSF and the detector response to model how sources will 
appear. The SIMPUT files produced by SOXS can be used as inputs to SIMX, and may be 
useful for simulating observations using other instruments, such as *Athena*, *Hitomi*, 
etc.

SIXTE
+++++

Website: http://www.sternwarte.uni-erlangen.de/research/sixte/index.php

SIXTE is a software package for X-ray telescope observation simulations developed 
at the Erlangen Centre for Astroparticle Physics (ECAP) under the leadership of Christian
Schmid. It allows to undertake instrument performance analyses and to produce simulated 
event files for mission and analysis studies. Its primary goal is to produce simulated
*Athena* observations, but it can produce observations of several other missions as
well. The SIMPUT files produced by SOXS can be used as inputs to SIXTE.

pyXSIM
++++++

Website: http://hea-www.cfa.harvard.edu/~jzuhone/pyxsim/

pyXSIM is a Python package for simulating X-ray observations from 3-D models of
astrophysical sources. pyXSIM makes it possible to generate synthetic X-ray 
observations of these sources from a wide variety of models, whether from grid-based 
simulation codes such as FLASH, Enzo, and Athena, to particle-based codes such as 
Gadget and AREPO, and even from datasets that have been created “by hand”, such as from
NumPy arrays. pyXSIM can be used to produce SIMPUT files which can be ingested by 
SOXS for making simulated X-ray Surveyor observations.