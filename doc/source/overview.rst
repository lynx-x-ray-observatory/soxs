.. _overview:

SOXS Overview
=============

This section provides a general overview of the capabilities of SOXS. SOXS has 
three main components:

1. Tools to generate simple models of astrophysical sources
2. I/O interface to SIMPUT files
3. An instrument simulator for X-ray telescopes

These are outlined in detail below. Figure 1 shows a flowchart of how one 
transforms models of astrophysical sources to event lists using SOXS.

.. figure:: images/flowchart.png
    :align: center
    :figclass: w
    :scale: 20 %

    Figure 1: Flowchart depicting how source models are written to SIMPUT files
    which then can be used to simulate X-ray observations using SOXS and compare 
    them with the outputs of other tools.

Generating Simplified Models of Sources
---------------------------------------

SOXS provides tools to build simplified models of sources on the sky. These consist of
``Spectrum`` and ``SpatialModel`` objects. The first creates models of spectra of a number
of different types and a number of different sources. Methods are also provided to

To find out more about ``Spectrum`` objects, visit one of the following links:

* Python interface: :ref:`spectra`
* Command-line interface: :ref:`cmd-spectra`

To learn about ``SpatialModel`` objects, visit one of the following links:

* Python interface: :ref:`spatial`
* Command-line interface: :ref:`cmd-spatial`

The purpose of these methods is to generate statistically representative samples of photons
which can be written to SIMPUT files, described next.

SIMPUT Files
------------

For storage and representation of source models, SOXS uses the SIMPUT file 
format. SIMPUT stands for "SIMulated inPUT", and was developed at the 
University of Erlangen-Nuremberg for use with the 
`SIXTE <http://www.sternwarte.uni-erlangen.de/research/sixte/index.php>`_
mock X-ray observation package. However, other similar packages, such as 
`SIMX <http://hea-www.cfa.harvard.edu/simx/>`_ and 
`MARX <http://space.mit.edu/CXC/MARX/>`_ use SIMPUT as a possible input format 
for source models, making it possible to used models produced in or for SOXS 
with these packages to compare the expected performance of *Lynx* with other 
instruments. 

SOXS currently provides facilities for both input and output of SIMPUT catalogs 
with spectra, images, and photon lists. A SIMPUT spectral model is simply a FITS 
table of energy and flux from a source. These may correspond to a point source, 
or if included with a FITS image it may be an extended source. A SIMPUT photon 
list model is a FITS table of RA, Dec, and energy of photons from the simulated 
source, generated assuming a large exposure time and a large collecting area so 
that the sample is large enough that the instrument simulator is able to 
"observe" a representative subset. The SIMPUT outputs can be used to compare 
simulations made with SOXS to simulations of other X-ray instruments past, 
current, and future, and likewise SIMPUT catalogs with photon lists not created 
with SOXS can be used as inputs to SOXS's instrument simulator.

For more information, visit :ref:`simput`.

The Instrument Simulator
------------------------

.. |instrument_simulator_cmd| replace:: ``instrument_simulator`` command-line script
.. _instrument_simulator_cmd: command_line/instrument.html#simulate-events

.. |instrument_simulator_py| replace:: ``instrument_simulator()`` Python function
.. _instrument_simulator_py: users_guide/instrument.html#running-the-instrument-simulator

The instrument simulator in SOXS takes unconvolved events in the form of a
SIMPUT file and performs the following operations:
 
1. Uses the effective area curve from an ARF to determine which events will 
   actually be detected.
2. Projects these events onto the detector plane and perform PSF blurring and 
   dithering of their positions.
3. Add particle/instrumental and astrophysical background events.
4. Convolves the event energies with the response matrix from an RMF to produce
   channels.
5. Writes everything to an event file.

The instrument simulator is called using either the |instrument_simulator_py|_ 
or the |instrument_simulator_cmd|_. 

Currently, the instrument simulator can simulate certain "default" instrument 
configurations for *Lynx*, *Athena*, *Chandra*, *AXIS*, *XRISM*, and *STAR-X*, 
but one can also supply a modified instrument configuration for use with the 
instrument simulator, which is laid out in more detail in :ref:`instrument`. 

Working with Other Tools
------------------------

SOXS is designed so that the SIMPUT files it produces may be used in other
X-ray simulation tools, and other tools which produce SIMPUT files may then be 
read in by SOXS. Here we list a few of the tools and software that you may have 
interest in using:

MARX
++++

Website: http://space.mit.edu/CXC/MARX/

MARX is a set of programs developed to provide a detailed ray-tracing simulation 
of the on-orbit performance of *Chandra*. The SIMPUT files produced by SOXS can 
be used as inputs to MARX to simulate *Chandra* observations to compare with 
those made by SOXS.

SIMX
++++

Website: http://hea-www.cfa.harvard.edu/simx/

SIMX simulates a photon-counting detector's response to an input source, 
including a simplified model of any telescope. SIMX is a "convolution tool" that
uses standard descriptions of telescope PSF and the detector response to model 
how sources will appear. The SIMPUT files produced by SOXS can be used as inputs
to SIMX, and may be useful for simulating observations using other instruments, 
such as *Athena*, *XRISM*, etc.

SIXTE
+++++

Website: http://www.sternwarte.uni-erlangen.de/research/sixte/index.php

SIXTE is a software package for X-ray telescope observation simulations 
developed at the Erlangen Centre for Astroparticle Physics (ECAP). It allows 
one to undertake instrument performance analyses and to produce simulated 
event files for mission and analysis studies. Its primary goal is to produce 
simulated *Athena* observations, but it can produce observations of several 
other missions as well. The SIMPUT files produced by SOXS can be used as 
inputs to SIXTE.

pyXSIM
++++++

Website: http://hea-www.cfa.harvard.edu/~jzuhone/pyxsim/

pyXSIM is a Python package for simulating X-ray observations from 3D models of
astrophysical sources. pyXSIM makes it possible to generate synthetic X-ray 
observations of these sources from a wide variety of models, whether from 
grid-based simulation codes such as FLASH, Enzo, and Athena, to particle-based
codes such as Gadget and AREPO, and even from datasets that have been created 
"by hand", such as from NumPy arrays. pyXSIM can be used to produce SIMPUT files
which can be ingested by SOXS for making simulated observations.