.. _spectra:

Creating and Using Spectra
==========================

XRStools provides a way to create common types of spectra that can then be
used in your scripts to create mock observations via the 
:class:`~xrs_tools.spectra.Spectrum` object.

Creating a Spectrum Object
--------------------------

A :class:`~xrs_tools.spectra.Spectrum` object is simply defined by a table 
of energies and photon fluxes. There are several ways to create a 
:class:`~xrs_tools.spectra.Spectrum`, depending on your use case. 

Reading a Spectrum from a File
++++++++++++++++++++++++++++++

If you have a spectrum tabulated in an ASCII text file, this can be read
in using the :meth:`~xrs_tools.spectra.Spectrum.from_file` method. The file
can be comprised of either two or three columns; the former 

Creating a Power-Law Spectrum
+++++++++++++++++++++++++++++


Generating Thermal Spectra
++++++++++++++++++++++++++

Applying Galactic Foreground Absorption to a Spectrum
-----------------------------------------------------

Generating Photons From a Spectrum
----------------------------------

