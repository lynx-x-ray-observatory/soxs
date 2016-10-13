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

A simple power-law spectrum can be created using the 
:meth:`~xrs_tools.spectra.Spectrum.from_powerlaw` method. This takes as input 

Generating Thermal Spectra
++++++++++++++++++++++++++

Generating a Spectrum from XSPEC
++++++++++++++++++++++++++++++++

Adding Spectra Together
-----------------------

Two :class:`~xrs_tools.spectra.Spectrum` objects can be co-added, provided that
they have the same energy binning:

.. code-block:: python

    total_spectrum = plaw_spectrum + thermal_spectrum
    
If they do not, an error will be thrown. 

Rescaling the Normalization of a Spectrum
-----------------------------------------


Applying Galactic Foreground Absorption to a Spectrum
-----------------------------------------------------

Generating Photons From a Spectrum
----------------------------------

