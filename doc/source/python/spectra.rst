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

The :meth:`~xrs_tools.spectra.Spectrum.apply_foreground_absorption` method
can be used to apply foreground absorption using the "wabs" model. It takes 
one parameter, the hydrogen column along the line of sight, in units of 
:math:`10^{22}~{\rm cm}^{-2}`:

.. code-block:: python

    spec = Spectrum.from_powerlaw(1.1, 0.05, 1.0e-9)
    n_H = 0.02
    spec.apply_foreground_absorption(n_H)

The flux in the energy bins will be reduced according to the absorption at a
given energy.

Generating Photon Energies From a Spectrum
------------------------------------------

Given a :class:`~xrs_tools.spectra.Spectrum`, a set of photon energies can be 
drawn from it using the :meth:`~xrs_tools.spectra.Spectrum.generate_energies`
method. This will most often be used to generate discrete samples for mock 
observations. For this method, an exposure time and a constant effective area 
must be supplied to convert the spectrum's flux to a number of photons. These
values need not be realistic--in fact, they both should be larger than the 
values for the mock observation that you want to simulate, to create a statistically
robust sample to draw photons from. 

An example using a :class:`~xrs_tools.spectra.Spectrum` created from a file:

.. code-block:: python

    spec = Spectrum.from_file("my_spec.dat")
    t_exp = 100000. # exposure time in seconds
    area = 30000. # constant effective area
    energies = spec.generate_energies(t_exp, area)

These photon energies can then be combined with sky positions at your discretion
and be written to SIMPUT files for use in mock observations. See :ref:`simput` for
more information.