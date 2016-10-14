.. _spectra:

Creating and Using Spectra
==========================

XRStools provides a way to create common types of spectra that can then be
used in your scripts to create mock observations via the 
:class:`~xrs_tools.spectra.Spectrum` object.

Creating Spectrum Objects
-------------------------

A :class:`~xrs_tools.spectra.Spectrum` object is simply defined by a table 
of energies and photon fluxes. There are several ways to create a 
:class:`~xrs_tools.spectra.Spectrum`, depending on your use case. 

Reading a Spectrum from a File
++++++++++++++++++++++++++++++

If you have a spectrum tabulated in an ASCII text file, this can be read
in using the :meth:`~xrs_tools.spectra.Spectrum.from_file` method. The file
must be comprised of two columns, the first being the energies of the bins
in keV and the second being the photon flux in units of 
:math:`{\rm photons}~{\rm cm}^{-2}~{\rm s}^{-1}`. The binning must be linear
and the bins must be equally spaced. For example:

.. code-block:: python

    my_spec = Spectrum.from_file("my_spec.dat")

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

    apec_model = ApecGenerator(0.05, 50.0, 10000)
    
    spec1 = Spectrum.from_powerlaw(1.1, 0.05, 1.0e-9)
    spec2 = apec_model.get_spectrum(6.0, 0.3, 0.05, 1.0e-3)

    total_spectrum = spec1 + spec2
    
If they do not, an error will be thrown. 

Rescaling the Normalization of a Spectrum
-----------------------------------------

You can rescale the normalization of the entire spectrum using the
:meth:`~xrs_tools.spectra.Spectrum.rescale_flux` method. This can be 
helpful when you want to set the normalization of the spectrum by the 
total flux within a certain energy band instead. 

.. code-block:: python

    spec = Spectrum.from_xspec()
    spec.rescale_flux(1.0e-9, emin=0.5, emax=7.0, flux_type="photons"):

``emin`` and ``emax`` can be used to set the band that the flux corresponds to. If they
are not set, they are assumed to be the bounds of the spectrum. The flux type can be 
``"photons"`` (the default) or ``"energy"``. In the former case, the units of the new 
flux must be :math:`{\rm photons}~{\rm cm}^{-2}~{\rm s}^{-1}`, and in the latter case 
the units must be :math:`{\rm erg}~{\rm cm}^{-2}~{\rm s}^{-1}`.

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