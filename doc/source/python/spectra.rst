.. _spectra:

Creating and Using Spectra
==========================

SOXS provides a way to create common types of spectra that can then be
used in your scripts to create mock observations via the 
:class:`~soxs.spectra.Spectrum` object.

Creating Spectrum Objects
-------------------------

A :class:`~soxs.spectra.Spectrum` object is simply defined by a table 
of energies and photon fluxes. There are several ways to create a 
:class:`~soxs.spectra.Spectrum`, depending on your use case. 

Reading a Spectrum from a File
++++++++++++++++++++++++++++++

If you have a spectrum tabulated in an ASCII text file, this can be read
in using the :meth:`~soxs.spectra.Spectrum.from_file` method. The file
must be comprised of two columns, the first being the energies of the bins
in keV and the second being the photon flux in units of 
:math:`{\rm photons}~{\rm cm}^{-2}~{\rm s}^{-1}~{\rm keV}^{-1}`. The binning 
must be linear and the bins must be equally spaced. For example:

.. code-block:: python

    from soxs import Spectrum
    my_spec = Spectrum.from_file("my_spec.dat")

Creating a Power-Law Spectrum
+++++++++++++++++++++++++++++

A simple power-law spectrum can be created using the 
:meth:`~soxs.spectra.Spectrum.from_powerlaw` method. This takes as input
a spectral index ``photon_index``, a redshift ``redshift``, and a normalization
of the source ``norm`` at 1 keV in the source frame, in units of 
:math:`{\rm photons}~{\rm cm}^{-2}~{\rm s}^{-1}~{\rm keV}^{-1}`. Mathematically, 
this is equivalent to:

.. math::

    F_E = K\left[\frac{E(1+z)}{{\rm 1~keV}}\right]^{-\alpha}
    
where :math:`\alpha` is the ``photon_index`` (note the sign convention). You can set
up a power-law spectrum like this:

.. code-block:: python

    alpha = 1.2
    zobs = 0.05
    norm = 1.0e-7
    spec = Spectrum.from_powerlaw(alpha, zobs, norm, emin=0.1, emax=10.0, nbins=20000)

The optional parameters ``emin``, ``emax``, and ``nbins`` can be used to control the
binning. 

Generating Thermal Spectra
++++++++++++++++++++++++++

Thermal spectra are generated in SOXS using the `AtomDB tables <http://www.atomdb.org>`_, 
and require special handling. The :class:`~soxs.spectra.ApecGenerator` class is a factory
class which generates new :class:`~soxs.spectra.Spectrum` objects. You start by initializing
an :class:`~soxs.spectra.ApecGenerator`:

.. code-block:: python

    from soxs import ApecGenerator
    agen = ApecGenerator(0.05, 50.0, 10000, apec_vers="2.0.2", broadening=True)

The ``broadening`` parameter sets whether or not spectral lines will be thermally and
velocity broadened. The ``apec_vers`` parameter sets the version of the AtomDB tables
to use. Versions 2.0.2 and 3.0.3 are built into SOXS. 

You may also supply another location for the AtomDB tables. For example, the following 
construction will look for the AtomDB tables in the current working directory:

.. code-block:: python

    agen = ApecGenerator(0.05, 50.0, 10000, apec_root=".")

Once you have an :class:`~soxs.spectra.ApecGenerator` object, you can use it to generate
thermal spectra: 

.. code-block:: python
    
    kT = 6.0 # in units of keV
    abund = 0.3 # solar units
    redshift = 0.05
    norm = 1.0e-3 # in units of 1.0e-14*EM/(4*pi*(1+z)**2*D_A**2)
    velocity = 100.0 # in units of km/s, optional
    spec1 = agen.get_spectrum(kT, abund, redshift, norm, velocity=velocity)

``spec1`` is just a standard :class:`~soxs.spectra.Spectrum` object.

Generating a Spectrum from XSPEC
++++++++++++++++++++++++++++++++

If you have XSPEC installed on your machine, you can use it with SOXS to create any 
spectral model that XSPEC supports. This is done by passing in a model string and a
list of parameters:

.. code-block:: python

    model_string = "phabs*(mekal+powerlaw)" # A somewhat complicated model
    params = [0.02, 6.0, 1.0, 0.3, 0.03, 1, 0.01, 1.2, 1.0e-3]
    spec = Spectrum.from_xspec(model_string, params, emin=0.1, emax=1.0, nbins=20000)
    
Note that the parameters must be in the same order that they would be if you were entering
them in XSPEC. The ``emin``, ``emax``, and ``nbins`` keyword arguments are used to control
the energy binning.

.. note::

    Generating spectra from XSPEC requires that the ``HEADAS`` environment is sourced
    before running the Python script, as it would be if you were using XSPEC to fit 
    spectra. 

Adding Spectra Together
-----------------------

Two :class:`~soxs.spectra.Spectrum` objects can be co-added, provided that
they have the same energy binning:

.. code-block:: python
 
    spec1 = Spectrum.from_powerlaw(1.1, 0.05, 1.0e-9)
    spec2 = agen.get_spectrum(6.0, 0.3, 0.05, 1.0e-3)

    total_spectrum = spec1 + spec2
    
If they do not, an error will be thrown. 

Rescaling the Normalization of a Spectrum
-----------------------------------------

You can rescale the normalization of the entire spectrum using the
:meth:`~soxs.spectra.Spectrum.rescale_flux` method. This can be 
helpful when you want to set the normalization of the spectrum by the 
total flux within a certain energy band instead. 

.. code-block:: python

    spec.rescale_flux(1.0e-9, emin=0.5, emax=7.0, flux_type="photons"):

``emin`` and ``emax`` can be used to set the band that the flux corresponds to. If they
are not set, they are assumed to be the bounds of the spectrum. The flux type can be 
``"photons"`` (the default) or ``"energy"``. In the former case, the units of the new 
flux must be :math:`{\rm photons}~{\rm cm}^{-2}~{\rm s}^{-1}`, and in the latter case 
the units must be :math:`{\rm erg}~{\rm cm}^{-2}~{\rm s}^{-1}`.

Applying Galactic Foreground Absorption to a Spectrum
-----------------------------------------------------

The :meth:`~soxs.spectra.Spectrum.apply_foreground_absorption` method
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

Given a :class:`~soxs.spectra.Spectrum`, a set of photon energies can be 
drawn from it using the :meth:`~soxs.spectra.Spectrum.generate_energies`
method. This will most often be used to generate discrete samples for mock 
observations. For this method, an exposure time and a constant effective area 
must be supplied to convert the spectrum's flux to a number of photons. These
values need not be realistic--in fact, they both should be larger than the 
values for the mock observation that you want to simulate, to create a statistically
robust sample to draw photons from. 

An example using a :class:`~soxs.spectra.Spectrum` created from a file:

.. code-block:: python

    spec = Spectrum.from_file("my_spec.dat")
    t_exp = 100000. # exposure time in seconds
    area = 30000. # constant effective area
    energies = spec.generate_energies(t_exp, area)

These photon energies can then be combined with sky positions at your discretion
and be written to SIMPUT files for use in mock observations. See :ref:`simput` for
more information.