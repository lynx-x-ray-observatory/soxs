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

Creating a Constant Spectrum
++++++++++++++++++++++++++++

A simple constant spectrum can be created using the :meth:`~soxs.spectra.Spectrum.from_constant` 
method. This takes as input the value of the flux ``const_flux``, which is in units of 
:math:`{\rm photons}~{\rm cm}^{-2}~{\rm s}^{-1}~{\rm keV}^{-1}`:

.. code-block:: python

    const_flux = 1.0e-7
    spec = Spectrum.from_constant(const_flux, emin=0.1, emax=10.0, nbins=20000)

The optional parameters ``emin``, ``emax``, and ``nbins`` can be used to control the
binning. 

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
to use. Versions 2.0.2 and 3.0.3 are built into SOXS, the default version is 3.0.3.

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
spectral model that XSPEC supports. You can do this in two ways. The first is by 
passing in a model string and a list of parameters to the :meth:`~soxs.spectra.Spectrum.from_xspec_model`
method:

.. code-block:: python

    model_string = "phabs*(mekal+powerlaw)" # A somewhat complicated model
    params = [0.02, 6.0, 1.0, 0.3, 0.03, 1, 0.01, 1.2, 1.0e-3]
    spec = Spectrum.from_xspec_model(model_string, params, emin=0.1, emax=1.0, nbins=20000)
    
Note that the parameters must be in the same order that they would be if you were entering
them in XSPEC. The ``emin``, ``emax``, and ``nbins`` keyword arguments are used to control
the energy binning.

The second way involves passing an XSPEC script file to the :meth:`~soxs.spectra.Spectrum.from_xspec_script`
method which defines an XSPEC model. For example, a script that creates a model spectrum
from a sum of two APEC models may look like this:

.. code-block:: text

    statistic chi
    method leven 10 0.01
    abund angr
    xsect bcmc
    cosmo 70 0 0.73
    xset delta 0.01
    systematic 0
    model  apec    +   apec
                0.2       0.01      0.008      0.008         64         64
                  1     -0.001          0          0          5          5
                  0      -0.01     -0.999     -0.999         10         10
        6.82251e-07       0.01          0          0      1e+24      1e+24
              0.099       0.01      0.008      0.008         64         64
                  1     -0.001          0          0          5          5
                  0      -0.01     -0.999     -0.999         10         10
        1.12328e-06       0.01          0          0      1e+24      1e+24

If it is contained within the file ``"two_apec.xcm"``, it can be used to create a :class:`~soxs.spectra.Spectrum`
like this:

.. code-block:: python

    spec = Spectrum.from_xspec_script("two_apec.xcm", emin=0.1, emax=1.0, nbins=20000)

.. note::

    Generating spectra from XSPEC requires that the ``HEADAS`` environment is sourced
    before running the Python script, as it would be if you were using XSPEC to fit 
    spectra. 

Math with ``Spectrum`` Objects
------------------------------

Two :class:`~soxs.spectra.Spectrum` objects can be co-added, provided that
they have the same energy binning:

.. code-block:: python
 
    spec1 = Spectrum.from_powerlaw(1.1, 0.05, 1.0e-9)
    spec2 = agen.get_spectrum(6.0, 0.3, 0.05, 1.0e-3)

    total_spectrum = spec1 + spec2
    
If they do not, an error will be thrown. 

You can also multiply a spectrum by a constant float number or divide it by one:

.. code-block:: python

    spec3 = 6.0*spec2
    spec4 = spec1/4.4

Getting the Flux of a Spectrum Within a Specific Energy Band
------------------------------------------------------------

The :meth:`~soxs.spectra.Spectrum.get_flux_in_band` method can be used
to quickly report on the total flux within a specific energy band:

.. code-block:: python

    emin = 0.5
    emax = 7.0
    print(spec.get_flux_in_band(emin, emax))

which returns a tuple of the photon flux and the energy flux, showing:

.. code-block:: pycon

    (<Quantity 2.2215588675210208e-07 ph / (cm2 s)>, <Quantity 7.8742710307246895e-16 erg / (cm2 s)>)
    
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
observations. For this method, an exposure time and a constant (energy-independent)
effective area must be supplied to convert the spectrum's flux to a number of photons. These
values need not be realistic--in fact, they both should be larger than the 
values for the mock observation that you want to simulate, to create a statistically
robust sample to draw photons from when we actually pass them to the instrument simulator.

An example using a :class:`~soxs.spectra.Spectrum` created from a file:

.. code-block:: python

    spec = Spectrum.from_file("my_spec.dat")
    t_exp = 100000. # exposure time in seconds
    area = 30000. # constant effective area
    energies = spec.generate_energies(t_exp, area)

The ``energies`` object :meth:`~soxs.spectra.Spectrum.generate_energies` returns is an augmented
NumPy array which also carries the unit information and the total flux of energies:

.. code-block:: python

    print(energies.unit)
    print(energies.flux)

.. code-block:: pycon

    Unit("keV")
    <Quantity 1.1256362913845828e-15 erg / (cm2 s)>

These photon energies can then be combined with sky positions at your discretion
and be written to SIMPUT files for use in mock observations. See :ref:`simput` for
more information.

"Convolved" Spectra
-------------------

One may want to examine a spectrum after it has been convolved with a particular effective
area curve. One can generate such a spectrum using :class:`~soxs.spectra.ConvolvedSpectrum`
from a :class:`~soxs.spectra.Spectrum` object and an ARF:

.. code-block:: python

    from soxs import ConvolvedSpectrum
    # Assuming one created an ApecGenerator agen...
    spec2 = agen.get_spectrum(6.0, 0.3, 0.05, 1.0e-3)
    cspec = ConvolvedSpectrum(spec2, "xrs_hdxi_3x10.arf")
    
The spectrum in this object has units of :math:`{\rm photons}~{\rm s}^{-1}~{\rm keV}^{-1}`,
and one can use all of :class:`~soxs.spectra.Spectrum`'s methods on it. For example, to 
determine the count and energy rate within a particular band:

.. code-block:: python

    cspec.get_flux_in_band(0.5, 7.0)

.. code-block:: python

    (<Quantity 6.802363401824924 ph / s>,
     <Quantity 1.2428592072628134e-08 erg / s>)

Or to generate an array of energies:

.. code-block:: python

    t_exp = 500000. # in seconds
    e = cspec.generate_energies(t_exp)

:class:`~soxs.spectra.ConvolvedSpectrum` objects are not used directly in the 
instrument simulator, but can be used for convenient when one wants to examine the
properties of a convolved spectrum.