.. _spectra:

Creating and Using Spectra
==========================

SOXS provides a way to create common types of spectra that can then be
used in your scripts to create mock observations via the
:class:`~soxs.spectra.Spectrum` object.

.. _spectrum-binning:

Spectrum Binning
----------------

The energy binning of spectral tables can be either linear or log--that is,
either the difference between the minimum and maximum energies of each bin is
constant across the spectrum (linear) or that the difference between the logarithm
of the minimum and maximum energies of each bin is constant across the spectrum
(log).

For most of the spectrum creation methods outlined below, there will be the following
keyword arguments to control the binning of spectral tables:

* ``emin``: The minimum energy of the spectral table in keV.
* ``emax``: The maximum energy of the spectral table in keV.
* ``nbins``: The number of bins in the spectrum.
* ``binscale``: An optional argument which takes either ``"linear"`` or ``"log"``.
  The default is always ``"linear"``.

Creating Spectrum Objects
-------------------------

A :class:`~soxs.spectra.Spectrum` object is simply defined by a table
of energies and photon fluxes. There are several ways to create a
:class:`~soxs.spectra.Spectrum`, depending on your use case.

Creating a Constant Spectrum
++++++++++++++++++++++++++++

A simple constant spectrum can be created using the
:meth:`~soxs.spectra.Spectrum.from_constant` method. This takes as input the
value of the flux ``const_flux``, which is in units of
:math:`{\rm photons}~{\rm cm}^{-2}~{\rm s}^{-1}~{\rm keV}^{-1}`. The parameters
``emin``, ``emax``, ``nbins``, and ``binscale`` are used to control the binning.

.. code-block:: python

    const_flux = 1.0e-7
    emin = 0.1
    emax = 10.0
    nbins = 20000
    spec = Spectrum.from_constant(const_flux, emin, emax, nbins, binscale="linear")

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

where :math:`\alpha` is the ``photon_index`` (note the sign convention). The parameters
``emin``, ``emax``, ``nbins``, and ``binscale`` are used to control the binning.

You can set up a power-law spectrum like this:

.. code-block:: python

    alpha = 1.2
    zobs = 0.05
    norm = 1.0e-7
    emin = 0.1
    emax = 10.0
    nbins = 20000
    spec = Spectrum.from_powerlaw(alpha, zobs, norm, emin, emax, nbins, binscale="log")

.. _xspec:

Generating a Spectrum from XSPEC
++++++++++++++++++++++++++++++++

If you have XSPEC installed on your machine, you can use it with SOXS to create
any spectral model that XSPEC supports. You can do this in two ways. The first
is by passing in a model string and a list of parameters to the
:meth:`~soxs.spectra.Spectrum.from_xspec_model` method:

.. code-block:: python

    model_string = "phabs*(mekal+powerlaw)" # A somewhat complicated model
    params = [0.02, 6.0, 1.0, 0.3, 0.03, 1, 0.01, 1.2, 1.0e-3]
    emin = 0.1
    emax = 5.0
    nbins = 20000
    spec = Spectrum.from_xspec_model(model_string, params, emin, emax, nbins)

Note that the parameters must be in the same order that they would be if you
were entering them in XSPEC. The parameters ``emin``, ``emax``, ``nbins``,
and ``binscale`` are used to control the binning.

The second way involves passing an XSPEC script file to the
:meth:`~soxs.spectra.Spectrum.from_xspec_script` method which defines an XSPEC
model. For example, a script that creates a model spectrum from a sum of two
APEC models may look like this:

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

If it is contained within the file ``"two_apec.xcm"``, it can be used to
create a :class:`~soxs.spectra.Spectrum` like this:

.. code-block:: python

    emin = 0.1
    emax = 5.0
    nbins = 20000
    spec = Spectrum.from_xspec_script("two_apec.xcm", emin, emax, nbins,
                                      binscale="log")

The parameters ``emin``, ``emax``, ``nbins``, and ``binscale`` are used to
control the binning.

.. note::

    Generating spectra from XSPEC requires that the ``HEADAS`` environment variable
    is defined within your shell before running the Python script, as it would be
    if you were using XSPEC to fit spectra. For example, for the ``zsh`` shell there
    should be a line like ``export HEADAS=${HOME}/heasoft-6.29/x86_64-apple-darwin21.1.0/``
    in your ``.zshrc`` file.

Math with ``Spectrum`` Objects
------------------------------

Two :class:`~soxs.spectra.Spectrum` objects can be co-added, provided that
they have the same energy binning:

.. code-block:: python

    spec1 = Spectrum.from_powerlaw(1.1, 0.05, 1.0e-9, 0.1, 10.0, 10000)
    spec2 = agen.get_spectrum(6.0, 0.3, 0.05, 1.0e-3)

    total_spectrum = spec1 + spec2

If they do not, an error will be thrown.

Or they can be subtracted:

.. code-block:: python

    diff_spectrum = spec1-spec2

You can also multiply a spectrum by a constant float number or divide it by one:

.. code-block:: python

    spec3 = 6.0*spec2
    spec4 = spec1/4.4

.. _band-ops:

Getting the Values and Total Flux or Luminosity of a Spectrum Within a Specific Energy Band
-------------------------------------------------------------------------------------------

A new :class:`~soxs.spectra.Spectrum` object can be created from a restricted
energy band of an existing one by calling the :meth:`~soxs.spectra.Spectrum.new_spec_from_band`
method:

.. code-block:: python

    emin = 0.5
    emax = 7.0
    subspec = spec.new_spec_from_band(emin, emax)

The :meth:`~soxs.spectra.Spectrum.get_flux_in_band` method can be used
to quickly report on the total flux within a specific energy band within
the observer frame:

.. code-block:: python

    emin = 0.5
    emax = 7.0
    print(spec.get_flux_in_band(emin, emax))

which returns a tuple of the photon flux and the energy flux, showing:

.. code-block:: pycon

    (<Quantity 2.2215588675210208e-07 ph / (cm2 s)>,
     <Quantity 7.8742710307246895e-16 erg / (cm2 s)>)

The :meth:`~soxs.spectra.Spectrum.get_lum_in_band` method can also be used
to quickly report on the total luminosity and count rate within a specific
energy band, where in this case the band in question is the rest frame of
the source. For this reason, either a redshift must be supplied, or for a
local source a distance must be given.

.. code-block:: python

    emin = 0.5
    emax = 7.0
    print(spec.get_lum_in_band(emin, emax, redshift=0.05))

which returns a tuple of the photon count rate and the luminosity, showing:

.. code-block:: pycon

    (<Quantity 1.35081761e+48 ph / s>, <Quantity 4.78819407e+39 erg / s>)

You can change the cosmology as well by supplying a :class:`~astropy.cosmology.Cosmology`
object to ``cosmology`` (otherwise the Planck 2018 cosmology is assumed):

.. code-block:: python

    from astropy.cosmology import WMAP9
    emin = 0.5
    emax = 7.0
    print(spec.get_lum_in_band(emin, emax, redshift=0.05, cosmology=WMAP9))

See the `AstroPy cosmology documentation <https://docs.astropy.org/en/stable/cosmology/index.html>`_
for more details.

You can supply a distance for a local source (redshift assumed zero) like this:

.. code-block:: python

    emin = 0.5
    emax = 7.0
    print(spec.get_lum_in_band(emin, emax, dist=(8.0, "kpc")))

Finally, :class:`~soxs.spectra.Spectrum` objects are "callable", and if one
supplies a single energy or array of energies, the values of the spectrum
at these energies will be returned. AstroPy :class:`~astropy.units.Quantity`
objects are detected and handled appropriately.

.. code-block:: python

    print(spec(3.0)) # energy assumed to be in keV

.. code-block:: pycon

    <Quantity 2.830468922349541e-10 ph / (cm2 keV s)>

.. code-block:: python

    from astropy.units import Quantity
    # AstroPy quantity, units will be converted to keV internally
    e = Quantity([1.6e-9, 3.2e-9, 8.0e-9], "erg")
    print(spec(e)) # energy assumed to be in keV

.. code-block:: pycon

    <Quantity [  9.47745587e-10,  4.42138950e-10,  1.61370731e-10] ph / (cm2 keV s)>

Rescaling the Normalization of a Spectrum
-----------------------------------------

You can rescale the normalization of the entire spectrum using the
:meth:`~soxs.spectra.Spectrum.rescale_flux` method. This can be
helpful when you want to set the normalization of the spectrum by the
total flux within a certain energy band instead.

.. code-block:: python

    spec.rescale_flux(1.0e-9, emin=0.5, emax=7.0, flux_type="photons"):

``emin`` and ``emax`` can be used to set the band that the flux corresponds to.
If they are not set, they are assumed to be the bounds of the spectrum. The flux
type can be ``"photons"`` (the default) or ``"energy"``. In the former case, the
units of the new flux must be :math:`{\rm photons}~{\rm cm}^{-2}~{\rm s}^{-1}`,
and in the latter case the units must be
:math:`{\rm erg}~{\rm cm}^{-2}~{\rm s}^{-1}`.

.. _galactic_abs:

Applying Galactic Foreground Absorption to a Spectrum
-----------------------------------------------------

The :meth:`~soxs.spectra.Spectrum.apply_foreground_absorption` method
can be used to apply foreground absorption using the ``"wabs"`` or
``"tbabs"`` models. It takes one required parameter, the hydrogen
column along the line of sight, in units of :math:`10^{22}~{\rm cm}^{-2}`.
Once can optionally specify which absorption model to use using the ``"model"``
parameter (default is ``"wabs"``):

.. code-block:: python

    spec = Spectrum.from_powerlaw(1.1, 0.05, 1.0e-9, 0.1, 10.0, 10000)
    n_H = 0.02
    spec.apply_foreground_absorption(n_H, model="tbabs")

The flux in the energy bins will be reduced according to the absorption at a
given energy. Optionally, to model absorption intrinsic to a source or
from a source intermediate between us and the source, one can supply an
optional ``redshift`` argument (default 0.0):

.. code-block:: python

    spec = Spectrum.from_powerlaw(1.1, 0.05, 1.0e-9, 0.1,
                                  10.0, 10000)
    n_H = 0.02
    spec.apply_foreground_absorption(n_H, model="tbabs", redshift=0.05)

Finally, the abundance table for the ``"tbabs"`` absorption model can be
specified (the default is ``"angr"``):

.. code-block:: python

    spec = Spectrum.from_powerlaw(1.1, 0.05, 1.0e-9, 0.1,
                                  10.0, 10000)
    n_H = 0.02
    spec.apply_foreground_absorption(n_H, model="tbabs", redshift=0.05,
                                     abund_table="wilm")

See :ref:`solar-abund-tables` for options for different abundance tables.

The current version for the ``"tbabs"`` model is 2.3.2.

.. _emiss_lines:

Adding Emission Lines to a Spectrum
-----------------------------------

The :meth:`~soxs.Spectrum.add_emission_line` method adds a single Gaussian
emission line to an existing :class:`~soxs.spectra.Spectrum` object. The
line energy, line width, and amplitude of the line (the line strength or
integral under the curve) must be specified. The formula for the emission
line is:

.. math::

    f(E) = \frac{A}{\sqrt{2\pi\sigma^2}}\exp{\left[-\frac{(E-E_0)^2}{2\sigma^2}\right]}

where :math:`E_0` is the line center and the line width is

.. math::

    {\rm FWHM} = 2\sqrt{2\ln{2}}\sigma

.. code-block:: python

    spec = Spectrum.from_powerlaw(1.1, 0.05, 1.0e-9, 0.1,
                                  10.0, 10000)
    line_center = (6.0, "keV") # "E_0" above
    line_width = (30.0, "eV") # "FWHM" above
    line_amp = (1.0e-7, "photon/s/cm**2") # "A" above
    spec.add_emission_line(line_center, line_width, line_amp)

The line width may also be specified in units of velocity, if that is more convenient:

.. code-block:: python

    spec = Spectrum.from_powerlaw(1.1, 0.05, 1.0e-9, 0.1,
                                  10.0, 10000)
    line_center = (6.0, "keV")
    line_width = (200.0, "km/s")
    line_amp = (1.0e-7, "photon/s/cm**2")
    spec.add_emission_line(line_center, line_width, line_amp)

Currently, this functionality only supports emission lines with a Gaussian shape.

.. _absorb_lines:

Adding Absorption Lines to a Spectrum
-------------------------------------

The :meth:`~soxs.Spectrum.add_absorption_line` method adds a single Gaussian
absorption line to an existing :class:`~soxs.spectra.Spectrum` object. The
line energy, line width, and equivalent width of the line must be specified.
The formula for the absorption line is given in terms of the optical depth
:math:`\tau(E)`:

.. math::

    \tau(E) = \frac{B}{\sqrt{2\pi\sigma^2}}\exp{\left[-\frac{(E-E_0)^2}{2\sigma^2}\right]}

where :math:`E_0` is the line center and the line width is

.. math::

    {\rm FWHM} = 2\sqrt{2\ln{2}}\sigma

and the strength of the absorption :math:`B` is

.. math::

    B = E_0^2\frac{\rm EW}{hc}

where :math:`{\rm EW}` is the equivalent width in angstroms. Then the unabsorbed
spectrum :math:`f_0(E)` is multiplied by the absorption like so to produce the
absorbed spectrum :math:`f(E)`:

.. math::

    f(E) = e^{-\tau(E)}f_0(E)

.. code-block:: python

    spec = Spectrum.from_powerlaw(1.1, 0.05, 1.0e-9, 0.1,
                                  10.0, 10000)
    line_center = (1.0, "keV") # "E_0" above
    line_width = (30.0, "eV") # "FWHM" above
    equiv_width = 2 # defaults to units of milli-Angstroms
    spec.add_absorption_line(line_center, line_width, equiv_width)

The line width may also be specified in units of velocity, if that is more convenient:

.. code-block:: python

    spec = Spectrum.from_powerlaw(1.1, 0.05, 1.0e-9, 0.1,
                                  10.0, 10000)
    line_center = (1.0, "keV")
    line_width = (500.0, "km/s")
    equiv_width = (3.0e-3, "Angstrom")
    spec.add_absorption_line(line_center, line_width, equiv_width)

Currently, this functionality only supports absorption lines with a Gaussian shape.

Generating Photon Energies From a Spectrum
------------------------------------------

Given a :class:`~soxs.spectra.Spectrum`, a set of photon energies can be
drawn from it using the :meth:`~soxs.spectra.Spectrum.generate_energies`
method. This will most often be used to generate discrete samples for mock
observations. For this method, an exposure time and a constant
(energy-independent) effective area must be supplied to convert the spectrum's
flux to a number of photons. These values need not be realistic--in fact, they
both should be larger than the values for the mock observation that you want to
simulate, to create a statistically robust sample to draw photons from when we
actually pass them to the instrument simulator.

An example using a :class:`~soxs.spectra.Spectrum` created from a file:

.. code-block:: python

    spec = Spectrum.from_file("my_spec.dat")
    t_exp = (100., "ks") # exposure time
    area = (3.0, "m**2") # constant effective area
    energies = spec.generate_energies(t_exp, area)

The ``energies`` object :meth:`~soxs.spectra.Spectrum.generate_energies` returns
is an augmented NumPy array which also carries the unit information and the total
flux of energies:

.. code-block:: python

    print(energies.unit)
    print(energies.flux)

.. code-block:: pycon

    Unit("keV")
    <Quantity 1.1256362913845828e-15 erg / (cm2 s)>

Normally, :meth:`~soxs.spectra.Spectrum.generate_energies` will not need to be
called by the end-user but will be used "under the hood" in the generation of
a :class:`~soxs.simput.PhotonList` as part of a :class:`~soxs.simput.SimputCatalog`.
See :ref:`simput` for more information.

.. _count-rate-spectra:

Count Rate Spectra
------------------

The :class:`~soxs.spectra.CountRateSpectrum` class is basically the same thing as a
the :class:`~soxs.spectra.Spectrum` class, except that it is in units of
:math:`\rm{counts}~\rm{s}^{-1}~\rm{keV}^{-1}`. This sort of spectrum makes the most
sense in the rest frame of a source. This object is usually not generated on its own,
but is the result of some other kind of operation (such as
`making source spectra in pyXSIM <https://hea-www.cfa.harvard.edu/~jzuhone/pyxsim/spectra.html>`_).

One important note about :class:`~soxs.spectra.CountRateSpectrum` objects is that you
can also call :meth:`~soxs.spectra.CountRateSpectrum.generate_energies` on them, except
that unlike :class:`~soxs.spectra.Spectrum` objects it is not necessary to specify an area,
but only an exposure time, to generate energies:

.. code-block:: python

    # here "spec" is a CountRateSpectrum object
    t_exp = (100., "ks") # exposure time
    energies = spec.generate_energies(t_exp)

.. _convolved-spectra:

"Convolved" Spectra
-------------------

One may want to examine a spectrum after it has been convolved with a particular
effective area curve. One can generate such a
:class:`~soxs.spectra.ConvolvedSpectrum` using the
:meth:`~soxs.spectra.ConvolvedSpectrum.convolve` method, feeding it a
:class:`~soxs.spectra.Spectrum` object and an ARF:

.. code-block:: python

    from soxs import ConvolvedSpectrum
    # Assuming one created an ApecGenerator agen...
    spec2 = agen.get_spectrum(6.0, 0.3, 0.05, 1.0e-3)
    cspec = ConvolvedSpectrum.convolve(spec2, "xrs_hdxi_3x10.arf")

The spectrum in this object has units of
:math:`{\rm photons}~{\rm s}^{-1}~{\rm keV}^{-1}`, and one can use many of
:class:`~soxs.spectra.Spectrum`'s methods on it. For example, to determine the
count and energy rate within a particular band:

.. code-block:: python

    cspec.get_flux_in_band(0.5, 7.0)

.. code-block:: python

    (<Quantity 6.802363401824924 ph / s>,
     <Quantity 1.2428592072628134e-08 erg / s>)

Or to generate an array of energies:

.. code-block:: python

    t_exp = (500.0, "ks")
    e = cspec.generate_energies(t_exp)

If one has already loaded a :class:`~soxs.instrument.AuxiliaryResponseFile`,
then one can also generate a :class:`~soxs.spectra.ConvolvedSpectrum` by simply
multiplying the ARF by a :class:`~soxs.spectra.Spectrum` object:

.. code-block:: python

    from soxs import AuxiliaryResponseFile
    arf = AuxiliaryResponseFile("xrs_hdxi_3x10.arf")
    # Assuming one created an ApecGenerator agen...
    spec2 = agen.get_spectrum(6.0, 0.3, 0.05, 1.0e-3)
    cspec = spec2*arf

To "deconvolve" a :class:`~soxs.spectra.ConvolvedSpectrum` object and return
a :class:`~soxs.spectra.Spectrum` object, simply call
:meth:`~soxs.spectra.ConvolvedSpectrum.deconvolve`:

.. code-block:: python

    spec_new = cspec.deconvolve()

.. _spectra-plots:

Plotting Spectra
----------------

All :class:`~soxs.spectra.Spectrum` objects and their associated subclasses have
a :meth:`~soxs.spectra.Spectrum.plot` method which can be used to make a
`Matplotlib <http://www.matplotlib.org>`_ plot. The :meth:`~soxs.spectra.Spectrum.plot`
method has no required arguments, but has a number of optional arguments for plot
customization. This method returns a tuple of the :class:`~matplotlib.figure.Figure` and
the :class:`~matplotlib.axes.Axes` objects to allow for further customization. This
example shows how to make a simple plot of an absorbed power-law spectrum:

.. code-block:: python

    spec = soxs.Spectrum.from_powerlaw(1.2, 0.02, 1.0e-3, 0.2, 9.0, 100000)
    spec.apply_foreground_absorption(0.1)
    fig, ax = spec.plot()

.. image:: ../images/plot_powerlaw.png

Here's another example of creating a plot of two thermal spectra with labels,
zooming in on a section of it, and setting the energy scale to linear:

.. code-block:: python

    agen = soxs.ApecGenerator(0.1, 10.0, 10000)
    spec1 = agen.get_spectrum(5.0, 0.3, 0.02, 1.0e-3)
    spec2 = agen.get_spectrum(3.0, 0.3, 0.02, 1.0e-3)
    fig, ax = spec1.plot(xmin=0.7, xmax=1.5, ymin=1.0e-4, ymax=3.0e-3,
                         xscale='linear', label="5 keV plasma")
    spec2.plot(fig=fig, ax=ax, label="3 keV plasma")

.. image:: ../images/plot_two_spectra.png

For other customizations, consult the :meth:`~soxs.spectra.Spectrum.plot` API.

.. _write-spectra:

Writing a Spectrum to Disk
--------------------------

:class:`~soxs.spectra.Spectrum` objects can be written to disk in three formats:
an ASCII text file in the ECSV format, a FITS file, or an HDF5 file. To write a
spectrum to an ASCII ECSV file, use the :meth:`~soxs.spectra.Spectrum.write_ascii_file`
method:

.. code-block:: python

    agen = soxs.ApecGenerator(0.1, 10.0, 10000)
    spec1 = agen.get_spectrum(5.0, 0.3, 0.02, 1.0e-3)
    spec1.write_ascii_file("my_spec.ecsv", overwrite=True)

To write a spectrum to an HDF5 file, use :meth:`~soxs.spectra.Spectrum.write_hdf5_file`:

.. code-block:: python

    agen = soxs.ApecGenerator(0.1, 10.0, 10000)
    spec1 = agen.get_spectrum(5.0, 0.3, 0.02, 1.0e-3)
    spec1.write_hdf5_file("my_spec.h5", overwrite=True)

To write a spectrum to a FITS file, use :meth:`~soxs.spectra.Spectrum.write_fits_file`:

.. code-block:: python

    agen = soxs.ApecGenerator(0.1, 10.0, 10000)
    spec1 = agen.get_spectrum(5.0, 0.3, 0.02, 1.0e-3)
    spec1.write_fits_file("my_spec.fits", overwrite=True)

In each case, the minimum and maximum energies for each bin in the table, the
flux in each bin (as well as its units), and the bin scaling (linear or log)
is written to the file. If writing a :class:`~soxs.spectrum.ConvolvedSpectrum`
object, the name of the ARF which was used to do the convolution is also stored.

.. _read-spectra:

Reading a Spectrum from Disk
----------------------------

:class:`~soxs.spectra.Spectrum` objects written using any of the writing methods
detailed above (ASCII ECSV, HDF5, or FITS) can be the spectrum can be read back
in again in, using :meth:`~soxs.spectra.Spectrum.from_file`:

.. code-block:: python

    from soxs import Spectrum
    my_spec = Spectrum.from_file("my_spec.ecsv")
