.. _thermal-spectra:

Generating Thermal Spectra
==========================

SOXS provides a number of models for generating thermal spectra. These models are
encapsulated in the notion of "generators", which are objects designed to set up
the general characteristics of spectra for a given model, and then allow one to
create spectra with particular parameters. These spectra are for hot plasmas in
collisional ionization equilibrium (CIE), non-equilibrium ionization (NEI), or
a combination of collisional and photoionization processes. The various thermal
spectrum generators available in SOXS are:

* :class:`~soxs.thermal_spectra.ApecGenerator`: APEC CIE and NEI spectra
* :class:`~soxs.thermal_spectra.SpexGenerator`: SPEX CIE spectra
* :class:`~soxs.thermal_spectra.MekalGenerator`: MeKaL CIE spectra
* :class:`~soxs.thermal_spectra.CloudyCIEGenerator`: Cloudy-based CIE spectra
* :class:`~soxs.thermal_spectra.IGMGenerator`: Cloudy-based collisional+photoionization
  spectra with optional resonant scattering from the CXB

Each of these are essentially "factory" classes which generate new
:class:`~soxs.spectra.Spectrum` objects.

.. _apec-spectra:

APEC Spectrum Generators
------------------------

APEC CIE and NEI thermal spectra are generated in SOXS using the
`AtomDB tables <http://www.atomdb.org>`_, using the
:class:`~soxs.thermal_spectra.ApecGenerator` object:

.. code-block:: python

    from soxs import ApecGenerator
    emin = 0.05
    emax = 50.0
    nbins = 10000
    binscale = "linear"
    agen = ApecGenerator(emin, emax, nbins, binscale=binscale,
                         apec_vers="2.0.2", broadening=True)

The parameters ``emin``, ``emax``, ``nbins``, and ``binscale`` are used to
control the binning.

The ``broadening`` parameter sets whether or not spectral lines will be
thermally and velocity broadened, and is ``True`` by default..The
``apec_vers`` parameter sets the version of the AtomDB tables to use.
Version 3.0.9 is the default, and the tables will be downloaded if necessary.

You may also supply another location for the AtomDB tables. For example, the
following construction will look for the AtomDB tables in the current working
directory:

.. code-block:: python

    agen = ApecGenerator(0.05, 50.0, 10000, apec_root=".")

If you set the ``apec_vers`` parameter but not the ``apec_root`` parameter, the
AtomDB table files will be looked for in (1) the current working directory and
(2) the location specified by ``soxs_data_dir`` in the :ref:`config`.

The range of temperatures supported by the APEC model is kT = .

.. _nolines:

APEC Spectra Without Lines
++++++++++++++++++++++++++

There is also an option to generate continuum spectra without line emission.
This is done by setting ``nolines=True`` in the constructor for
:class:`~soxs.thermal_spectra.ApecGenerator`:

.. code-block:: python

    agen = ApecGenerator(0.05, 50.0, 10000, nolines=True)

.. _spex-spectra:

SPEX Spectrum Generators
------------------------

Thermal spectra using the
`CIE emission model <https://spex-xray.github.io/spex-help/models/cie.html>`_
provided in `SPEX <https://www.sron.nl/astrophysics-spex>`_ can be generated
using the :class:`~soxs.thermal_spectra.SpexGenerator` class. The same
underlying machinery as the APEC model is used, as the SPEX model has been
converted to the APEC table format using the code at
https://github.com/jeremysanders/spex_to_xspec. As such, this class takes the
same arguments as :class:`~soxs.thermal_spectra.ApecGenerator`, with the
exception that the version and file location arguments are named ``spex_vers``
and ``spex_root``, respectively.

.. code-block:: python

    from soxs import SpexGenerator

    sgen = SpexGenerator(0.05, 50.0, 10000, binscale="log")

If you set the ``spex_vers`` parameter but not the ``spex_root`` parameter, the
AtomDB table files will be looked for in (1) the current working directory and
(2) the location specified by ``soxs_data_dir`` in the :ref:`config`. The current
default version of the SPEX thermal model in SOXS is 3.06.01.

The range of temperatures supported by the SPEX model is kT = .

.. _mekal-spectra:

MeKaL Spectrum Generators
-------------------------

MeKaL is an X-ray emission model for a hot, diffuse, thermal plasma
in CIE based on calculations by Mewe and Kaastra with Fe L calculations
by Liedahl. Relevant references are:

* https://ui.adsabs.harvard.edu/abs/1985A%26AS...62..197M
* https://ui.adsabs.harvard.edu/abs/1986A%26AS...65..511M
* https://ui.adsabs.harvard.edu/abs/1995ApJ...438L.115L

Thermal CIE spectra using the MeKaL model can be generated using the
:class:`~soxs.thermal_spectra.MekalGenerator` object, in a manner
similar to the APEC and SPEX generators:

.. code-block:: python

    emin = 0.1
    emax = 10.0
    nbins = 3000
    mgen = soxs.MekalGenerator(emin, emax, nbins, binscale="linear")

Note that it is not possible to create MeKaL spectra without lines
or with thermal broadening. The range of temperatures supported by the
MeKaL model is kT = .

.. _cloudy-spectra:

Cloudy CIE Spectrum Generators
------------------------------

The :class:`~soxs.thermal_spectra.CloudyCIEGenerator` generates CIE
thermal spectra using the emission model from
`Cloudy <https://gitlab.nublado.org/cloudy/cloudy/-/wikis/home>`_, which
are interpolated from an
`XSPEC atable <https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/general/ogip_92_009/ogip_92_009.html>`_.
The sequence of Cloudy commands used to generate the XSPEC atable is
as follows:

.. code-block::

    #########
    title cie_v4_lo_tgrid
    #
    #database stout level MAX
    database chianti level MAX
    no molecules
    no grain physics
    set phfit 1996
    abundances "./feld.abn"
    ###
    metals 0.0
    hden 0
    #
    coronal equil 5 vary
    grid range 4.0 to 9.0 in 0.025 dex steps sequential
    set continuum resolution 0.1
    stop column density 1.5032e+18 linear
    iterate to convergence
    save xspec atable reflected spectrum "cie_v4_lo_tgrid_n1z1.fits" range 0.05 10.
    #########

This sequence of commands is repeated for solar and low abundances so that
the abundance parameter can be taken into account via a linear combination of
two tables. For the individual abundances, they are obtained by setting e.g.
"element neon off" in the run and doing the appropriate arithmetic.

Thermal CIE spectra using the Cloudy model can be generated using the
:class:`~soxs.thermal_spectra.CloudyCIEGenerator` object, in a manner
similar to the other generators:

.. code-block:: python

    emin = 0.1
    emax = 9.0
    nbins = 3000
    cgen = soxs.CloudyCIEGenerator(emin, emax, nbins, binscale="linear")

The energy range of the generated table at a redshift of 0 is 0.05-10.0 keV. The
default resolution of Cloudy in this range is dE/E = 0.005 for E < 8.16 keV, and
dE/E = 0.03 for higher energies. SOXS provides two different spectral resolutions
for the Cloudy CIE tables, the lower of which (and the default) is 10 times finer
than the above and the higher is 40 times finer. Which is used can be controlled
by the *model_vers* argument. To specify the higher resolution table:

.. code-block:: python

    emin = 0.1
    emax = 9.0
    nbins = 3000
    cgen = soxs.CloudyCIEGenerator(emin, emax, nbins, binscale="log",
                                   model_vers="4_hi")

Note that it is not possible to create Cloudy CIE spectra without lines, and
thermal broadening is automatically included and cannot be turned off. The
range of temperatures supported by the Cloudy CIE model is :math:`T = 10^4-10^9` K.

.. _igm-spectra:

IGM Spectrum Generators
-----------------------

The :class:`~soxs.thermal_spectra.IGMGenerator` generates thermal
X-ray emission spectra from a photoionized and collisionally ionized
plasma, as well as resonant scattering by the CXB, based on
`Khabibullin & Churazov 2019 <https://ui.adsabs.harvard.edu/abs/2019MNRAS.482.4972K/>`_
and `Churazov et al. 2001 <https://ui.adsabs.harvard.edu/abs/2001MNRAS.323...93C/>`_.
Because this model includes photoionization and (optionally) resonant
scattering of the CXB, it is density-dependent. It is intended to be used
primarily for simulations of spectra from low-density, sub-keV temperature
plasmas such as the warm-hot intergalactic medium (WHIM), or the low-density
parts of the circumgalactic medium (CGM). For resonant scattering, it is assumed
that a fraction of CXB photons are scattering off of heavy ions, enhancing line
emission.

IGM spectra using the Cloudy model can be generated using the
:class:`~soxs.thermal_spectra.IGMGenerator` object, in a manner
similar to the other generators:

.. code-block:: python

    emin = 0.1
    emax = 5.0
    nbins = 3000
    igen = soxs.IGMGenerator(emin, emax, nbins, binscale="log")

If you want to include the effects of resonant scattering off of CXB photons,
you must set ``resonant_scattering=True``. Optionally, you may also
change the fraction of the CXB that is scattered from the ions using the
``cxb_factor`` parameter, but ideally this should remain at the default value
of 0.5:

.. code-block:: python

    emin = 0.1
    emax = 10.0
    nbins = 3000
    igen = soxs.IGMGenerator(emin, emax, nbins, binscale="linear",
                             resonant_scattering=True, cxb_factor=0.3)

The energy range of the generated table at a redshift of 0 is 0.05-10.0 keV. The
default resolution of Cloudy in this range is dE/E = 0.005 for E < 8.16 keV, and
dE/E = 0.03 for higher energies. SOXS provides two different spectral resolutions
for the IGM tables, the lower of which (and the default) is 10 times finer
than the above and the higher is 40 times finer. Which is used can be controlled
by the *model_vers* argument. To specify the higher resolution table:

.. code-block:: python

    emin = 0.1
    emax = 9.0
    nbins = 3000
    igen = soxs.IGMGenerator(emin, emax, nbins, binscale="log",
                             model_vers="4_hi")

Generating Spectra
------------------

Once you have a generator object, you can use it to generate thermal
spectra using ``get_spectrum`` method which is available for each model.
For the CIE generators, the general signature looks like this:

.. code-block:: python

    spec = sgen.get_spectrum(kT, abund, redshift, norm, velocity=velocity)

and the parameters are:

* ``kT``: The temperature of the plasma, with default units of keV
* ``abund``: The metal abundance, in solar units.
  See :ref:`var-abund` below for more fine-grained control of abundances.
* ``redshift``: The redshift of the plasma
* ``norm``: The normalization of the model, assuming the standard prescription of
  :math:`10^{-14}\int{n_en_p}dV/[4\pi(1+z)^2D_A^2]` where :math:`n_e` and
  :math`n_p` are the electron and proton number densities, :math:`z` is the
  redshift, and :math:`D_A` is the angular diameter distance to the source. All
  units are in cgs.
* ``velocity``: The (optional) velocity broadening parameter, in units of km/s.
  If not zero, this broadens spectral lines using a Gaussian model assuming the
  ``velocity`` parameter is the velocity dispersion :math:`\sigma_v`. If not set,
  there is no velocity broadening. Currently only available for the
  :class:`~soxs.thermal_spectra.ApecGenerator` and
  :class:`~soxs.thermal_spectra.SpexGenerator` classes.

A more specific invocation may look like this:

.. code-block:: python

    kT = (6.0, "keV")
    abund = 0.3 # solar units
    redshift = 0.05
    norm = 1.0e-3
    velocity = (100.0, "km/s") # optional
    spec1 = agen.get_spectrum(kT, abund, redshift, norm, velocity=velocity)

``spec1`` is just a standard :class:`~soxs.spectra.Spectrum` object.

For the IGM generator, because it includes the effects of photoionization,
it also depends on the hydrogen number density, and the signature looks
like this:

.. code-block:: python

    spec = igen.get_spectrum(kT, nH, abund, redshift, norm)

Where the ``nH`` parameter is the number density of hydrogen atoms in units
of :math:`cm^{-3}`.

.. _var-abund:

Abundance Settings
------------------

By default, the various generators handle abundances greater than H and
He using the ``abund`` parameter in the various ``get_spectrum`` methods.
Exactly what abundances are set by this parameter depends on the model used:

* APEC and SPEX: Includes C, N, O, Ne, Mg, Al, Si, S, Ar, Ca, Fe, Ni (He
  fixed at abundance table value, Li, Be, B, F, Na, P, Cl, K, Sc, Ti,
  V, Cr, Mn, Co, Cu, Zn fixed at solar).
* MeKaL: Includes C, N, O, Ne, Na, Mg, Al, Si, S, Ar, Ca, Fe, Ni (He
  fixed at abundance table value, other elements not modeled)
* Cloudy CIE and IGM: He fixed at abundance table value, all higher
  elements up Zn to included.

More fine-grained control of individual elements is possible. All of the
generators accept a ``var_elem`` optional argument to specify which
elements should be allowed to vary freely:

.. code-block:: python

    var_elem = ["O", "Ca"] # allow oxygen and calcium to vary freely
    agen = ApecGenerator(0.05, 50.0, 10000, var_elem=var_elem, binscale="log")

Whatever elements are not specified here are assumed to be set as normal,
whether they are H, He, trace elements, or metals covered by the ``abund``
parameter.

.. note::

    For the APEC, SPEX, and MeKaL models, any of the elements listed above
    can be specified as variable. For the Cloudy CIE and IGM models, only
    the elements C, N, O, Ne, Fe, S, Si, Ca, and Mg can be variable.

Now, spectra which are created from this generator object using its
``get_spectrum`` method should set values for the abundances of these elements
in solar units. This is done by supplying the ``elem_abund`` dict like so:

.. code-block:: python

    kT = 6.0
    abund = 0.3 # for all other metals
    redshift = 0.05
    norm = 1.0e-3
    O_abund = 0.5
    Ca_abund = 0.4
    spec = agen.get_spectrum(kT, abund, redshift, norm,
                             elem_abund={"O": O_abund, "Ca": Ca_abund})

Note that setting the ``abund`` parameter is still necessary for the other
metals.

All abundances are relative to the
`Anders & Grevesse 1989 <http://adsabs.harvard.edu/abs/1989GeCoA..53..197A>`_
tables. See :ref:`changing-abund-tables` to see how to use a different table.

.. _nei:

Non-Equilibrium Ionization Spectra with APEC
++++++++++++++++++++++++++++++++++++++++++++

A variation on specifying variable abundances in SOXS allows one to construct
non-equilibrium ionization (NEI) spectra. In this case, all ions one desires to
contribute to the spectrum must be put in by hand, with the exception of H and
He, which may be specified, but if they are not they are assumed to be fully
ionized at their Solar abundances.

To create an :class:`~soxs.thermal_spectra.ApecGenerator` object which produces
NEI spectra, one must specify not only the elements one wants but also their
ionization states. The notation is to represent an ion by the element first,
followed by the ``^`` symbol, followed by its ionization state. So for oxygen,
:math:`O^{+1}` would correspond to ``"O^1"``, and so on. The keyword argument
``nei=True`` must also be set. An example using four oxygen ions and two
nitrogen ions is shown below:

.. code-block:: python

    var_elem = ["O^1", "O^2", "O^3", "O^4", "N^4", "N^5"]
    agen = ApecGenerator(0.05, 10.0, 10000, var_elem=var_elem, nei=True)

Once this has been created, we use a special method for NEI spectra,
:meth:`~soxs.thermal_spectra.ApecGenerator.get_nei_spectrum`

.. code-block:: python

    kT = 5.0
    norm = 1.0e-3
    redshift = 0.0
    elem_abund = {"O^1": 0.3, "O^2": 0.5, "O^3": 0.2, "O^4": 0.5,
                  "N^4": 0.2, "N^5": 0.4}
    spec = agen.get_nei_spectrum(kT, elem_abund, redshift, norm)

.. warning::

    SOXS does not make any assumptions about the correctness of the
    relative ion abundances which you input into
    :meth:`~soxs.thermal_spectra.ApecGenerator.get_nei_spectrum`. It
    assumes you have run a NEI code to determine the correct abundances,
    and only computes the spectrum.

.. warning::

    Generating NEI spectra is not currently possible for any model other
    than the APEC model.

.. _changing-abund-tables:

Changing Abundance Tables
-------------------------

The abundance parameters discussed so far assume abundance of a particular
element or a number of elements relative to the Solar value. Underlying this
are the values of the Solar abundances themselves. By default, SOXS uses the
abundance table from
`Anders & Grevesse 1989 <http://adsabs.harvard.edu/abs/1989GeCoA..53..197A>`_.
However, it is possible to change the Solar abundance table in SOXS for the
:class:`~soxs.thermal_spectra.ApecGenerator`,
:class:`~soxs.thermal_spectra.SpexGenerator`,
and :class:`~soxs.thermal_spectra.MekalGenerator` classes.

The abundance tables included with SOXS are:

* ``"angr"``: `Anders & Grevesse 1989 <http://adsabs.harvard.edu/abs/1989GeCoA..53..197A>`_
* ``"aspl"``: `Asplund et al. 2009 <http://adsabs.harvard.edu/abs/2009ARA%26A..47..481A>`_
* ``"wilm"``: `Wilms et al. 2000 <http://adsabs.harvard.edu/abs/2000ApJ...542..914W>`_
* ``"lodd"``: `Lodders 2003 <http://adsabs.harvard.edu/abs/2003ApJ...591.1220L>`_
* ``"feld"``: `Feldman 1992 <https://ui.adsabs.harvard.edu/abs/1992PhyS...46..202F>`_
* ``"cl17.03"``: The abundances used by default in Cloudy 17.03.

The easiest way to ensure that you always use a particular abundance table
is to set the ``abund_table`` element in the :ref:`config`, like so:

.. code-block:: parsed-literal
    [soxs]
    soxs_data_dir = /Users/jzuhone/Data/soxs
    abund_table = lodd

However, the Solar abundance table can also be changed on-the-fly for the APEC,
SPEX, or MeKaL models like this:

.. code-block:: python

    agen = ApecGenerator(0.05, 50.0, 10000, abund_table="aspl")

Alternatively, one can supply their own abundance table by providing a NumPy
array, list, or tuple of abundances 30 elements in length corresponding to
the Solar abundances relative to hydrogen in the order of H, He, Li, Be, B,
C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar, K, Ca, Sc, Ti, V, Cr, Mn, Fe,
Co, Ni, Cu, and Zn. An example:

.. code-block:: python

    my_abund = np.array([1.00E+00, 8.51E-02, 1.12E-11, 2.40E-11, 5.01E-10,
                         2.69E-04, 6.76E-05, 4.90E-04, 3.63E-08, 8.51E-05,
                         1.74E-06, 3.98E-05, 2.82E-06, 3.24E-05, 2.57E-07,
                         1.32E-05, 3.16E-07, 2.51E-06, 1.07E-07, 2.19E-06,
                         1.41E-09, 8.91E-08, 8.51E-09, 4.37E-07, 2.69E-07,
                         3.16E-05, 9.77E-08, 1.66E-06, 1.55E-08, 3.63E-08])

    agen = ApecGenerator(0.05, 50.0, 10000, abund_table=my_abund)

.. warning::

    It is currently not possible to change the abundance table for either the
    Cloudy CIE or IGM models, as they always use ``"feld"``.

.. warning::

    Although it is possible to specify a custom table of abundances from a
    file for the simulation of thermal spectra, this is not possible for the
    TBabs abundance model used in SOXS--one must instead use one of the
    included options mentioned above. See :ref:`galactic_abs`.

.. _downloading-thermal-tables:

Downloading Thermal Spectra Tables
----------------------------------

SOXS will download the thermal spectra tables necessary to create spectra
using the different models, and store them in the location for SOXS data
specified in the SOXS configuration (see :ref:`config` for more information).
However, if you would like to download data for a specific model to a
particular location, SOXS provides the
:func:`~soxs.thermal_spectra.download_spectrum_tables` function. Specific
model files can be downloaded using the following syntax:

* ``"apec"``: Downloads the latest version of the APEC tables for
  :class:`~soxs.thermal_spectra.APECGenerator`. To download a particular
  version, specify (e.g.) ``"apec_v3.0.9"``, or to download the NEI version
  of the tables, specify ``"apec_v3.0.9_nei"``.
* ``"spex"``: Downloads the latest version of the SPEX tables for
  :class:`~soxs.thermal_spectra.SPEXGenerator`. To download a particular
  version, specify (e.g.) ``"spex_v3.06.01"``.
* ``"cie"``: Downloads the latest version of the low-resolution tables for
  :class:`~soxs.thermal_spectra.CloudyCIEGenerator`. Specify ``"cie_v4_hi"``
  to get the high-resolution version of the tables.
* ``"igm"``: Downloads the latest version of the low-resolution tables for
  :class:`~soxs.thermal_spectra.IGMGenerator`. Specify ``"igm_v4_hi"`` to get
  the high-resolution version of the tables.
* ``"mekal"``: Downloads the MeKaL table.

Download files to the current working directory:

.. code-block:: python

    soxs.download_spectrum_tables("mekal")

or to a location ``loc`` of your choice:

.. code-block:: python

    soxs.download_spectrum_tables("cie_v4_hi", loc="/Users/jzuhone/Data")
