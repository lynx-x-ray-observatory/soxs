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
* :class:`~soxs.thermal_spectra.CloudyCIEGenerator`: Cloudy CIE spectra
* :class:`~soxs.thermal_spectra.IGMGenerator`: Cloudy-based collisional+photoionization 
  spectra with optional resonant scattering from the CXB

Each of these are essentially "factory" classes which generate new 
:class:`~soxs.spectra.Spectrum` objects.

.. _apec-spectra:

APEC Spectra
------------

APEC CIE and NEI thermal spectra are generated in SOXS using the 
`AtomDB tables <http://www.atomdb.org>`_, using the 
:class:`~soxs.apec.ApecGenerator` object:

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
thermally and velocity broadened. The ``apec_vers`` parameter sets the version 
of the AtomDB tables to use. Version 3.0.9 is the default, and the tables will
be downloaded if necessary. 

You may also supply another location for the AtomDB tables. For example, the 
following construction will look for the AtomDB tables in the current working 
directory:

.. code-block:: python

    agen = ApecGenerator(0.05, 50.0, 10000, apec_root=".")

If you set the ``apec_vers`` parameter but not the ``apec_root`` parameter, the
AtomDB table files will be looked for in (1) the current working directory and
(2) the location specified by ``soxs_data_dir`` in the :ref:`config`.

Once you have an :class:`~soxs.apec.ApecGenerator` object, you can use it to
generate thermal spectra using the :meth:`~soxs.apec.ApecGenerator.get_spectrum`
method. The parameters are:

* ``kT``: The temperature of the plasma, with default units of keV
* ``abund``: The metal abundance, in solar units. Includes C, N, O, Ne, Mg, Al, 
  Si, S, Ar, Ca, Fe, Ni (He fixed at cosmic, other trace elements fixed at solar). 
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
  there is no velocity broadening. 

.. code-block:: python
    
    kT = 6.0 (6.0, "keV")
    abund = 0.3 # solar units
    redshift = 0.05
    norm = 1.0e-3 
    velocity = (100.0, "km/s") # optional
    spec1 = agen.get_spectrum(kT, abund, redshift, norm, velocity=velocity)

``spec1`` is just a standard :class:`~soxs.spectra.Spectrum` object.

.. _var-abund:

Variable Abundances
+++++++++++++++++++

By default, :class:`~soxs.apec.ApecGenerator` assumes all abundances besides
H, He, and the trace elements are set to the value provided by the ``abund``
parameter. However, more fine-grained control is possible. 
:class:`~soxs.apec.ApecGenerator` accepts a ``var_elem`` optional argument
to specify which elements should be allowed to vary freely:

.. code-block:: python

    var_elem = ["O", "Ca"] # allow oxygen and calcium to vary freely 
    agen = ApecGenerator(0.05, 50.0, 10000, var_elem=var_elem, binscale="log")
    
Whatever elements are not specified here are assumed to be set as normal, whether
they are H, He, trace elements, or metals covered by the ``abund`` parameter. 
Now, spectra which are created from this :class:`~soxs.apec.ApecGenerator`
object using the :meth:`~soxs.apec.ApecGenerator.get_spectrum` method should 
set values for the abundances of these elements in solar units. This is done by 
supplying the ``elem_abund`` dict like so:

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

.. _nei:

Non-Equilibrium Ionization Spectra with APEC
++++++++++++++++++++++++++++++++++++++++++++

A variation on specifying variable abundances in SOXS allows one to construct
non-equilibrium ionization (NEI) spectra. In this case, all ions one desires to
contribute to the spectrum must be put in by hand, with the exception of H and
He, which may be specified, but if they are not they are assumed to be fully
ionized at their Solar abundances.

To create an :class:`~soxs.apec.ApecGenerator` object which produces NEI
spectra, one must specify not only the elements one wants but also their 
ionization states. The notation is to represent an ion by the element first, 
followed by the ``^`` symbol, followed by its ionization state. So for oxygen,
:math:`O^{+1}` would correspond to ``"O^1"``, and so on. The keyword argument 
``nei=True`` must also be set. An example using four oxygen ions and two 
nitrogen ions is shown below:

.. code-block:: python

    var_elem = ["O^1", "O^2", "O^3", "O^4", "N^4", "N^5"]
    agen = ApecGenerator(0.05, 10.0, 10000, var_elem=var_elem, nei=True)

Once this has been created, we use a special method for NEI spectra, 
:meth:`~soxs.apec.ApecGenerator.get_nei_spectrum`

.. code-block:: python

    kT = 5.0 
    norm = 1.0e-3 
    redshift = 0.0
    elem_abund = {"O^1": 0.3, "O^2": 0.5, "O^3": 0.2, "O^4": 0.5,
                  "N^4": 0.2, "N^5": 0.4}
    spec = agen.get_nei_spectrum(kT, elem_abund, redshift, norm)
    
.. warning::

    SOXS does not make any assumptions about the correctness of the relative ion
    abundances which you input into :meth:`~soxs.apec.ApecGenerator.get_nei_spectrum`.
    It assumes you have run a NEI code to determine the correct abundances, and
    only computes the spectrum.

.. _nolines:

APEC Spectra Without Lines
++++++++++++++++++++++++++

There is also an option to generate continuum spectra only from the AtomDB
or SPEX tables. This is done by setting ``nolines=True`` in the constructor for
:class:`~soxs.apec.ApecGenerator`:

.. code-block:: python

    agen = ApecGenerator(0.05, 50.0, 10000, nolines=True)

.. _solar-abund-tables:

Changing Abundance Tables
+++++++++++++++++++++++++

The abundance parameters discussed so far assume abundance of a particular 
element or a number of elements relative to the Solar value. Underlying this
are the values of the Solar abundances themselves. It is possible to change the
Solar abundance table in SOXS via the optional ``abund_table`` argument to 
:class:`~soxs.apec.ApecGenerator`. By default, SOXS uses the abundance table
set in the :ref:`config`, which by default are the
`Anders & Grevesse 1989 <http://adsabs.harvard.edu/abs/1989GeCoA..53..197A>`_ 
abundances. This corresponds to a setting of ``"angr"`` for this parameter, but it 
is possible to use other tables of solar abundances. The other tables included 
with SOXS are:

* ``"aspl"``: `Asplund et al. 2009 <http://adsabs.harvard.edu/abs/2009ARA%26A..47..481A>`_
* ``"wilm"``: `Wilms et al. 2000 <http://adsabs.harvard.edu/abs/2000ApJ...542..914W>`_
* ``"lodd"``: `Lodders 2003 <http://adsabs.harvard.edu/abs/2003ApJ...591.1220L>`_
* ``"feld"``: `Feldman 1992 <https://ui.adsabs.harvard.edu/abs/1992PhyS...46..202F>`_
* ``"cl17.03"``: The abundances used by default in Cloudy 17.03.

The easiest way to ensure that you always use a particular abundance table is to
set it in the :ref:`config`. However, the Solar abundance table can be changed 
on-the-fly like this:

.. code-block:: python

    agen = ApecGenerator(0.05, 50.0, 10000, abund_table="aspl")

Alternatively, one can supply their own abundance table by providing a NumPy array, list,
or tuple of abundances 30 elements in length corresponding to the Solar abundances
relative to hydrogen in the order of H, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P,
S, Cl, Ar, K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, and Zn. An example:

.. code-block:: python

    my_abund = np.array([1.00E+00, 8.51E-02, 1.12E-11, 2.40E-11, 5.01E-10,
                         2.69E-04, 6.76E-05, 4.90E-04, 3.63E-08, 8.51E-05,
                         1.74E-06, 3.98E-05, 2.82E-06, 3.24E-05, 2.57E-07,
                         1.32E-05, 3.16E-07, 2.51E-06, 1.07E-07, 2.19E-06,
                         1.41E-09, 8.91E-08, 8.51E-09, 4.37E-07, 2.69E-07,
                         3.16E-05, 9.77E-08, 1.66E-06, 1.55E-08, 3.63E-08])

    agen = ApecGenerator(0.05, 50.0, 10000, abund_table=my_abund)

.. warning::

    Although it is possible to specify a custom table of abundances from a 
    file for the simulation of thermal spectra, this is not possible for the 
    TBabs abundance model used in SOXS--one must instead use one of the
    included options mentioned above. See :ref:`galactic_abs`.

.. _spex-spectra:

SPEX Spectra
------------

Thermal CIE spectra 
.. _mekal-spectra:

MeKaL Spectra
-------------

.. _cloudy-spectra

Cloudy Spectra
--------------
