.. _general-info:

General Information Regarding the Python Interface to SOXS
==========================================================

.. _response-path:

Path to SOXS Data Files
-----------------------

To use either :func:`~soxs.instrument.instrument_simulator` or 
:func:`~soxs.instrument.simulate_spectrum`, data files such as the instrumental
responses, background models, and PSF models are required. In versions of SOXS
previous to v3.0.0, it was necessary to download these files on your own and
place them either in the current working directory, or in a location specified
by the :ref:`config`. Now, whenever an instrument is used, SOXS will first 
check the current working directory for the necessary files, and then will 
check the location specified by the ``soxs_data_dir`` entry in the configuration
file. If the files are not found in either location, they will be downloaded
automatically. If ``soxs_data_dir`` is not set in the configuration file, or is
set to an invalid directory, a default directory will be chosen:

.. code-block:: pycon

    soxs : [WARNING  ] 2021-04-14 22:05:49,790 Setting 'soxs_data_dir' to /Users/jzuhone/Library/Caches/soxs for this session. Please update your configuration if you want it somewhere else.

See :ref:`config` for more information about the location of the configuration 
file and how to set its parameters.

.. _units:

Special Argument Handling for Quantities with Units
---------------------------------------------------

Many arguments to functions and class defintions which have units can 
take a special format which allows one to specify that particular
quantity in the units desired by the user. For example, the 
:func:`~soxs.cosmology.make_cosmological_sources_file` function has
several arguments which accept units. If one supplies floating-point
numbers, they will be in a default set of units:

.. code-block:: python

    import soxs
    filename = "cosmo.simput"
    name = "cosmo_srcs"
    sky_center = [30.0, 45.0]
    exp_time = 500000.0 # seconds
    fov = 40.0 # arcmin
    area = 40000.0 # cm^2
    nH = 0.02 # atoms/cm^2
    soxs.make_cosmological_sources_file(filename, name, exp_time, fov, 
                                        sky_center, nH=nH, area=area):

However, these same arguments accept values with unit information, either in the
form of ``(value, unit)`` tuples, :class:`~astropy.units.Quantity`, or
:class:`~yt.units.yt_array.YTQuantity` objects:

.. code-block:: python

    import soxs
    from astropy.units import Quantity
    filename = "cosmo.simput"
    name = "cosmo_srcs"
    sky_center = [30.0, 45.0]
    exp_time = (500.0, "ks")
    fov = Quantity(0.666667, "deg")
    area = (4.0, "m**2") 
    nH = Quantity(2.0e20, "cm**-2") 
    soxs.make_cosmological_sources_file(filename, name, exp_time, fov, 
                                        sky_center, nH=nH, area=area):

Since the quantities are the same but in different units, these two calls would
be equivalent. Check the :ref:`api` for any given function or class definition 
to see which of them have arguments which can take values with specified units, 
and what the default units are.

.. _random-numbers:

Random Number Generation
------------------------

Many routines in SOXS require generating random numbers for energies, sky
positions, spectral channels, etc. By default, for every SOXS run this will
be a different set of random numbers. It is often the case, however, that one
wants to use a consistent, repeatable set of random numbers to reproduce results
exactly. For this, many functions in SOXS take a ``prng`` optional argument, 
which has a default of ``None``, but if set to an integer will use this value as
a random seed. 

For example, to generate photon energies from a :class:`~soxs.spectra.Spectrum`
object using the :meth:`~soxs.spectra.Spectrum.generate_energies` method, one 
would set the random seed in this way:

.. code-block:: python

    t_exp = (50.0, "ks")
    area = (3.0, "m**2")
    prng = 24
    e = spec.generate_energies(t_exp, area, prng=prng)

Check the :ref:`api` to see which functions or methods allow for the input of 
random seeds. 