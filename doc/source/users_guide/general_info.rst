.. _general-info:

General Information Regarding the Python Interface to SOXS
==========================================================

.. _response-path:

Path to the Response Files
--------------------------

To use either the :func:`~soxs.instrument.instrument_simulator` or 
:func:`~soxs.instrument.simulate_spectrum`, it is necessary to download the 
response files from the :ref:`responses` page and place them in an appropriate
location, of which there are two. The first is simply to place the response 
files needed for the instrument simulator in the current working directory
from which you run SOXS. However, it is probably more convenient to place the
response files in a default path, which can be specified in the SOXS
configuration file like this:

.. code-block:: text

    [soxs]
    response_path = /Users/jzuhone/Data/soxs_responses

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
    simput_prefix = "cosmo"
    phlist_prefix = "cosmo"
    sky_center = [30.0, 45.0]
    exp_time = 500000.0 # seconds
    fov = 40.0 # arcmin
    area = 40000.0 # cm^2
    nH = 0.02 # atoms/cm^2
    make_cosmological_sources_file(simput_prefix, phlist_prefix, exp_time, fov, 
                                   sky_center, nH=nH, area=area):

However, these same arguments accept values with unit information, either in the
form of ``(value, unit)`` tuples, :class:`~astropy.units.Quantity`, or
:class:`~yt.units.yt_array.YTQuantity` objects:

.. code-block:: python

    import soxs
    from astropy.units import Quantity
    simput_prefix = "cosmo"
    phlist_prefix = "cosmo"
    sky_center = [30.0, 45.0]
    exp_time = (500.0, "ks")
    fov = Quantity(0.666667, "deg")
    area = (4.0, "m**2") 
    nH = Quantity(2.0e20, "cm**-2") 
    make_cosmological_sources_file(simput_prefix, phlist_prefix, exp_time, fov, 
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