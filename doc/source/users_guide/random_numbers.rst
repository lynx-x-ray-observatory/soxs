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

    t_exp = 50000.0
    area = 30000.0
    prng = 24
    e = spec.generate_energies(t_exp, area, prng=prng)

Check the :ref:`api` to see which functions or methods allow for the input of 
random seeds. 