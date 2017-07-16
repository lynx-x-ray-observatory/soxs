.. _random-numbers-cmd:

Random Number Generation
------------------------

Many routines in SOXS require generating random numbers for energies, sky
positions, spectral channels, etc. By default, for every SOXS run this will
be a different set of random numbers. It is often the case, however, that one
wants to use a consistent, repeatable set of random numbers to reproduce results
exactly. For this, many of the command-line scripts in SOXS take a 
``random_seed`` optional argument, which has a default of ``None``, but if set 
to an integer will use this value as a random seed. 

For example, to use a consistent random seed in generating positions and
energies for an annulus source, one would set the random seed
like this:

.. code-block:: bash

    make_annulus_source my_cat annulus 30.0 45.0 10.0 30.0 thermal_spec.dat 100.0,ks --random_seed=24

Check the documentation for the various command line scripts to see which
functions have the ``random_seed`` argument. 