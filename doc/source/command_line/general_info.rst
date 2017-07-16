.. _general-info-cmd:

General Information Regarding the Command-Line Interface to SOXS
================================================================

Special Argument Handling for Quantities with Units
---------------------------------------------------

Many arguments in the command line scripts which have units can 
take a special format which allows one to specify that particular
quantity in the units desired by the user. For example, the 
:ref:`cmd-make-point-source` script has the arguments ``exp_time``
and ``area``, which assume the default units of seconds and :math:`\rm{cm^2}`,
respectively, if one supplies floating-point numbers:

.. code-block:: bash

    [~]$ make_point_source pt_src src1 20.0 -32.0 pt_src_spectrum.dat 100000.0 --overwrite --area=30000.0

but can take other units, like ks and :math:`\rm{m^2}`, in this format:

.. code-block:: bash

    [~]$ make_point_source pt_src src1 20.0 -32.0 pt_src_spectrum.dat 100.0,ks --overwrite --area=3.0,m**2

The following arguments used in the command line scripts accept values with a 
unit specification:

* ``exp_time``: Exposure time, used in many scripts, default units of seconds
* ``area``: Collecting area, used in many scripts, default units of :math:`\rm{cm}^2`
* ``emin``: Minimum energy, used in many scripts, default units of keV
* ``emax``: Minimum energy, used in many scripts, default units of keV
* ``nH``: Foreground galactic absorption column, used in many scripts, 
  default units of :math:`10^{22} \rm{atoms/cm^2}`
* ``kT``: Temperature, used in :ref:`cmd-make-thermal-spectrum`, default
  units of keV
* ``velocity``: Velocity broadening parameter, used in :ref:`cmd-make-thermal-spectrum`, 
  default units of km/s

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