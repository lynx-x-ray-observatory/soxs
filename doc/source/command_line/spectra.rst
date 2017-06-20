.. _cmd-spectra:

Command Line Scripts for Spectra
================================

These are scripts that create ASCII tables of spectra for use by other 
modules in SOXS. For details on what's going on under the hood, see :ref:`spectra`.

``make_thermal_spectrum``
-------------------------

This script creates an ASCII table of an optionally absorbed thermal spectrum 
using the AtomDB tables.

.. code-block:: text

    usage: make_thermal_spectrum [-h] [--velocity VELOCITY] [--emin EMIN]
                                 [--emax EMAX] [--nbins NBINS]
                                 [--apec_vers APEC_VERS] [--absorb ABSORB] 
                                 [--nh NH] [--overwrite] [--broadening | --no_broadening]

                                 kT abund redshift norm specfile
    
    Create a thermal spectrum and write it to a file.
    
    positional arguments:
      kT                    The temperature in keV.
      abund                 The metal abundance in solar units.
      redshift              The redshift of the source.
      norm                  The normalization of the model, in the standard Xspec
                            units of 1.0e-14*EM/(4*pi*(1+z)**2*D_A**2).
      specfile              The filename to write the spectrum to.
    
    optional arguments:
      -h, --help            show this help message and exit
      --velocity VELOCITY   The velocity broadening parameter, in units of km/s.
                            Default: 0.0
      --emin EMIN           The minimum energy in keV. Default: 0.01
      --emax EMAX           The maximum energy in keV. Default: 50.0
      --nbins NBINS         The number of bins in the spectrum. Default: 10000
      --apec_vers APEC_VERS
                            The version of the AtomDB tables to use. Default: 3.0.3
      --absorb ABSORB       Apply foreground Galactic absorption. Default is
                            "wabs".
      --nh NH               The hydrogen column in units of 10**22 atoms/cm**2.
                            Default: 0.02
      --overwrite           Overwrite an existing file with the same name.
      --broadening          Turn thermal and velocity broadening on. On by default.
      --no_broadening       Turn thermal and velocity broadening off. On by default.

Examples
++++++++

Make a basic spectrum for a thermal plasma. 

.. code-block:: bash

    [~]$ make_thermal_spectrum 6.0 0.3 0.05 1.0e-4 my_thermal_spectrum.dat --overwrite

The same spectrum, but with velocity broadening.

.. code-block:: bash

    [~]$ make_thermal_spectrum 6.0 0.3 0.05 1.0e-4 my_thermal_spectrum.dat --velocity=200.0 --overwrite

The same spectrum, but with velocity and thermal broadening turned off.

.. code-block:: bash

    [~]$ make_thermal_spectrum 6.0 0.3 0.05 1.0e-4 my_thermal_spectrum.dat --no_broadening --overwrite

The same spectrum, but with foreground galactic absorption using the "wabs" model
with :math:`N_H = 0.04~10^{22}~\rm{atoms~cm^{-2}}`.

.. code-block:: bash

    [~]$ make_thermal_spectrum 6.0 0.3 0.05 1.0e-4 my_thermal_spectrum.dat --absorb "wabs" --nh 0.04 --overwrite

The same spectrum, but with different binning.

.. code-block:: bash

    [~]$ make_thermal_spectrum 6.0 0.3 0.05 1.0e-4 my_thermal_spectrum.dat --emin=0.1 --emax=10.0 --nbins=20000 --overwrite

The same spectrum, but with a different APEC version.

.. code-block:: bash

    [~]$ make_thermal_spectrum 6.0 0.3 0.05 1.0e-4 my_thermal_spectrum.dat --apec_vers=2.0.2 --overwrite

``make_powerlaw_spectrum``
--------------------------

This script creates an ASCII table of an optionally absorbed power-law spectrum. This spectrum has the
form:

.. math::

    F_E = K\left[\frac{E(1+z)}{{\rm 1~keV}}\right]^{-\alpha}

.. code-block:: text

    usage: make_powerlaw_spectrum [-h] [--emin EMIN] [--emax EMAX] [--nbins NBINS]
                                  [--absorb ABSORB] [--nh NH] [--overwrite]
                                  photon_index redshift norm specfile
    
    Create a power-law spectrum and write it to a file.
    
    positional arguments:
      photon_index   The spectral index of the power law.
      redshift       The redshift of the source.
      norm           The normalization of the source in units of
                     photons/s/cm**2/keV at 1 keV in the source frame.
      specfile       The filename to write the spectrum to.
    
    optional arguments:
      -h, --help     show this help message and exit
      --emin EMIN    The minimum energy in keV. Default: 0.01
      --emax EMAX    The maximum energy in keV. Default: 50.0
      --nbins NBINS  The number of bins in the spectrum. Default: 10000
      --absorb ABSORB       Apply foreground Galactic absorption. Default is
                            "wabs".
      --nh NH        The hydrogen column in units of 10**22 atoms/cm**2. Default:
                     0.02
      --overwrite    Overwrite an existing file with the same name.
                 
Examples
++++++++

Make a basic power-law spectrum. 

.. code-block:: bash

    [~]$ make_powerlaw_spectrum 1.1 0.05 1.0e-4 my_powerlaw_spectrum.dat --overwrite

The same spectrum, but with foreground galactic absorption using the "tbabs" model
with :math:`N_H = 0.04~10^{22}~\rm{atoms~cm^{-2}}`.

.. code-block:: bash

    [~]$ make_powerlaw_spectrum 1.1 0.05 1.0e-4 my_powerlaw_spectrum.dat --absorb "tbabs" --nh 0.04 --overwrite

The same spectrum, but with different binning.

.. code-block:: bash

    [~]$ make_powerlaw_spectrum 1.1 0.05 1.0e-4 my_powerlaw_spectrum.dat --emin=0.1 --emax=10.0 --nbins=20000 --overwrite
