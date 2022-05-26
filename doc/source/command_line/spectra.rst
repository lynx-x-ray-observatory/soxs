.. _cmd-spectra:

Command Line Scripts for Spectra
================================

These are scripts that create ASCII tables of spectra for use by other 
modules in SOXS. For details on what's going on under the hood, see :ref:`spectra`.

.. _cmd-make-thermal-spectrum:

``make_thermal_spectrum``
-------------------------

This script creates an ASCII table of an optionally absorbed thermal spectrum 
using the AtomDB tables.

.. code-block:: text

    usage: make_thermal_spectrum [-h] [--velocity VELOCITY]
                                 [--apec_vers APEC_VERS] [--binscale BINSCALE]
                                 [--absorb_model ABSORB_MODEL] [--nh NH]
                                 [--overwrite] [--nolines]
                                 [--abund_table ABUND_TABLE]
                                 [--broadening | --no_broadening]
                                 kT abund redshift norm specfile emin emax nbins
    
    Create a thermal spectrum and write it to a file. The abundances of individual
    elements can be set by supplying optional arguments in the form of --O=0.5,
    --Mg=0.6, etc.
    
    positional arguments:
      kT                    The temperature in keV.
      abund                 The metal abundance in solar units.
      redshift              The redshift of the source.
      norm                  The normalization of the model, in the standard Xspec
                            units of 1.0e-14*EM/(4*pi*(1+z)**2*D_A**2).
      specfile              The filename to write the spectrum to.
      emin                  The minimum energy in keV.
      emax                  The maximum energy in keV.
      nbins                 The number of bins in the spectrum.
    
    optional arguments:
      -h, --help            show this help message and exit
      --velocity VELOCITY   The velocity broadening parameter, in units of km/s.
                            Default: 0.0
      --apec_vers APEC_VERS
                            The version of the AtomDB tables to use. Default is
                            to use the version currently included with this
                            version of SOXS.
      --binscale BINSCALE   The scale of the energy binning: "linear" or "log".
                            Default: "linear"
      --absorb_model ABSORB_MODEL
                            Model for applying foreground Galactic absorption.
      --nh NH               The hydrogen column in units of 10**22 atoms/cm**2.
      --overwrite           Overwrite an existing file with the same name.
      --nolines             Make a spectrum without lines.
      --abund_table ABUND_TABLE
                            The abundance table to be used for solar abundances.
                            Either a string corresponding to a built-in table or
                            an ASCII fiele containing a column of 30 floats
                            corresponding to the abundances of each element
                            relative to the abundance of H. Default is set in the 
                            SOXS configuration file, the default for which is "angr".
      --broadening          Turn thermal and velocity broadening on. On by
                            default.
      --no_broadening       Turn thermal and velocity broadening off. On by
                            default.

Examples
++++++++

Make a basic spectrum for a thermal plasma. 

.. code-block:: bash

    [~]$ make_thermal_spectrum 6.0 0.3 0.05 1.0e-4 my_thermal_spectrum.dat 0.1 10.0 10000 --overwrite

The same spectrum, but with velocity broadening.

.. code-block:: bash

    [~]$ make_thermal_spectrum 6.0 0.3 0.05 1.0e-4 my_thermal_spectrum.dat 0.1 10.0 10000 --velocity=200.0 --overwrite

The same spectrum, but with velocity and thermal broadening turned off.

.. code-block:: bash

    [~]$ make_thermal_spectrum 6.0 0.3 0.05 1.0e-4 my_thermal_spectrum.dat 0.1 10.0 10000 --no_broadening --overwrite

The same spectrum, but with foreground galactic absorption using the "wabs" model
with :math:`N_H = 0.04~\times~10^{22}~\rm{atoms~cm^{-2}}`.

.. code-block:: bash

    [~]$ make_thermal_spectrum 6.0 0.3 0.05 1.0e-4 my_thermal_spectrum.dat 0.1 10.0 10000 --absorb_model="wabs" --nh 0.04 --overwrite

The same spectrum, but with a different APEC version.

.. code-block:: bash

    [~]$ make_thermal_spectrum 6.0 0.3 0.05 1.0e-4 my_thermal_spectrum.dat 0.1 10.0 10000 --apec_vers=2.0.2 --overwrite

The same spectrum, but without emission lines. 

.. code-block:: bash

    [~]$ make_thermal_spectrum 6.0 0.3 0.05 1.0e-4 my_thermal_spectrum.dat 0.1 10.0 10000 --nolines --overwrite

The same spectrum, but setting the abundances of elements oxygen and calcium separately.

.. code-block:: bash

    [~]$ make_thermal_spectrum 6.0 0.3 0.05 1.0e-4 my_thermal_spectrum.dat 0.1 10.0 10000 --O=0.5 --Ca=0.7 --overwrite

The same spectrum, but using Asplund abundances instead of Anders & Grevesse.

.. code-block:: bash

    [~]$ make_thermal_spectrum 6.0 0.3 0.05 1.0e-4 my_thermal_spectrum.dat 0.1 10.0 10000 --abund_table=aspl --overwrite

The same spectrum, but with log-spaced binning.

.. code-block:: bash

    [~]$ make_thermal_spectrum 6.0 0.3 0.05 1.0e-4 my_thermal_spectrum.dat 0.1 10.0 10000 --binscale=log --overwrite

The same spectrum, but using abundances drawn from an ASCII table file instead of Anders & Grevesse.

.. code-block:: bash

    [~]$ make_thermal_spectrum 6.0 0.3 0.05 1.0e-4 my_thermal_spectrum.dat 0.1 10.0 10000 --abund_table=my_abund.dat --overwrite

``make_powerlaw_spectrum``
--------------------------

This script creates an ASCII table of an optionally absorbed power-law spectrum. This spectrum has the
form:

.. math::

    F_E = K\left[\frac{E(1+z)}{{\rm 1~keV}}\right]^{-\alpha}

.. code-block:: text

    usage: make_powerlaw_spectrum [-h] [--binscale BINSCALE]
                                  [--absorb_model ABSORB_MODEL] [--nh NH]
                                  [--overwrite]
                                  photon_index redshift norm specfile emin emax
                                  nbins
    
    Create a power-law spectrum and write it to a file.
    
    positional arguments:
      photon_index          The spectral index of the power law.
      redshift              The redshift of the source.
      norm                  The normalization of the source in units of
                            photons/s/cm**2/keV at 1 keV in the source frame.
      specfile              The filename to write the spectrum to.
      emin                  The minimum energy in keV.
      emax                  The maximum energy in keV.
      nbins                 The number of bins in the spectrum.
    
    optional arguments:
      -h, --help            show this help message and exit
      --binscale BINSCALE   The scale of the energy binning: "linear" or "log".
                            Default: "linear"
      --absorb_model ABSORB_MODEL
                            Model for applying foreground Galactic absorption.
      --nh NH               The hydrogen column in units of 10**22 atoms/cm**2.
      --overwrite           Overwrite an existing file with the same name.

Examples
++++++++

Make a basic power-law spectrum. 

.. code-block:: bash

    [~]$ make_powerlaw_spectrum 1.1 0.05 1.0e-4 my_powerlaw_spectrum.dat 0.1 10.0 10000 --overwrite

The same spectrum, but with foreground galactic absorption using the "tbabs" model
with :math:`N_H = 0.04~10^{22}~\rm{atoms~cm^{-2}}`.

.. code-block:: bash

    [~]$ make_powerlaw_spectrum 1.1 0.05 1.0e-4 my_powerlaw_spectrum.dat 0.1 10.0 10000 --absorb_model="tbabs" --nh 0.04 --overwrite

The same spectrum, but with log-spaced binning.

.. code-block:: bash

    [~]$ make_powerlaw_spectrum 1.1 0.05 1.0e-4 my_powerlaw_spectrum.dat 0.1 10.0 10000 --binscale=log --overwrite
