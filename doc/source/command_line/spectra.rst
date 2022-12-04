.. _cmd-spectra:

Command Line Scripts for Spectra
================================

These are scripts that create ASCII tables of spectra for use by other
modules in SOXS. For details on what's going on under the hood, see :ref:`thermal-spectra`.

.. _cmd-make-cie-spectrum:

``make_cie_spectrum``
---------------------

This script creates an ASCII table of an optionally absorbed thermal spectrum
using a CIE model. Options are APEC (the default), SPEX, MeKaL, and Cloudy. For more
details on these options, see :ref:`thermal-spectra`.

.. code-block:: text

    usage: make_cie_spectrum [-h] [--velocity VELOCITY] [--model_vers MODEL_VERS]
                             [--binscale BINSCALE] [--absorb_model ABSORB_MODEL]
                             [--nH_abs NH_ABS] [--overwrite] [--nolines]
                             [--abund_table ABUND_TABLE] [--model MODEL]
                             [--broadening | --no_broadening]
                             kT abund redshift norm specfile emin emax nbins

    Create a thermal CIE spectrum and write it to a file. The abundances of individual
    elements can be set by supplying optional arguments in the form of --O=0.5, --Mg=0.6,
    etc.

    positional arguments:
      kT                    The temperature in keV.
      abund                 The metal abundance in solar units.
      redshift              The redshift of the source.
      norm                  The normalization of the model, in the standard Xspec units
                            of 1.0e-14*EM/(4*pi*(1+z)**2*D_A**2).
      specfile              The filename to write the spectrum to.
      emin                  The minimum energy in keV.
      emax                  The maximum energy in keV.
      nbins                 The number of bins in the spectrum.

    options:
      -h, --help            show this help message and exit
      --velocity VELOCITY   The velocity broadening parameter, in units of km/s. Default:
                            0.0 Only available for 'apec' and 'spex' models.
      --model_vers MODEL_VERS
                            The version of the CIE tables to use. Default is to use the
                            version currently included with this version of SOXS.
      --binscale BINSCALE   The scale of the energy binning: "linear" or "log". Default:
                            "linear"
      --absorb_model ABSORB_MODEL
                            Model for applying foreground Galactic absorption.
      --nH_abs NH_ABS       The hydrogen column in units of 10**22 atoms/cm**2.
      --overwrite           Overwrite an existing file with the same name.
      --nolines             Make a spectrum without lines. Only available for 'apec' and
                            'spex' models.
      --abund_table ABUND_TABLE
                            The abundance table to be used for solar abundances. Either a
                            string corresponding to a built-in table or an ASCII file
                            containing a column of 30 floats corresponding to the
                            abundances of each element relative to the abundance of H.
                            Default is set in the SOXS configuration file, the default
                            for which is 'angr'.
      --model MODEL         The CIE model to use when generating the spectrum, either
                            'apec', 'spex', 'mekal', or 'cloudy'. Default: 'apec'
      --broadening          Turn thermal and velocity broadening on. On by default. Only
                            available for 'apec' and 'spex' models.
      --no_broadening       Turn thermal and velocity broadening off. On by default. Only
                            available for 'apec' and 'spex' models.

Examples
++++++++

Make a basic spectrum for a thermal plasma.

.. code-block:: bash

    [~]$ make_cie_spectrum 6.0 0.3 0.05 1.0e-4 my_thermal_spectrum.dat 0.1 10.0 10000 --overwrite

The same spectrum, but with velocity broadening.

.. code-block:: bash

    [~]$ make_cie_spectrum 6.0 0.3 0.05 1.0e-4 my_thermal_spectrum.dat 0.1 10.0 10000 --velocity=200.0 --overwrite

The same spectrum, but with velocity and thermal broadening turned off.

.. code-block:: bash

    [~]$ make_cie_spectrum 6.0 0.3 0.05 1.0e-4 my_thermal_spectrum.dat 0.1 10.0 10000 --no_broadening --overwrite

The same spectrum, but with foreground galactic absorption using the "wabs" model
with :math:`N_H = 0.04~\times~10^{22}~\rm{atoms~cm^{-2}}`.

.. code-block:: bash

    [~]$ make_cie_spectrum 6.0 0.3 0.05 1.0e-4 my_thermal_spectrum.dat 0.1 10.0 10000 --absorb_model="wabs" --nH_abs 0.04 --overwrite

The same spectrum, but with a different APEC version.

.. code-block:: bash

    [~]$ make_cie_spectrum 6.0 0.3 0.05 1.0e-4 my_thermal_spectrum.dat 0.1 10.0 10000 --model_vers=2.0.2 --overwrite

The same spectrum, but without emission lines.

.. code-block:: bash

    [~]$ make_cie_spectrum 6.0 0.3 0.05 1.0e-4 my_thermal_spectrum.dat 0.1 10.0 10000 --nolines --overwrite

The same spectrum, but setting the abundances of elements oxygen and calcium separately.

.. code-block:: bash

    [~]$ make_cie_spectrum 6.0 0.3 0.05 1.0e-4 my_thermal_spectrum.dat 0.1 10.0 10000 --O=0.5 --Ca=0.7 --overwrite

The same spectrum, but using Asplund abundances instead of Anders & Grevesse.

.. code-block:: bash

    [~]$ make_cie_spectrum 6.0 0.3 0.05 1.0e-4 my_thermal_spectrum.dat 0.1 10.0 10000 --abund_table=aspl --overwrite

The same spectrum, but with log-spaced binning.

.. code-block:: bash

    [~]$ make_cie_spectrum 6.0 0.3 0.05 1.0e-4 my_thermal_spectrum.dat 0.1 10.0 10000 --binscale=log --overwrite

The same spectrum, but using abundances drawn from an ASCII table file instead of Anders & Grevesse.

.. code-block:: bash

    [~]$ make_cie_spectrum 6.0 0.3 0.05 1.0e-4 my_thermal_spectrum.dat 0.1 10.0 10000 --abund_table=my_abund.dat --overwrite

The same spectrum, but using the SPEX model instead of APEC.

.. code-block:: bash

    [~]$ make_cie_spectrum 6.0 0.3 0.05 1.0e-4 my_thermal_spectrum.dat 0.1 10.0 10000 --model=spex --overwrite

.. _cmd-make-igm-spectrum:

``make_igm_spectrum``
---------------------

This script creates an ASCII table of an optionally absorbed thermal spectrum
using the SOXS IGM model. For more details on what's going on under the hood,
see :ref:`igm-spectra`.

.. code-block:: text

    usage: make_igm_spectrum [-h] [--binscale BINSCALE] [--resonant_scattering]
                             [--cxb_factor CXB_FACTOR]
                             [--absorb_model ABSORB_MODEL] [--nH_abs NH_ABS]
                             [--overwrite]
                             kT nH abund redshift norm specfile emin emax nbins

    Create a thermal spectrum using the SOXS IGM model and write it to a file. The
    abundances of individual elements can be set by supplying optional arguments in the
    form of --O=0.5, --Mg=0.6, etc.

    positional arguments:
      kT                    The temperature in keV.
      nH                    The hydrogen number density in cm**-3.
      abund                 The metal abundance in solar units.
      redshift              The redshift of the source.
      norm                  The normalization of the model, in the standard Xspec units
                            of 1.0e-14*EM/(4*pi*(1+z)**2*D_A**2).
      specfile              The filename to write the spectrum to.
      emin                  The minimum energy in keV.
      emax                  The maximum energy in keV.
      nbins                 The number of bins in the spectrum.

    options:
      -h, --help            show this help message and exit
      --binscale BINSCALE   The scale of the energy binning: "linear" or "log". Default:
                            "linear"
      --resonant_scattering
                            Whether or not to include the effects of resonant scattering
                            from CXB photons. Default: False
      --cxb_factor CXB_FACTOR
                            The fraction of the CXB photons that are resonant scattered
                            to enhance the lines. Default: 0.5
      --absorb_model ABSORB_MODEL
                            Model for applying foreground Galactic absorption.
      --nH_abs NH_ABS       The hydrogen column in units of 10**22 atoms/cm**2.
      --overwrite           Overwrite an existing file with the same name.

Examples
++++++++

Make a basic IGM spectrum for a thermal plasma.

.. code-block:: bash

    [~]$ make_igm_spectrum 6.0 1.0e-3 0.3 0.05 1.0e-4 my_thermal_spectrum.dat 0.1 10.0 10000 --overwrite

The same spectrum, but with foreground galactic absorption using the "tbabs" model
with :math:`N_H = 0.04~\times~10^{22}~\rm{atoms~cm^{-2}}`.

.. code-block:: bash

    [~]$ make_igm_spectrum 6.0 1.0e-3 0.3 0.05 1.0e-4 my_thermal_spectrum.dat 0.1 10.0 10000 --absorb_model="tbbs" --nH_abs 0.04 --overwrite

The same spectrum, but with log-spaced binning.

.. code-block:: bash

    [~]$ make_igm_spectrum 6.0 1.0e-3 0.3 0.05 1.0e-4 my_thermal_spectrum.dat 0.1 10.0 10000 --binscale=log --overwrite

The same spectrum, but with resonant scattering.

.. code-block:: bash

    [~]$ make_igm_spectrum 6.0 1.0e-3 0.3 0.05 1.0e-4 my_thermal_spectrum.dat 0.1 10.0 10000 --resonant_scattering --overwrite

The same spectrum, but with resonant scattering and only 0.3 of the CXB scattered.

.. code-block:: bash

    [~]$ make_igm_spectrum 6.0 1.0e-3 0.3 0.05 1.0e-4 my_thermal_spectrum.dat 0.1 10.0 10000 --cxb_factor=0.3 --resonant_scattering --overwrite

The same spectrum, but with abundances of O, Ne, and Fe specified.

.. code-block:: bash

    [~]$ make_igm_spectrum 6.0 1.0e-3 0.3 0.05 1.0e-4 my_thermal_spectrum.dat 0.1 10.0 10000 --O=0.7 --Ne=0.6 --Fe=0.8 --overwrite

``make_powerlaw_spectrum``
--------------------------

This script creates an ASCII table of an optionally absorbed power-law spectrum. This spectrum has the
form:

.. math::

    F_E = K\left[\frac{E(1+z)}{{\rm 1~keV}}\right]^{-\alpha}

.. code-block:: text

    usage: make_powerlaw_spectrum [-h] [--binscale BINSCALE]
                                  [--absorb_model ABSORB_MODEL] [--nH_abs NH_ABS]
                                  [--abund_table ABUND_TABLE] [--overwrite]
                                  photon_index redshift norm specfile emin emax nbins

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

    options:
      -h, --help            show this help message and exit
      --binscale BINSCALE   The scale of the energy binning: "linear" or "log". Default:
                            "linear"
      --absorb_model ABSORB_MODEL
                            Model for applying foreground Galactic absorption.
      --nH_abs NH_ABS       The hydrogen column in units of 10**22 atoms/cm**2.
      --abund_table ABUND_TABLE
                            The abundance table to be used if the absorption model is
                            TBabs. Takes a string corresponding to a built-in table.
                            Default is set in the SOXS configuration file, the default
                            for which is 'angr'.
      --overwrite           Overwrite an existing file with the same name.

Examples
++++++++

Make a basic power-law spectrum.

.. code-block:: bash

    [~]$ make_powerlaw_spectrum 1.1 0.05 1.0e-4 my_powerlaw_spectrum.dat 0.1 10.0 10000 --overwrite

The same spectrum, but with foreground galactic absorption using the "tbabs" model
with :math:`N_H = 0.04~10^{22}~\rm{atoms~cm^{-2}}`.

.. code-block:: bash

    [~]$ make_powerlaw_spectrum 1.1 0.05 1.0e-4 my_powerlaw_spectrum.dat 0.1 10.0 10000 --absorb_model="tbabs" --nH_abs 0.04 --overwrite

The same spectrum, but with log-spaced binning.

.. code-block:: bash

    [~]$ make_powerlaw_spectrum 1.1 0.05 1.0e-4 my_powerlaw_spectrum.dat 0.1 10.0 10000 --binscale=log --overwrite

The same spectrum, but switching the abundance table used for the "tbabs" model.

.. code-block:: bash

    [~]$ make_powerlaw_spectrum 1.1 0.05 1.0e-4 my_powerlaw_spectrum.dat 0.1 10.0 10000 --absorb_model="tbabs" --nH_abs 0.04 --abund_table="feld" --overwrite
