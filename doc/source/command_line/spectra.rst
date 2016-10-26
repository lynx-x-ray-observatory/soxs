.. _cmd-spectra:

Command Line Scripts for Spectra
================================

These are scripts that create spectra for use by other modules in SOXS.

``make_thermal_spectrum``
-------------------------

.. code-block:: text

    usage: make_thermal_spectrum [-h] [--velocity VELOCITY] [--emin EMIN]
                                 [--emax EMAX] [--nbins NBINS]
                                 [--apec_vers APEC_VERS] [--no_broadening]
                                 [--absorb] [--nh NH] [--clobber]
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
                            The version of the AtomDB tables to use. Default:
                            3.0.3
      --no_broadening       Set this to turn off thermal and velocity broadening.
      --absorb              Whether or not to apply foreground galactic
                            absorption.
      --nh NH               The hydrogen column in units of 10**22 atoms/cm**2.
                            Default: 0.02
      --clobber             Whether or not to clobber an existing file with the
                            same name.
                            
Examples
++++++++

``make_powerlaw_spectrum``
--------------------------

.. code-block:: text

    usage: make_powerlaw_spectrum [-h] [--emin EMIN] [--emax EMAX] [--nbins NBINS]
                                  [--absorb] [--nh NH] [--clobber]
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
      --absorb       Whether or not to apply foreground galactic absorption.
      --nh NH        The hydrogen column in units of 10**22 atoms/cm**2. Default:
                     0.02
      --clobber      Whether or not to clobber an existing file with the same
                     name.
                 
Examples
++++++++