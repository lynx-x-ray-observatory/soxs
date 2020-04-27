.. _changelog:

ChangeLog
=========

Version 2.3.0
-------------

This version of SOXS fixes a few bugs, updates response files, and adds a speedup
for certain types of RMF convolutions.

* Fixed an issue in cosmological sources generation where comoving units were
  assumed when they should have been proper.
* Fixed an issue where the ``make_thermal_spectrum`` command-line script had the
  wrong default version of APEC specified. 
* *Chandra* imaging responses for ACIS-I and ACIS-S have been updated to Cycle 22. 
* The *XRISM*/Resolve ARF has been updated. 
* RMF convolutions are now faster in most situations.
* The field of view for the *Lynx*/HDXI has been increased from 20 arcmin to 22
  arcmin. 


Version 2.2.0
-------------

This version of SOXS focuses on new instrument modes and response files, as well
as containing bugfixes and improvements. This version supports Python 2.7, 3.5,
3.6, and 3.7.

* Fixed an issue where an invalid APEC version was being found when the user 
  did not specify a path to the AtomDB tables. Thanks to David Turner for this
  bugfix. 
* The *Lynx* microcalorimeter, now named "LXM", has been split into three 
  subarrays, currently corresponding to three different instrument modes. 
* The only *Lynx* mirror configuration currently available is the 
  :math:`d = 3~m, f = 10~m` system. All other confiugrations have been removed
  in this version of SOXS.
* A new naming scheme has been adopted for many instruments for clarity, but
  the old names will be accepted with a warning. 
* The *Chandra* Cycle 19 responses have been replaced by their Cycle 20 
  versions.
* The *Athena* response files have been updated to their latest versions. 
* The *Hitomi* SXS instrument mode has been replaced by the *XRISM* Resolve
  instrument mode, and the response files have been updated accordingly. 

For more information on the new instrument configurations, consult the
:ref:`instrument` section of the User's Guide. 

Version 2.1.0
-------------

This release of SOXS provides new features, bugfixes, optimizations, and other
improvements.

* The 2.1.x series of SOXS will be the last to support Python 2.7.
* Support for non-equilibrium ionization plasma emission using AtomDB has been
  added to SOXS. see :ref:`nei` for more details.
* The default AtomDB/APEC version provided with SOXS is now v3.0.9.
* Generating spectra without imaging using (see :ref:`simulate-spectrum`) is now
  faster, especially for high-resolution instruments such as microcalorimeters 
  and gratings.
* The default abundance table used when generating thermal spectra can now be set in
  the SOXS configuration file. See :ref:`solar-abund-tables` for more information.
* Absorption lines can now be added to spectra. See :ref:`absorb_lines` for more
  information.
* A new function for generating a simple imaging instrument based on an existing
  instrument specification has been added. See :ref:`simple-instruments` for more
  information. 
* A bug that prevented the multiplication of a
  :class:`~soxs.background.spectra.BackgroundSpectrum` object by a constant has
  been fixed.
* New convenience methods for generating :class:`~soxs.instrument.AuxiliaryResponseFile`
  and :class:`~soxs.instrument.RedistributionMatrixFile` objects directly from
  existing instrument specification names has been added.
* A new keyword argument, ``plot_counts``, has been added to the
  :func:`~soxs.events.plot_spectrum` function which allows the counts instead of
  the count rate to be plotted.
* The response files and instrumental background for the 
  `AXIS <http://axis.astro.umd.edu>`_ mission have been updated to their latest 
  versions.

Version 2.0.0
-------------

This is a major new release with a number of important new features and some bugfixes.

Most Important New Features and Changes
+++++++++++++++++++++++++++++++++++++++

* Beginning with this version and going forward, response files will not be included
  when SOXS is installed, primarily due to space considerations. Response files should
  be downloaded from the :ref:`responses` page either separately or as a whole.
  Instrument simulation can be performed with the response files located in the current
  working directory or in the default ``response_path`` specified in the SOXS 
  configuration file. See :ref:`config`, :ref:`response-path`, and :ref:`cmd-response-path`
  for more details.
* A configuration file can now be used with SOXS. See :ref:`config` for more details.
* The ability to simulate gratings spectra with :func:`~soxs.instrument.simulate_spectrum`
  and the ``simulate_spectrum`` command-line tool has been added. See :ref:`gratings` for 
  more information and :ref:`custom-non-imaging` for instructions on how to make a custom
  gratings instrument specification. Special thanks to `Lia Corrales <http://www.liacorrales.com/>`_
  for useful discussions and advice with respect to gratings spectra. 
* The :class:`~soxs.simput.SimputCatalog` and :class:`~soxs.simput.PhotonList` classes
  have been added for improved SIMPUT catalog handling, which greatly simplifies the 
  simulation of sources. See :ref:`simput` for more information. 
* A bug that prevented backgrounds from being added from a file properly to simulations
  with a non-zero roll angle has been fixed. 

Changes to Simulation of Spectra
++++++++++++++++++++++++++++++++

* A number of class methods for :class:`~soxs.spectra.Spectrum` and their associated
  command-line scripts now have ``emin``, ``emax``, and ``nbins`` as required arguments.
  Previously these were optional arguments. More information can be found at :ref:`spectra`
  and :ref:`cmd-spectra`. These are backwards-incompatible changes.
* The interpolating spline which allowed :class:`~soxs.spectra.Spectrum` objects to
  be called with an energy argument to get the values of the spectrum for arbitrary
  energies was not being regenerated if the spectrum was changed, say by foreground
  absorption. This has been fixed.
* The ability to apply intrinsic foreground absorption to a :class:`~soxs.spectra.Spectrum`
  has been added by adding an optional ``redshift`` argument to 
  :meth:`~soxs.spectra.Spectrum.apply_foreground_absorption`. 
* A method to easily plot :class:`~soxs.spectra.Spectrum` objects, 
  :meth:`~soxs.spectra.Spectrum.plot`, has been added. See :ref:`spectra-plots` for details.
* For APEC spectra created using :class:`~soxs.spectra.ApecGenerator`, it is now possible to
  use Solar abundance tables other than the implicitly assumed Anders & Grevesse 1989. See
  :ref:`solar-abund-tables` and :ref:`cmd-spectra` for details.
* The accuracy of the ``TBabs`` absorption model interpolation in SOXS has been improved.
* A method to add individual Gaussian-shaped lines to a :class:`~soxs.spectra.Spectrum`, 
  :meth:`~soxs.spectra.Spectrum.add_emission_line`, has been added. 
* The ability to write :class:`~soxs.spectra.Spectrum` objects to HDF5 files has
  been added via the :meth:`~soxs.spectra.Spectrum.write_h5_file` method. See
  :ref:`write-spectra` for details.

Changes to Instrument Simulation
++++++++++++++++++++++++++++++++

* :func:`~soxs.events.plot_spectrum` has been given more options. see :ref:`plot-spectrum`
  for details.
* A ``reblock`` optional argument has been added to :func:`~soxs.events.write_image` and
  :func:`~soxs.events.make_exposure_map` to allow the binning of images and exposure maps to
  be changed. See :ref:`event-tools` for details.
* Small improvements were made to reading parameters from RMFs, improving consistency
  and allowing more corner cases to be supported.
* If a ``COUNT_RATE`` column is not in a FITS table file containing a spectrum, the count 
  rate will be generated automatically in :func:`~soxs.events.plot_spectrum`.
* The ability to simulate background components has been added to 
  :func:`~soxs.instrument.simulate_spectrum`. See :ref:`simulate-spectrum` and
  :ref:`cmd-simulate-spectrum` for more details.
* The :meth:`~soxs.instrument.AuxiliaryResponseFile.plot` method of 
  :class:`~soxs.instrument.AuxiliaryResponseFile` now returns both a 
  :class:`~matplotlib.figure.Figure` and :class:`~matplotlib.axes.Axes` objects.

Changes to Instrument Specifications
++++++++++++++++++++++++++++++++++++

* An instrument specification for the *Lynx* gratings has been added to the instrument registry.
* Instrument specifications for *Chandra*/ACIS-S have been added to the instrument registry.
  Special thanks to Andrea Botteon for supplying the model for the ACIS-S particle background.
* Instrument specifications for *Chandra*/ACIS-S with the HETG have been added to the instrument
  registry. The instrument models correspond to the MEG and HEG :math:`\pm` first order.
* The *Chandra*/ACIS-I instrument specifications for Cycle 18 have been replaced with Cycle 19 
  specifications.
* When defining instrument specifications, it is now possible to specify a per-chip
  particle background model. See :ref:`custom-instruments` for more details.
* An instrument specification for the `AXIS <http://axis.astro.umd.edu>`_ mission
  concept has been added.

Version 1.3.0
-------------

This is a release with important new features and some bugfixes.

* SOXS now includes the ability to implement instruments with more than one chip
  with gaps in between, and chips which are not square in size. See :ref:`instrument`
  for more information.
* The *Chandra* ACIS-I instrument specifications have been changed so that they
  implement 4 chips in a 2x2 array, using the new SOXS chip functionality.
  The old specifications still exist in the instrument registry as ``"acisi_cy0_old"``
  and ``"acisi_cy18_old"``.
* The *Athena* WFI and X-IFU instrument specifications have been changed so that
  they more closely match the current models, using the new SOXS chip functionality.
  The old specifications still exist in the instrument registry as ``"athena_wfi_old"``
  and ``"athena_xifu_old"``.
* SOXS now has the ability to create exposure maps for SOXS simulations and use them
  when making images and radial profiles. See :ref:`event-tools` and :ref:`cmd-events` 
  for more information.
* Many arguments to functions and command line scripts which have units (such as 
  exposure time, field of view, area, temperature, etc.) now accept arguments with
  units. See :ref:`units` and :ref:`cmd-units` for more information.
* The "square" and "circle" dither pattern options have been replaced with a single
  option, a Lissajous pattern like that used by *Chandra*. This is a backwards-incompatible
  change.
* New methods have been added to create :class:`~soxs.spectra.ConvolvedSpectrum` objects
  and deconvolve them to :class:`~soxs.spectra.Spectrum` objects. See 
  :ref:`convolved-spectra` for more details.
* A method to extract a subset of a spectrum and create a new one, 
  :meth:`~soxs.spectra.Spectrum.new_spec_from_band`, has been added. 
* :class:`~soxs.spectra.Spectrum` objects are now "callable", taking an energy
  or an array of energies, at which the flux values will be interpolated.
* :class:`~soxs.spectra.ApecGenerator` objects can now generate spectra that 
  vary the elemental abundances separately. See :ref:`thermal-spectra` and 
  :ref:`cmd-make-thermal-spectrum` for more details.
* :class:`~soxs.spectra.ApecGenerator` objects can now generate spectra without 
  line emission. See :ref:`thermal-spectra` and :ref:`cmd-make-thermal-spectrum` 
  for more details.
* A bug that prevented one from adding new instrumental background spectra to the
  instrumental background spectrum registry has been fixed. 
* A bug that resulted in spectra being plotted with the incorrect energies in 
  :func:`~soxs.events.plot_spectrum` has been fixed.

Version 1.2.0
-------------

This is a release with three new features, a change in AtomDB version, and some
fixes to the documentation.

* An instrument specification for the *Hitomi*/SXS has been added. Thanks to
  Eric Miller of MIT for generating the response files.
* There are now two options for absorption models, "wabs" and "tbabs". All tools
  which take a parameter for the Galactic hydrogen column ``nH`` now take an
  optional parameter which can be set to ``"wabs"`` or ``"tbabs"``. The default 
  is still ``"wabs"``.
* SOXS now bundles only one version of the AtomDB tables, v3.0.8. It is still
  possible to point to your own directory containing a different version. 
* The :meth:`~soxs.spectra.Spectrum.from_file` method now accepts HDF5 files as
  input. 
* Various minor corrections to the documentation were made.

Version 1.1.1
-------------

This is a release with a single minor feature addition, which allows the foreground
galactic absorption parameter ``nH`` to be supplied to 
:func:`~soxs.instrument.make_background_file`, which is applied to the point-source
background.

Version 1.1.0
-------------

This is an important release that contains new features and bugfixes.

* The ability to provide an ASCII table of point source properties to re-use
  the same distribution of point sources has been added to 
  :func:`~soxs.background.point_sources.make_point_sources_file` and 
  :func:`~soxs.instrument.make_background_file`. 
* A new function, :func:`~soxs.background.point_sources.make_point_source_list`, has been
  added to provide a way to generate an ASCII table of point source properties
  for input into making background files and point source catalogs without
  having to create the events.
* For the point-source background, the photon spectral index for the galaxies is
  now :math:`\alpha = 2`, and the photon spectral index for the AGN is drawn
  from a fit to Figure 13a from 
  `Hickox & Markevitch 2006 <http://adsabs.harvard.edu/abs/2006ApJ...645...95H>`_.
* The *Athena* instrument models have been updated to more accurately reflect
  the current design parameters.
* A bug that prevented one from using an instrument model that did not have
  an instrumental background has been fixed.
* An experimental feature to turn off uniform randomization of events within
  pixels has been added.
* Dithering now occurs in detector coordinates instead of sky coordinates.

Version 1.0.1
-------------

This is a bugfix release to fix the fact that the ``soxs.background`` submodule
was not being imported properly. 

Version 1.0.0
-------------

This version is a major new release with a complete revamp of the way that
SOXS handles backgrounds, as well as a number of other new features and 
bugfixes.

* Backgrounds will now either be added when running the instrument simulator
  or can be created separately for a particular instrument, saved to an event
  file, and then used for multiple observations. This enables one to avoid having 
  to create a background for every observation, which can be prohibitive for 
  long exposures. 
* Added a point-source component to the astrophysical background. 
* The background keyword arguments for :func:`~soxs.instrument.instrument_simulator`
  are now ``instr_bkgnd``, ``foreground``, and ``ptsrc_bkgnd``. ``astro_bkgnd``
  has been removed. This is a backwards-incompatible change. 
* Added the capability to create a source composed of cosmological halos drawn
  from a cosmological simulation. 
* Instrument specifications for *Chandra*/ACIS-I have been added, with responses
  from Cycle 0 and Cycle 18. 
* SOXS now has the new dependencies of `h5py <http://www.h5py.org>`_ and 
  `SciPy <http://www.scipy.org>`_, as well as `AstroPy <http://www.astropy.org>`_ 
  version 1.3. 
* Added the ability to specify a name for a source in a SIMPUT catalog when
  writing a photon list file.
* Test coverage has been improved, especially for backgrounds. 
* Tests are now performed on Python versions 2.7, 3.5, and 3.6.
* In the Python interface, integers may now be provided for random seeds as
  arguments to functions. 
* An argument to provide a random seed to generate a consistent set of random
  numbers has been added to all of the command line scripts which make use of
  random numbers. 
* Fixed a bug in determining the detector and chip coordinates of events when
  creating an event file. 
* The ``clobber`` argument for overwriting files has been replaced by 
  ``overwrite``. This is a backwards-incompatible change.

Version 0.5.1
-------------

This version is a bugfix release. 

* Fixed a big when writing FITS table files when AstroPy 1.3 is installed. 

Version 0.5.0
-------------

This version contains new features and bugfixes.

* The PSF can now be set to ``None`` (or ``null`` in JSON files) in an 
  instrument specification for no PSF scattering of events.
* The particle background can be set to ``None`` (or ``null`` in JSON files) in
  an instrument specification for no particle background.
* A faster progress bar, `tqdm <https://github.com/tqdm/tqdm>`_, is now in use 
  in SOXS.
* Fixed a minor bug in the interpolation of APEC tables for thermal spectra. The
  difference in the generated spectra is small, at around the fifth decimal 
  place.
* Added a constant spectrum generator: :meth:`~soxs.spectra.Spectrum.from_constant`.
* Added ellipticity and angle parameters to :class:`~soxs.spatial.RadialFunctionModel` 
  objects to create models with ellipticity.
* Added flat-field coordinates to :class:`~soxs.spatial.SpatialModel` objects.
* Made public subclass of :class:`~soxs.spectra.Spectrum` objects, 
  :class:`~soxs.spectra.ConvolvedSpectrum`, which is a :class:`~soxs.spectra.Spectrum` 
  convolved with an ARF.
* Small internal changes designed to provide a more seamless interface to 
  `pyXSIM <http://hea-www.cfa.harvard.edu/~jzuhone/pyxsim>`_.
* Three new tools have been included to produce derivative products from event 
  files:

  * :func:`~soxs.events.write_image`: Bins events into an image and writes it to
    a FITS file.
  * :func:`~soxs.events.write_spectrum`: Bins events into a spectrum and writes it
    to a FITS file.
  * :func:`~soxs.events.write_radial_profiles`: Bins events into a radial 
    profile and writes it to a FITS file.

Version 0.4.0
-------------

This version contains new features and bugfixes. Some changes are not 
backwards-compatible. 

* SOXS has been re-branded as "Simulating Observations of X-ray Sources".
* Instrument specifications for the *Athena* WFI and X-IFU have been added to 
  the instrument registry.
* A test suite infrastructure has been added to SOXS, which runs automatically 
  on GitHub when changes are made to the source code. 
* Simulating backgrounds without an input source is now possible by providing 
  ``None`` to :func:`~soxs.instrument.instrument_simulator` or ``"None"`` to the
  ``instrument_simulator`` command line script (see :ref:`cmd-instrument`).
* The default astrophysical background in SOXS was not identical to the 
  advertised input spectrum, which has been fixed.
* The options for dealing with background have been restricted. Backgrounds can 
  now only be turned on and off. The keyword arguments to 
  :func:`~soxs.instrument.instrument_simulator` for dealing with background have
  been correspondingly modified (see :ref:`instrument` and 
  :ref:`cmd-instrument`). This is a backwards-incompatible change.
* The default version of APEC in :class:`~soxs.spectra.ApecGenerator` is now 
  version 2.0.2, to match XSPEC. 
* A new option has been added to the instrument specification to turn dithering 
  on and off by default for a given instrument. Please change instrument 
  specification JSON files accordingly.
* Instead of the plate scale, the instrument field of view is specified in the 
  instrument specification, and the plate scale is calculated from this and the 
  number of pixels. Please change instrument specification JSON files 
  accordingly.

Version 0.3.1
-------------

This is a bugfix release.

* The RMF for the HDXI was updated so that the binning between it and the HDXI 
  ARFs is consistent.
* Various small edits to the documentation were made.

Version 0.3.0
-------------

This version contains new features and bugfixes.

* An *Athena*-like microcalorimeter background is now the default particle 
  background for all microcalorimeter models.
* All instrumental backgrounds now have a dependence on the focal length. The 
  focal length is now an element of the instrument specification. 
* The names of the instruments in the instrument registry were made consistent 
  with their associated keys.
* A convenience function, :meth:`~soxs.spectra.Spectrum.get_flux_in_band`, has 
  been added. 
* A new method of generating a spectrum from an XSPEC script, 
  :meth:`~soxs.spectra.Spectrum.from_xspec_script`, has been added.
* The :meth:`~soxs.spectra.Spectrum.from_xspec` method has been renamed to 
  :meth:`~soxs.spectra.Spectrum.from_xspec_model`. 
* Removed unnecessary commas between coordinate values from the examples in 
  :ref:`cmd-spatial`. 
* Added a new capability to create a SIMPUT file from an ASCII table of RA, Dec,
  and energy, in the ``make_phlist_from_ascii`` command-line script.
* Added a new class for creating rectangle/line-shaped sources, 
  :class:`~soxs.spatial.RectangleModel`, and a corresponding command-line 
  script, ``make_rectangle_source``. 
* The signature of ``write_photon_list`` has changed to accept a ``flux`` 
  argument instead of exposure time and area.

Version 0.2.1
-------------

This is a bugfix release.

* The supporting files (ARFs, RMFs, spectral files, etc.) were not being bundled
  properly in previous versions. 

Version 0.2.0
-------------

This version contains new features.

* New ARFs corresponding to various configurations of the mirrors have been 
  added and the old ARFs have been removed (November 1st, 2016).
* Documentation now includes references to ways of getting help and the license.

Version 0.1.1
-------------

This is solely a bugfix release.

* Fixed a bug where the dither did not have the correct width.
* Fixed a bug for cases with no dithering.
* Various minor improvements to the documentation