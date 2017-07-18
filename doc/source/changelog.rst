.. _changelog:

ChangeLog
=========

Version 1.3.0
-------------

This is a release with important new features and some bugfixes.

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