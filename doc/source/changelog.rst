.. _changelog:

ChangeLog
=========

Version 0.3.0
-------------

This version contains new features and bugfixes.

* An *Athena*-like microcalorimeter background is now the default particle background for all microcalorimeter models.
* All instrumental backgrounds now have a dependence on the focal length. The focal length is now an element of the
  instrument specification. 
* The ``name``s of the instruments in the instrument registry were made consistent with their associated keys.
* A convenience function, :meth:`~soxs.spectra.Spectrum.get_flux_in_band`, has been added. 
* A new method of generating a spectrum from an XSPEC script, :meth:`~soxs.spectra.Spectrum.from_xspec_script`, has been added.
* The :meth:`~soxs.spectra.Spectrum.from_xspec` method has been renamed to :meth:`~soxs.spectra.Spectrum.from_xspec_model`. 
* Removed unnecessary commas between coordinate values from the examples in :ref:`cmd-spatial`. 
* Added a new capability to create a SIMPUT file from an ASCII table of RA, Dec, and energy, 
  in the ``make_phlist_from_ascii`` command-line script.
* Added a new class for creating rectangle/line-shaped sources, :class:`~soxs.spatial.RectangleModel`, and a corresponding
  command-line script, ``make_rectangle_source``. 
* The signature of ``write_photon_list`` has changed to accept a ``flux`` argument instead of exposure time and area.

Version 0.2.1
-------------

This is a bugfix release.

* The supporting files (ARFs, RMFs, spectral files, etc.) were not being bundled properly in previous versions. 

Version 0.2.0
-------------

This version contains new features.

* New ARFs corresponding to various configurations of the mirrors have been added and the old ARFs have been
  removed (November 1st, 2016).
* Documentation now includes references to ways of getting help and the license.

Version 0.1.1
-------------

This is solely a bugfix release.

* Fixed a bug where the dither did not have the correct width.
* Fixed a bug for cases with no dithering.
* Various minor improvements to the documentation