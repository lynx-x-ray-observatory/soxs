.. _changelog:

ChangeLog
=========

Version 4.8.5
-------------

This version of SOXS fixes two bugs.

* For certain instruments, in particular the *Lynx* LXM varieties, instrumental
  background simulations were yielding no events because of how the files were
  being read. This bug has now been fixed.
* Custom instrument files that were not part of the official SOXS instrument file
  registry, but placed in the ``soxs_data_dir`` location, were not being properly
  detected. This bug has been fixed. Thanks to `xshaokun <https://github.com/xshaokun>`_
  for pointing this out in
  `GitHub Issue #38 <https://github.com/lynx-x-ray-observatory/soxs/issues/38>`_.

Version 4.8.4
-------------

This version of SOXS fixes two bugs and one documentation error.

* Spatial region handling in the context of :func:`~soxs.events.filter_events` and
  :func:`~soxs.events.write_spectrum` has been fixed for composite regions which
  mixed regions which include and exclude data.
* When creating spectra using :func:`~soxs.instrument.simulate_spectrum` and
  including backgrounds, the square root of the ``bkgnd_area`` parameter was
  used to normalize the backgrounds instead of the ``bkgnd_area`` parameter itself.
  This has been fixed. Thanks to `liuguanfu1120 <https://github.com/liuguanfu1120>`_
  for pointing this out in
  `GitHub Issue #36 <https://github.com/lynx-x-ray-observatory/soxs/issues/36>`_.
* The :ref:`background` page incorrectly reported ``"wabs"`` as the default
  absorption model for the Milky Way foreground, instead of ``"tbabs"``. This has
  been corrected. Thanks to `liuguanfu1120 <https://github.com/liuguanfu1120>`_
  for pointing this out in
  `Github Issue #35 <https://github.com/lynx-x-ray-observatory/soxs/issues/35>`_.

Version 4.8.3
-------------

This version of SOXS fixes three LEM RMFs, corresponding to the 2.5 eV, 1.3 eV, and
1.2 eV resolution versions. The first two are used in the the ``lem_outer_array``
and ``lem_inner_array`` instrument specifications, respectively. The problem was that the
``EBOUNDS`` extension of these RMFs had the ``CHANNEL`` column 0-indexed, when it should
have been 1-indexed. This does not affect any spectral simulations, but it may affect
plotting routines in other packages that use spectra produced by SOXS. Thanks to Adam
Foster for raising this to our attention.

Version 4.8.2
-------------

This version of SOXS:

* Fixes a bug in which it was impossible to create an exposure map for instruments
  without dither.
* Drops support for Python 3.8.

Version 4.8.1
-------------

This version of SOXS contains a bugfix, where a spurious warning message about
using instrument files in the SOXS data directory is removed.

Version 4.8.0
-------------

This version of SOXS contains a bugfix, a small change in behavior, and new
response files for the XRISM instrument models.

* A bug that occurred when SOXS was not able to find the location of the
  ``"CHANTYPE"`` header keyword in RMFs when making backgrounds has now been
  fixed. Thanks to Charles Romero for pointing this out.
* When plotting an instrument spectrum using :func:`~soxs.events.plot_spectrum`,
  if ``plot_counts`` is set to True, the y-axis will now be in units of counts/bin
  instead of counts/keV regardless of whether ``plot_energy`` is ``True`` or
  ``False``.
* New response files for the XRISM instrument models have been added based on
  those provided for XRISM AO-1, which in turn has resulted in new XRISM
  instrument models. See :ref:`xrism` for details.

Version 4.7.1
-------------

This version of SOXS contains a critical bugfix. Occasionally, instrument files
such as the ARF, RMF, and instrumental background files would be updated to
improved versions with the same filename (this is particularly the case for the
ACIS-I/S Cycle 0 files). What this should mean is that if a user had previously
downloaded the files, they would be updated to the new version. However, the
code was not checking for the existence of new files with different checksums
properly. This has now been fixed. Users are encouraged to upgrade, and if they
have used the ACIS-I/S Cycle 0 files, to check that they have the latest versions.
The simplest way to do this is to do the following for each instrument of interest:

.. code-block:: python

    import soxs
    arf = soxs.AuxiliaryResponseFile.from_instrument("chandra_acisi_cy0")
    rmf = soxs.RedistributionMatrixFile.from_instrument("chandra_acisi_cy0")

Version 4.7.0
-------------

This version of SOXS contains new features and bugfixes.

* Python 3.12 is now supported.
* The foreground model normalization used in :func:`~soxs.simput.make_bkgnd_simput`
  was not being scaled appropriately by the field of view size. This has been
  fixed.
* New useful attributes for :class:`~soxs.spectra.Spectrum` objects have been
  added. See :ref:`spec-attribs` for more details.
* The default SPEX version for CIE spectra has been updated to 3.07.03.
* When loading an RMF, SOXS now checks the ``EBOUNDS`` header for the
  ``CHANTYPE`` keyword if it is not present in the ``MATRIX`` header.
  This fixes issues with the new LEM instruments released in version 4.6.0.
* It is now possible to vary the abundance of the hot halo components of
  the foreground model. See :ref:`foreground` and :ref:`config` for more
  details.
* It is now possible to append the SIMPUT sources produced by
  :func:`~soxs.simput.make_bkgnd_simput` to an existing SIMPUT catalog.

Version 4.6.0
-------------

This version of SOXS contains new features and bugfixes.

* It is now possible, in combination with pyXSIM version 4.3.0 or later, to
  use pyXSIM event lists written to HDF5 files as inputs to
  :func:`~soxs.instrument.instrument_simulator`. See :ref:`instrument` for
  details.
* A bug affecting the convolution of spectra with RMFs with multiple channel
  groups in the context of :func:`~soxs.instrument.simulate_spectrum` has been
  fixed.
* A new mode for :func:`~soxs.instrument.simulate_spectrum`, which allows for
  the instrument specification to be a 2 or 3-tuple specifying the ARF, RMF,
  and particle background, has been added. See :ref:`simulate-spectrum` for
  more details. This has not yet been enabled in the command-line interface
  version of ``simulate_spectrum``.
* It is now possible to include the simulation of the MW foreground and the
  CXB in :func:`~soxs.instrument.simulate_spectrum` for gratings instruments.
* It is now possible to use multiple, separated energy bands when extracting events for
  the construction of an image using :func:`~soxs.events.write_image`. See
  :ref:`write-image` for more details.
* A bug that prevented the making of exposure maps for *LEM* instruments has been fixed.
  Thanks to Arash Bodaghee for reporting this bug.
* *LEM* instrument configurations for 2.5 eV and 1.3 eV spectral resolution have been added,
  and should be considered the defaults, as these represent the baseline instrument. The
  old configurations are still present. See :ref:`lem` for more details.
* The response files used for the ``"chandra_acisi_cy0"`` and ``"chandra_aciss_cy0"``
  instruments suffered from the effects of poor calibration due to a high ACIS focal
  plane temperature. They have been replaced with response files that do not suffer
  from this issue.
* The particle background files for the *Chandra*/ACIS imaging instruments and the
  *Lynx* imaging instruments have been updated with minor changes.
* A new function to fill regions in an image where point sources have been removed,
  :func:`~soxs.events.fill_regions`, has been added. See :ref:`fill-regions` for
  more details.

Version 4.5.3
-------------

This version of SOXS contains two bugfixes related to generating mock observations
from SIMPUT catalogs which use FITS images for modeling photon positions:

* The input image did not have an extension name in the SIMPUT spectrum file, so
  SOXS was not able to read it. This has been fixed.
* The input image was not being rotated to the correct orientation. This has been
  fixed.

Version 4.5.2
-------------

This version of SOXS contains three bugfixes:

* For instruments with image-based PSFs, the PSF image was incorrectly transposed.
  Thankfully, this only affected the *XRISM*/Resolve instrument, since its image is
  rectangular and asymmetric. This has now been fixed.
* Default aimpoint coordinates corresponding to the detector center have now been
  added to simple, square-shaped instruments created with
  :meth:`~soxs.instrument_registry.make_simple_instrument`.
* SIMPUT filenames are now no longer limited to 80 characters inside the SIMPUT
  catalog, and better handling is provided for filenames with relative paths. Thanks
  to Chang-Goo Kim for submitting `PR 19 <https://github.com/lynx-x-ray-observatory/soxs/pull/19>`_ which fixes this.


Version 4.5.1
-------------

This version of SOXS contains two bugfixes:

* A critical downstream bug in pyXSIM has been fixed, where normalizations of X-ray
  fields, spectra, and mocks that used the IGM thermal spectrum model were overestimated.
  Users who need this functionality are also encouraged to upgrade to pyXSIM v4.2.0.
* Inputting a file name as the ``imhdu`` argument to
  :meth:`~soxs.simuput.SimputSource.from_spectrum` was not working, and has now been
  fixed.

Version 4.5.0
-------------

This version of SOXS contains a minor bugfix and a number of small new features.

* More corner cases of SIMPUT catalogs made using the SIMPUT library
  which caused errors in SOXS are now supported.
* It is now possible to supply a :class:`soxs.simput.SimputCatalog` instance
  as the ``input_events`` argument to :func:`~soxs.instrument.instrument_simulator`.
* It is now possible to specify values of the ``reblock`` parameter that are less
  than 1 to :func:`soxs.events.write_image`.
* It is now possible to filter events on time in :func:`soxs.events.filter_events`,
  :func:`soxs.events.write_image`, and :func:`soxs.events.write_spectrum`.
* It is now possible to exclude events with region filters in :func:`soxs.events.filter_events`
  and :func:`soxs.events.write_spectrum`.
* A new function to merge source and background event files,
  :func:`soxs.events.merge_event_files`, has been added.

Version 4.4.0
-------------

This version of SOXS contains critical bugfixes and one new feature.

* There was an `off-by-one` indexing error in the production of energies for diffuse
  background spectra, as well as any spectra produced with
  :func:`~soxs.instrument.simulate_spectrum`, which results in a small energy shift
  (almost always below the energy resolution). This bug has been fixed.
* The ``"ENERGY"`` column in event files produced by SOXS now represent the energies that
  are approximated by the instrument response based on their channel. Effectively, this
  now means that these energies are at the instrument resolution. This is in line with
  what is present in real data. A new column in the event files, ``"SOXS_ENERGY"``, contains
  the energies incident on the detector derived from the source, which previously were
  in the ``"ENERGY"`` column.
* Region files or expressions with multiple regions inside them are now correctly
  parsed when using :func:`~soxs.events.filter_events` or :func:`~soxs.events.write_spectrum`.
* It is now possible to create a spectrum without Poisson noise using
  :func:`~soxs.instrument.simulate_spectrum` or the ``simulate_spectrum`` command-line
  script. See :ref:`simulate-spectrum` or :ref:`cmd-simulate-spectrum` for more details.
* The ``"CHANNEL"`` field in the ``"EBOUNDS"`` data in the LEM RMFs was 0-indexed when it
  should have been 1-indexed. This has been fixed.

Version 4.3.0
-------------

This version of SOXS contains new features.

* A new version of the spectral model used in the
  :class:`~soxs.thermal_spectra.CloudyCIEGenerator` class has been provided, with
  improved energy resolution. See :ref:`cloudy-spectra` for more details.
* A new version of the spectral model used in the
  :class:`~soxs.thermal_spectra.IGMGenerator` class has been provided, with
  improved energy resolution. See :ref:`igm-spectra` for more details.
* A new function to download table files for the thermal spectra models has been
  provided. See :ref:`downloading-thermal-tables` for more details.

Version 4.2.1
-------------

This update to SOXS contains bugfixes.

* The *AXIS* instrument specification was not working properly due to an issue
  with the implementation of the PSF file. This has now been fixed.
* In several places, data from FITS files is now converted to the native byteorder
  of the system upon reading.
* The minimum AstroPy version is now 4.0 and the minimum h5py version is now 3.0.

Version 4.2.0
-------------

This update to SOXS contains new features and a bugfix.

* Installation and use on Windows 64-bit platforms is now supported.
* New PSF models using encircled energy fraction (EEF) files are now supported.
  See :ref:`psf-models` for more details.
* The *XRISM* *Resolve* instrument specification has been updated, and a new
  instrument specification for *Xtend* has been added. See :ref:`xrism` for
  more details.
* The *AXIS* instrument specification has been updated. See :ref:`axis-probe` for
  more details.
* If one had not binned a :class:`~soxs.spectra.Spectrum` object more finely
  than the instrument's ARF/RMF when using :func:`~soxs.instrument.simulate_spectrum`,
  then gaps would appear in the resulting convolved spectrum. This is now
  handled by linearly interpolating the spectral model into the ARF energy
  bins.
* The *LEM* instrumental background has been boosted to 1 counts/s/keV/(30'x30')
  from the previous value of 0.07 counts/s/keV/(30'x30').
* A new function for creating a SIMPUT catalog including models for the Galactic
  foreground and the CXB point sources, :func:`~soxs.simput.make_bkgnd_simput`,
  has been included. See :ref:`bkgnd-simput` for more details.

Version 4.1.0
-------------

This update to SOXS contains bug fixes and two new features.

* A bug that scaled the flux of :class:`~soxs.simput.SimputSpectrum` sources
  incorrectly has been corrected.
* Bugs that prevented :class:`~soxs.simput.SimputSpectrum` sources from being
  used in SIXTE, SIMX, and MARX have been fixed.
* It is now possible to specify a region file with creating a spectrum with
  :func:`~soxs.events.write_spectrum`, to select a subset of events based on
  spatial region. See :ref:`write-spectrum` for more details.
* The method :meth:`~soxs.spectrum.Spectrum.get_lum_in_band` to compute the
  rest-frame luminosity of a :class:`~soxs.spectrum.Spectrum` within an energy
  band has been added.

Version 4.0.0
-------------

This update to SOXS contains a large number of new features, mostly related to
the generation of spectra.

* New options have been added for the simulation of thermal spectra, including
  from `SPEX <https://www.sron.nl/astrophysics-spex>`_, MeKaL, a CIE model based
  on `Cloudy <https://gitlab.nublado.org/cloudy/cloudy/-/wikis/home>`_,
  and a model for emission from the IGM including photoionization and resonant
  scattering off of the CXB based on Cloudy and provided by Ildar Khabibullin.
  See :ref:`thermal-spectra` for details.
* The option to create :class:`~soxs.spectra.Spectrum` objects with log-spaced
  energy binning has been added. See :ref:`spectrum-binning` for details.
* The option to create a new spectrum from an old one by rebinning has been added
  to the :class:`~soxs.spectra.Spectrum` class. See :ref:`spectrum-binning` for details.
* It is no longer necessary to source the HEADAS environment before creating a
  :class:`~soxs.spectra.Spectrum` object using either the
  :meth:`~soxs.spectra.Spectrum.from_xspec_script` or
  :meth:`~soxs.spectra.Spectrum.from_xspec_model`. See :ref:`xspec` for more details.
* Reading and writing of :class:`~soxs.spectra.Spectrum` objects has been refactored,
  so that the tables use the min and max of each energy bin instead of the middle
  energy of the bin. This allows for log-spaced energy binning (mentioned above) to
  be supported. Also, :class:`~soxs.spectra.Spectrum` objects can now be written to
  FITS table files as well as ASCII and HDF5. See :ref:`read-spectra` and
  :ref:`write-spectra` for details.
* An option to create a mosaicked event file in addition to an image file has been
  added to the :func:`~soxs.mosaic.make_mosaic_image` function. See :ref:`mosaic`
  for more details.
* The default absorption model for the galactic foreground has been changed to TBabs.
* The accuracy of the TBabs absorption model has been improved.
* It is now possible to specify different abundance tables in the construction of the
  TBabs absorption model.
* The galactic foregroud model now includes thermal broadening of emission lines, and
  it is also now possible to optionally add velocity broadening. See :ref:`foreground`
  for more details.
* The LEM ARF has been updated.
* Instrumental background models have been added to the LEM instrument models.
* The abundance table from `Feldman (1992) <https://ui.adsabs.harvard.edu/abs/1992PhyS...46..202F>`_
  has been added to the options for abundance tables for the
  :class:`~soxs.thermal_spectra.ApecGenerator` and :class:`~soxs.thermal_spectra.SpexGenerator`.
* The default abundance table from Cloudy v17.03 has been added to the options for abundance
  tables for the :class:`~soxs.thermal_spectra.ApecGenerator` and
  :class:`~soxs.thermal_spectra.SpexGenerator`.
* The command-line script ``make_thermal_spectrum`` has been changed to ``make_cie_spectrum`` and
  has many more options for computing CIE spectra. See :ref:`cmd-make-cie-spectrum` for details.
* The command-line script ``make_igm_spectrum`` has been added for making thermal spectra with
  photoionization and resonant scattering. See :ref:`cmd-make-igm-spectrum` for details.
* In the command-line scripts ``make_cie_spectrum``, ``make_igm_spectrum``, and
  ``make_powerlaw_spectrum``, the parameter for foreground Galactic absorption ``nh`` has been
  renamed to ``nH_abs``.

Version 3.4.0
-------------

* The LEM response files have been updated.
* The XRISM response files have been updated.
* A bug that prevented multi-image PSF types to be used has been fixed.
* The astrophysical foreground calculation method has been updated so that
  its spectral bins always match the binning of the RMF for the instrument
  being simulated, which results in more accurate spectral fits for this
  component.
* The point-source background :math:`\log N-\log S` distribution has been
  extended to fluxes up to :math:`S \sim 10^{-12}~\rm{erg}~\rm{s}~\rm{cm}^{-2}`.
  See :ref:`ptsrc-bkgnd` for more details.
* A diffuse component to the point-source background has been added, to model
  completely unresolved sources at low flux. See :ref:`ptsrc-bkgnd` for more
  details.

Version 3.3.0
-------------

* New instrument specifications for the
  `LEM probe concept <https://lem.physics.wisc.edu>`_ have been added, for
  spectral resolutions of 0.9 eV and 2 eV.
* A new function for filtering event files, :func:`~soxs.events.filter_events`,
  has been added. See :ref:`filtering-events` for more details.
* A number of small bugs have been fixed.


Version 3.2.0
-------------

* More customizations to the astrophysical X-ray foreground are now available.
  See :ref:`foreground` and :ref:`config` for details.
* A new function, :func:`~soxs.utils.set_soxs_config`, for setting configuration
  values, has been added. See :ref:`config` for details.
* A new function, :func:`~soxs.utils.set_mission_config`, for setting
  mission-specific configuration values, has been added. See :ref:`mission-config`
  for details.


Version 3.1.0
-------------

* Python 3.10 is now officially supported. The minimum supported Python version
  is now 3.8.
* For :func:`~soxs.instrument.simulate_spectrum`, the power-law index for the
  unresolved point-source component of the astrophysical background is now
  :math:`\alpha = 1.52`.
* An instrument specification for the
  `LEM probe concept <https://lem.physics.wisc.edu>`_ has been added.
* A bug which prevented the use of the ``xrism_resolve`` instrument has been
  fixed.
* The default neutral hydrogen column for the astrophysical background components
  is now :math:`n_H = 0.018 \times 10^{22}~\rm{atoms}~\rm{cm}^{-2}`
* The default value of the neutral hydrogen column and the absorption model for
  astrophysical backgrounds can now be set in :ref:`config`. These can no longer
  be fine tuned in :func:`~soxs.instrument.make_background_file` or
  :func:`~soxs.instrument.simulate_spectrum`.
* The default APEC version can now be set in the :ref:`config`.
* The keyword argument ``input_pt_sources`` has been added to the
  :func:`~soxs.instrument.instrument_simulator`, to allow a consistent set of
  point sources to be simulated. See :ref:`point-source-list` for information
  on how to create this file. The keyword argument to do the same for
  :func:`~soxs.instrument.make_background_file` is now renamed to
  ``input_pt_sources`` from ``input_sources`` for consistency.

Version 3.0.2
-------------

This version of SOXS contains bug fixes and a minor new feature.

* Relative paths are now handled correctly in SIMPUT catalogs.
* A number of problems in parsing instrument specifications have been fixed.
* A bug which caused a crash when an RMF with ``N_CHAN`` = 0 in columns has
  been fixed.
* :class:`~soxs.spectra.ConvolvedSpectrum` objects can now be added and
  subtracted.
* Doc examples which use pyXSIM now use pyXSIM 3.0.0.

Version 3.0.1
-------------

This bugfix update to SOXS contains bug fixes and a minor new feature.

* A bug which prevented SIMPUT photon lists written by SOXS to be read in by
  SIXTE has been fixed.
* A bug which prevented the use of instrumental background files which do not
  contain the ``"EXPOSURE"`` keyword in the header has been fixed.
* :func:`~soxs.instrument_registry.add_instrument_to_registry` now catches
  more errors in the setup of custom instruments and flags them informatively.
* Subtraction of two :class:`~soxs.spectra.Spectrum` objects is now possible.

Version 3.0.0
-------------

This major version update of SOXS contains new features and optimizations.
**NOTE: there are some backwards-incompatible changes in this release.**

* SOXS now supports two new PSF model types, ``"image"``, which uses a single
  FITS image for the PSF model, and ``"multi_image"``, which can use a number
  of FITS images corresponding to different incident photon energies and
  different off-axis angles. See :ref:`psf-models` for details.
* SOXS now uses standard PHA files with FITS tables of channel and count rate
  to create instrumental/particle background. See :ref:`instr-bkgnd` for more
  details.
* SOXS now supports "spectrum" SIMPUT sources, with and without images, for
  generating mock observations. See :ref:`simput` for details.
* SOXS now uses the "spectrum" SIMPUT sources in all of the command line scripts
  which create spatial models, so the signatures of those scripts have changed.
  See :ref:`cmd-spatial` for details.
* The Python function :meth:`~soxs.background.point_sources.make_point_source_list`
  and the command line script :ref:`cmd-make-point-source-list` no longer require
  the ``exp_time`` and ``area`` arguments.
* The capability to create mosaics of multiple SOXS event files into a single
  FITS image, with options for exposure correction, has been added. See
  :ref:`mosaic` for details.
* It is no longer necessary to download response files manually, as response
  files will now be downloaded automatically if they are needed to simulate
  an observation or create a spectrum. See :ref:`response_files` for details.
* Similarly, the latest version of the APEC tables is no longer packaged with
  SOXS, both the CIE and NEI versions of the APEC tables will be downloaded
  automatically if they are needed. See :ref:`thermal-spectra` for details.
* SOXS now uses the `new 201-bin temperature files from AtomDB
  <http://www.atomdb.org/download.php>`_ for thermal spectrum models.
* The SOXS configuration option ``"response_path"`` has been changed to
  ``"soxs_data_dir"`` and the former is deprecated. See :ref:`config`
  for details.
* All instrument specifications must now specifically implement at least one
  chip explicitly, so it is no longer permissible to specify the ``"chips"``
  argument to be ``None``. **This is a backwards-incompatible change.**
* Generating the galactic foreground and the instrumental background is now
  faster and uses less memory.
* Exposure map calculation now uses far less memory and is slightly faster.
* New options have been added to the :func:`~soxs.events.plot_spectrum` function.
  See :ref:`plot-spectrum` for details.
* *Chandra* grating responses for ACIS-S have been updated to Cycle 22.
* SOXS now uses the
  `AstroPy Regions package <https://astropy-regions.readthedocs.io/en/latest/>`_
  for region handling.
* An option for writing ds9 regions corresponding to the sky positions and sizes
  of the halos from events created from the cosmology source catalog has been
  added. See :ref:`cosmo-source-catalog` for more details.
* The :class:`~soxs.spatial.DoubleBetaModel` spatial source model has been
  added (see :ref:`double-beta-model`), along with the command-line script
  :ref:`cmd-make-double-beta-model-source`.
* The ``make_beta_model`` command line script has been renamed to
  :ref:`cmd-make-beta-model-source`.
* An instrument specification for the
  `*STAR-X* mission concept <https://ui.adsabs.harvard.edu/abs/2017SPIE10399E..08M/abstract>`_
  has been added.

Version 2.3.0
-------------

This version of SOXS fixes a few bugs, updates instrument specifications, and
adds a speedup for certain types of RMF convolutions. **Support for Python 2.7
has been dropped in this release.**

* Fixed an issue in cosmological sources generation where comoving units were
  assumed when they should have been proper.
* Fixed an issue where the ``make_thermal_spectrum`` command-line script had the
  wrong default version of APEC specified.
* *Chandra* imaging responses for ACIS-I and ACIS-S have been updated to Cycle 22.
* The *XRISM*/Resolve ARF has been updated to use a version with higher effective
  area. See the :ref:`xrism` section of :ref:`instrument` for more information.
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
