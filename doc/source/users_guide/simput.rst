.. _simput:

Working with SIMPUT Catalogs
============================

The default storage format for creating 2D models of X-ray sources in SOXS
is `SIMPUT <https://www.sternwarte.uni-erlangen.de/research/sixte/simput.php>`_, 
which is fast becoming a standard for making mock X-ray observations. 

SOXS provides three classes to handle SIMPUT I/O: 
:class:`~soxs.simput.SimputSpectrum`, :class:`~soxs.simput.SimputPhotonList`, 
and :class:`~soxs.simput.SimputCatalog`. These utilize the spectral and spatial 
model classes described in the :ref:`spectra` and :ref:`spatial` sections of the
SOXS documentation. 

.. _simput-spectra:

SIMPUT Spectral Models
----------------------

The simplest SIMPUT source model in SOXS is that of a point-source located at a 
particular location on the sky, with a given energy spectrum. A 
:class:`~soxs.simput.SimputSpectrum` can be created using the 
:meth:`~soxs.simput.SimputSpectrum.from_spectrum` method, using a
:class:`~soxs.spectra.Spectrum` object:

.. code-block:: python

    import soxs
    
    # Create a thermal Spectrum
    agen = soxs.ApecGenerator(0.1, 10.0, 1000)
    spec1 = agen.get_spectrum(6.0, 0.3, 0.01, 1.0e-3)
    spec1.apply_foreground_absorption(0.02)

    # Pick some coordinates
    ra1 = 22.0 # degrees
    dec1 = -30.0 # degrees
    
    name1 = "ptsrc" # name of source

    # Create a SimputSpectrum
    src1 = soxs.SimputSpectrum.from_spectrum(name1, spec1, ra1, dec1)

An extended source can be created from a FITS image if you have one. Note that
the image must have a header with some coordinate information, where 
:math:`n \in {1,2}`:

* ``"CRPIXn"``: reference pixel x,y coordinates, usually the image
  center
* ``"CUNITn"``: both should contain ``"deg"``
* ``"CDELTn"``: width of each pixel in the x and y directions in 
  units of ``"CUNITn"``, ``"CDELT1"`` should be negative
* ``"CRVALn"``: reference celestial x,y coordinates--required, but 
  these will generally not be used. 
* ``"CTYPEn"``: must be a projection type, typically ``"RA---TAN"``
  and ``"DEC--TAN"`` are used. 

.. code-block:: python

    from astropy.io import fits
    
    # Create a thermal Spectrum
    spec2 = agen.get_spectrum(5.0, 1.0, 0.03, 1.0e-4)
    spec2.apply_foreground_absorption(0.02)
    
    # Pick some coordinates
    ra2 = 22.01 # degrees
    dec2 = -29.98 # degrees

    imhdu "cluster_image.fits[0]" # this specifies the name and extension of the 
                                  # image
    
    name2 = "cluster1" # name of source

    # Create a thermal SimputSpectrum
    src2 = soxs.SimputSpectrum.from_spectrum(name2, spec2, ra2, dec2, 
                                             imhdu=imhdu)
    
Note in this case that the entire extended source will have the same spectrum.

Alternatively, if you have a :class:`~soxs.spatial.SpatialModel`, you can use
the :meth:`~soxs.simput.SimputSpectrum.from_models` method to create an extended
source, where you use the :class:`~soxs.spatial.SpatialModel` to create an image
with a specific ``width`` and resolution ``nx``:

.. code-block:: python

    from astropy.io import fits
    
    # Create a Spectrum
    spec3 = agen.get_spectrum(2.2, 0.5, 0.05, 2.0e-2)
    spec3.apply_foreground_absorption(0.02)
        
    # Create an AnnulusModel
    ra3 = 22.03 # degrees
    dec3 = -30.03 # degrees
    r_in = 5.0 # arcseconds
    r_out = 20.0 # arcseconds
    ann = soxs.AnnulusModel(ra3, dec3, r_in, r_out)
    
    width = 20.0 # of the image, in arcminutes
    nx = 4000 # resolution of the image

    name3 = "cluster2" # name of source

    # Create a SimputSpectrum
    src3 = soxs.SimputSpectrum.from_models(name3, spec3, ann, width, nx)

In this case the whole extended source has the same spectrum as well. 

.. _photon-lists:

SIMPUT Photon List Models
-------------------------

Spectral and spatial models for X-ray sources can be combined to produce a list 
of photon coordinates and energies using the 
:class:`~soxs.simput.SimputPhotonList` class. Specifically, one can generate a 
:class:`~soxs.simput.SimputPhotonList` using the 
:meth:`~soxs.simput.SimputPhotonList.from_models` method. This requires a 
:class:`~soxs.spectra.Spectrum` for the spectral model, a 
:class:`~soxs.spatial.SpatialModel` for modeling the spatial extent of the 
source, an exposure time, and a flat effective area.

.. code-block:: python

    import soxs
     
    # Create the spectral model
    spec4 = soxs.Spectrum.from_powerlaw(1.0, 0.01, 1.0e-2, 0.1, 10.0, 100000)
    spec4.apply_foreground_absorption(0.04)
    
    # Create a RectangleModel
    ra4 = 21.97 # degrees
    dec4 = -30.0 # degrees
    width = 100.0 # in arcseconds
    height = 4.0 # in arcseconds
    theta = 30.0
    rect = soxs.RectangleSource(ra4, dec4, width, height, theta=theta)
    
    # Set the parameters
    exp_time = (500.0, "ks")
    area = (3.0, "m**2")

    name4 = "jet" # name of source

    # Create the photon list
    src4 = soxs.SimputPhotonList.from_models(name4, spec4, rect, exp_time, area)
                         
Plotting Photon Lists
+++++++++++++++++++++

The event positions from a :class:`~soxs.simput.SimputPhotonList` can be plotted
using the :meth:`~soxs.simput.SimputPhotonList.plot` method. This will make a 
scatter plot of the photon RA and Dec on the sky, optionally filtered within an 
energy band. For an example of how to use this method, see the 
:ref:`two-clusters` cookbook example.

.. _simput-catalogs:

SIMPUT Catalogs
---------------

A SIMPUT catalog can be worked with using the :class:`~soxs.simput.SimputCatalog`
class. A :class:`~soxs.simput.SimputCatalog` object associated with a single
:class:`~soxs.simput.SimputSource` can be created using the 
:meth:`~soxs.simput.SimputCatalog.from_source` method:

.. code-block:: python

    import soxs
        
    # Create the SIMPUT catalog
    sim_cat = SimputCatalog.from_source("my_sources.simput", src1, overwrite=True)

which writes both the catalog and the source to the same file 
``"my_sources.simput"``. If you want to write the 
:class:`~soxs.simput.SimputSource` to a separate file, use the ``src_filename``
keyword argument:

.. code-block:: python

    import soxs
        
    # Create the SIMPUT catalog
    sim_cat = SimputCatalog.from_source("my_sources.simput", src1, 
                                        src_filename="ptsrc.fits", 
                                        overwrite=True)

To add more sources to an existing :class:`~soxs.simput.SimputCatalog`, 
use the :meth:`~soxs.simput.SimputCatalog.append` method:

.. code-block:: python

    # This adds src2 to the catalog in "my_sources.simput" and the same file
    sim_cat.append(src2)
    
    # This adds src3 to the catalog in "my_sources.simput" and the file
    # cluster1.fits
    sim_cat.append(src3, src_filename="cluster1.fits", overwrite=True)

An existing SIMPUT catalog can be read in from disk using
:meth:`~soxs.simput.SimputCatalog.from_file`:

.. code-block:: python

    import soxs
    sim_cat = soxs.SimputCatalog.from_file("my_sources_simput.fits")
