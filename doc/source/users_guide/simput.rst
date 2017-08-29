.. _simput:

Working with SIMPUT Catalogs
============================

The default storage format for unconvolved events in SOXS is SIMPUT, which is 
fast becoming a standard for making mock X-ray observations. SOXS provides 
two classes to handle SIMPUT I/O, :class:`~soxs.simput.PhotonList` and 
:class:`~soxs.simput.SimputCatalog`.

.. _photon-lists:

Photon Lists
------------

Spectral and spatial models for X-ray sources can be combined to produce a list 
of photon coordinates and energies using the :class:`~soxs.simput.PhotonList` class. 
Specifically, one can generate a :class:`~soxs.simput.PhotonList` using the 
:meth:`~soxs.simput.PhotonList.from_models` method. This requires a 
:class:`~soxs.spectra.Spectrum` for the spectral model, a 
:class:`~soxs.spatial.SpatialModel` for modeling the spatial extent of the source,
an exposure time, and a flat effective area. The goal of producing a 
:class:`~soxs.simput.PhotonList` is to generate a large sample of candidate events
which can be used as a Monte-Carlo sample by either the
:func:`~soxs.instrument.instrument_simulator` or another tool such as MARX, SIMX, or
SIXTE to produce a mock X-ray observation. 

.. code-block:: python

    from soxs import PhotonList, PointSource, Spectrum
    
    # Create the spectral model
    spec = Spectrum.from_powerlaw(1.0, 0.01, 1.0e-2, 0.1, 10.0, 100000)
    spec.apply_foreground_absorption(0.04)
    
    # Create the spatial model
    pt_src = PointSource(30.0, 45.0)
    
    # Set the parameters
    exp_time = (500.0, "ks")
    area = (3.0, "m**2")

    # Create the photon list
    phlist = PhotonList.from_models('pt_src', spec, pt_src, exp_time, area)
                         
In this example, we've given the :class:`~soxs.simput.PhotonList` the name
`'pt_src'`, which will be used as the prefix for any file that is written
from this :class:`~soxs.simput.PhotonList`. To write the photon list to a
FITS file and the corresponding SIMPUT catalog file, use 
:meth:`~soxs.simput.PhotonList.write_photon_list`. In this case, no SIMPUT
catalog file yet exists, so a new one will be created:

.. code-block:: python

    simput_prefix = "my_sources"
    phlist.write_photon_list(simput_prefix, overwrite=True)

Alternatively, you may already have a SIMPUT file associated with a photon 
list file, but want to add another source to the same SIMPUT catalog. You can
accomplish this by making the same call to 
:meth:`~soxs.simput.PhotonList.write_photon_list`, but setting ``append=True``:

.. code-block:: python

    simput_prefix = "my_sources"
    phlist.write_photon_list(simput_prefix, append=True, overwrite=True)

The files written are ``"my_sources_simput.fits"`` and ``"pt_src_phlist.fits"``.

Plotting Photon Lists
+++++++++++++++++++++

The event positions from a :class:`~soxs.simput.PhotonList` can be plotted using
the :meth:`~soxs.simput.PhotonList.plot` method. This functionality requires the
`WCSAxes <http://wcsaxes.readthedocs.io/>`_ Python package to be installed. This
will make a scatter plot of the photon RA and Dec on the sky, optionally filtered
within an energy band. For an example of how to use this method, see the 
:ref:`two-clusters` cookbook example.

.. _simput-catalogs:

SIMPUT Catalogs
---------------

A SIMPUT catalog can be worked with directly using the :class:`~soxs.simput.SimputCatalog`
class. A :class:`~soxs.simput.SimputCatalog` object associated with a single
:class:`~soxs.simput.PhotonList` can be generated using the 
:meth:`~soxs.simput.SimputCatalog.from_models` method in the same way as the
:meth:`~soxs.simput.PhotonList.from_models` method:

.. code-block:: python

    from soxs import SimputCatalog, PointSource, Spectrum
    
    # Create the spectral model
    spec = Spectrum.from_powerlaw(1.0, 0.01, 1.0e-2, 0.1, 10.0, 100000)
    spec.apply_foreground_absorption(0.04)
    
    # Create the spatial model
    pt_src = PointSource(30.0, 45.0)
    
    # Set the parameters
    exp_time = (500.0, "ks")
    area = (3.0, "m**2")

    # Create the SIMPUT catalog
    sim_cat = SimputCatalog.from_models("my_sources", 'pt_src', spec, pt_src, 
                                        exp_time, area)

You can write this catalog and its photon list file to disk using 
:meth:`~soxs.simput.SimputCatalog.write_catalog`:

.. code-block:: python

    sim_cat.write_catalog(overwrite=True)

The files written are ``"my_sources_simput.fits"`` and ``"pt_src_phlist.fits"``.

An existing SIMPUT catalog can be read in from disk using
:meth:`~soxs.simput.SimputCatalog.from_file`:

.. code-block:: python

    import soxs
    sim_cat = soxs.SimputCatalog.from_file("my_sources_simput.fits")

If you have one or more :class:`~soxs.simput.PhotonList` objects, they
can be used to create a :class:`~soxs.simput.SimputCatalog` object on the
fly:

.. code-block:: python

    import soxs
    ...
    clusters = [cluster1, cluster2] # A list of PhotonList objects
    sim_cat = soxs.SimputCatalog("clusters", clusters)

Finally, an existing :class:`~soxs.simput.PhotonList` can be appended to an 
existing :class:`~soxs.simput.SimputCatalog` using 
:meth:`~soxs.simput.SimputCatalog.append`:

.. code-block:: python

    import soxs
    
    # Create the spectral model
    agen = soxs.ApecGenerator(0.05, 20.0, 200000)
    spec = agen.get_spectrum(4.0, 0.3, 0.05, 1.0e-3)
    spec.apply_foreground_absorption(0.04)
    
    # Create the spatial model
    beta_src = BetaModel(30.0, 45.0, 20.0, 1.666667)
        
    # Set the parameters
    exp_time = (500.0, "ks")
    area = (3.0, "m**2")

    # Create the photon list
    cluster = PhotonList.from_models('cluster', spec, beta_src, exp_time, area)

    sim_cat.append(cluster)