.. _source-catalogs:

Simulating Source Catalogs
==========================

SOXS has the capability of simulating two different types of source catalogs:
X-rays from a population of halos (galaxies, galaxy groups, and galaxy clusters),
and point sources. 

.. _cosmo-source-catalog:

Cosmological Source Catalog
---------------------------

SOXS provides a routine to generate X-ray photons from cosmologically distant 
sources. This is made possible using a halo catalog from 

The halo catalog is extracted from a light cone simulation of 100 square degrees
and a maximum redshift of z = 3. A low-mass cut has been made at 
:math:`M_{500c} = 3 \times 10^{12}~M_\odot`.

Using scaling relations from , the halo temperature and flux are derived from the
halo mass. The halos are then given random ellipticities and orientations. 

The cosmological parameters for this halo catalog are:

* :math:`H_0 = 0.7`
* :math:`\Omega_m = 0.279`
* :math:`\Omega_\Lambda = 0.721`
* :math:`\Omega_b = 0.0463`
* :math:`w_{\rm DE} = -1`
* :math:`\sigma_8 = 0.823`

:func:`~soxs.cosmology.make_cosmological_sources_file` generates a photon list
file for a SIMPUT catalog using the cosmological sources model:

.. code-block:: python

    simput_prefix = "my_bkgnd"
    phlist_prefix = "cosmo"
    exp_time = 500000.0 # seconds
    fov = 20.0 # arcmin
    sky_center = [30.0, 45.0] # RA, Dec in degrees
    nH = 0.02 # Foreground galactic absorption, optional
    area = 40000.0 # Flat collecting area to generate photon sample
    make_cosmological_sources_file(simput_prefix, phlist_prefix, exp_time, 
                                   fov, sky_center, nH=nH, area=area)

By default, a random position will be chosen within the halo catalog. If you 
would prefer to simulate a specific region within the catalog, set the keyword
argument ``cat_center`` to a particular coordinate between [-5, 5] degrees in 
either direction:

.. code-block:: python

    cat_center = [-0.2, 3.0]
    make_cosmological_sources_file(simput_prefix, phlist_prefix, exp_time, 
                                   fov, sky_center, nH=nH, area=area, 
                                   cat_center=cat_center, append=True)

.. _point-source-catalog:

Point Source Catalog
--------------------

SOXS also provides a function to create a SIMPUT catalog of point-sources. 
It is not necessary to do this for including point sources as a background
component in SOXS, as this will be done automatically, but it may be useful 
if you would like to tweak parameters of the sources, store the positions and
fluxes of the sources generated, or use the SIMPUT catalog in another simulation
program such as MARX or SIMX. 

:func:`~soxs.background.point_sources.make_point_sources_file` generates a
photon list file for a SIMPUT catalog using the point-source background model
described in :ref:`ptsrc-bkgnd`:

.. code-block:: python

    simput_prefix = "my_bkgnd"
    phlist_prefix = "pt_src"
    exp_time = 500000.0 # seconds
    fov = 20.0 # arcmin
    sky_center = [30.0, 45.0] # RA, Dec in degrees
    nH = 0.02 # Foreground galactic absorption, optional
    area = 40000.0 # Flat collecting area to generate photon sample
    make_point_sources_file(simput_prefix, phlist_prefix, exp_time, fov, 
                            sky_center, nH=nH, area=area)


.. note::

    For both cosmological and point sources, As with other SIMPUT catalogs, if you
    supply a value for ``simput_prefix`` to this function that refers to an existing
    catalog and set ``append=True``, the photon list file will be appended to an 
    existing SIMPUT catalog.

