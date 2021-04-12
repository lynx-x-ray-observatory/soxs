.. _spatial:

Spatial Models in SOXS
======================

The ``SpatialModel`` class can be used to create RA and Dec positions of photons, 
which can be combined with the energies from a ``Spectrum`` object to create a 
source that can be written to a SIMPUT photon list. Several ``SpatialModel`` 
derivatives are available, which are documented below.

In general, each ``SpatialModel`` takes the following information:

1. A central RA and Dec for the source
2. Some prescription for how the photons should be distributed on the sky in
   terms of parameters and/or a model function

Each ``SpatialModel`` can be used to generate sky coordinates for events, as
described below in :ref:`generate-coords`.

``PointSourceModel``
--------------------

The :class:`~soxs.spatial.PointSourceModel` generates photon positions for a 
point source.

.. code-block:: python

    from soxs import PointSourceModel
    ra0 = 30.0 # source RA in degrees
    dec0 = 45.0 # source Dec in degrees
    pt_src = PointSourceModel(ra0, dec0)

Though this model is trivial, it is constructed in the same way as the other 
models below for consistency.

Radial Models
-------------

The following classes generate azimuthally symmetric models (though see 
:ref:`ellipticity`) from functions or lookup tables for a surface 
brightness profile as a function of radius.

``BetaModel``
+++++++++++++

The :class:`~soxs.spatial.BetaModel` generates photon positions for a 
:math:`\beta`-model profile, often used to model galaxy clusters. The 
functional form of the :math:`\beta`-model for a surface brightness 
profile is:

.. math::

    S(r) = S_0\left[1+\left(\frac{r}{r_c}\right)^2\right]^{(-3\beta+1/2)}

where :math:`S_0` is the central surface brightness, :math:`\beta` is 
the slope parameter, and :math:`r_c` is the core radius. To construct one:

.. code-block:: python

    from soxs import BetaModel
    ra0 = 30.0 # center RA in degrees
    dec0 = 45.0 # center Dec in degrees
    r_c = 20.0 # the core radius in arc seconds
    beta = 2./3. # the beta slope parameter
    beta_src = BetaModel(ra0, dec0, r_c, beta)

The normalization of the :class:`~soxs.spatial.BetaModel` will be determined
by the :class:`~soxs.spectra.Spectrum` object it is combined with, so the 
:math:`S_0` parameter is not specified. 

.. _double-beta-model:

``DoubleBetaModel``
+++++++++++++++++++

The :class:`~soxs.spatial.DoubleBetaModel` generates photon positions for a 
sum of two :math:`\beta`-model profiles, often used to model cool-core galaxy 
clusters. This sum is parameterized as:

.. math::

    S(r) = S_{0,1}\left\{\left[1+\left(\frac{r}{r_{c,1}}\right)^2\right]^{(-3\beta_1+1/2)} +
           \frac{S_{0,2}}{S_{0,1}}\left[1+\left(\frac{r}{r_{c,2}}\right)^2\right]^{(-3\beta_2+1/2)}\right\}

where :math:`S_{0,1}` and :math:`S_{0,2}` are the central surface brightness 
parameters of the two profiles, :math:`\beta_1` and :math:`\beta_2` are the 
slope parameters of the two profiles, and :math:`r_{c,1}` and :math:`r_{c,2}` are 
the core radius parameters. The ratio :math:`S_{0,2}/S_{0,1}` is parameterized by 
``sb_ratio`` in the example below. To construct a :class:`~soxs.spatial.DoubleBetaModel` 
object:

.. code-block:: python

    from soxs import DoubleBetaModel
    ra0 = 30.0 # center RA in degrees
    dec0 = 45.0 # center Dec in degrees
    r_c1 = 20.0 # the inner core radius in arc seconds
    beta1 = 2./3. # the inner beta slope parameter
    r_c2 = 100.0 # the outer core radius in arc seconds
    beta2 = 1. # the outer beta slope parameter
    sb_ratio = 0.5 # the ratio of the outer to the inner SB peak value
    beta_src = DoubleBetaModel(ra0, dec0, r_c1, beta1, r_c2, beta2,
                               sb_ratio)

``AnnulusModel``
++++++++++++++++

The :class:`~soxs.spatial.AnnulusModel` can be used to generate photon 
positions for a annulus or disk with uniform surface brightness:

.. code-block:: python

    from soxs import AnnulusModel
    ra0 = 30.0 # center RA in degrees
    dec0 = 45.0 # center Dec in degrees
    r_in = 0.0 # inner radius of shell in arcseconds
    r_out = 10.0 # outer radius of shell in arcseconds
    ann_src = AnnulusModel(ra0, dec0, r_in, r_out)


``RadialFunctionModel``
+++++++++++++++++++++++

:class:`~soxs.spatial.RadialFunctionModel` takes as input a central RA, 
Dec, and a Python function or callable object to generate an azimuthally 
symmetric distribution of photon positions:

.. code-block:: python

    from soxs import RadialFunctionModel
    # A simple inverse square-law surface brightness profile.
    # There is no need to normalize it properly, since that 
    # will be taken care of by the accompanying spectral 
    # model. r is in arcseconds.
    def S_r(r):
        return 1.0/(r*r)
    ra0 = 100.0 # center RA in degrees
    dec0 = -30.0 # center Dec in degrees
    my_src = RadialFunctionModel(ra0, dec0, S_r)

``RadialArrayModel``
++++++++++++++++++++

:class:`~soxs.spatial.RadialArrayModel` takes as input a central RA, 
Dec, and two NumPy arrays of radius and surface brightness to generate 
an azimuthally symmetric distribution of photon positions:

.. code-block:: python

    from soxs import RadialArrayModel
    ra0 = 100.0 # center RA in degrees
    dec0 = -30.0 # center Dec in degrees
    r = np.linspace(0.0, 100.0, 10000) # binned array of radii in arcseconds
    r_s = 100.0 # scale radius of arcseconds
    S_r = 1.0/((1.0+r/r_s)**2*(r/r_s)) # the surface brightness array
    my_src = RadialArrayModel(ra0, dec0, r, S_r)

``RadialFileModel``
+++++++++++++++++++

:class:`~soxs.spatial.RadialFileModel` takes as input a central RA, Dec, 
and an ASCII table of two columns, radius and surface brightness, to 
generate an azimuthally symmetric distribution of photon positions:

.. code-block:: python

    from soxs import RadialFileModel
    ra0 = 100.0 # center RA in degrees
    dec0 = -30.0 # center Dec in degrees
    my_src = RadialFileModel(ra0, dec0, "my_profile.dat")

.. _ellipticity:

Ellipticity of Radial Source Models
+++++++++++++++++++++++++++++++++++

Any of the radial source models listed above take two parameters, 
``ellipticity`` and ``theta``, which define the ellipticity of the 
model and the orientation of the ellipse, respectively. For example, 
to make an elliptical annulus source tilted 45 degrees from the horizontal:

.. code-block:: python

    from soxs import AnnulusModel
    ra0 = 30.0 # center RA in degrees
    dec0 = 45.0 # center Dec in degrees
    r_in = 10.0 # inner radius of shell in arcseconds
    r_out = 30.0 # outer radius of shell in arcseconds
    ellipticity = 0.5
    theta = 45.0
    ann_src = AnnulusModel(ra0, dec0, r_in, r_out, ellipticity=ellipticity)

where ``ellipticity`` will shrink the annulus (or other shape) in 
the y-direction if < 1 or will expand it in the y-direction if > 1. 

``RectangleModel``
------------------

The :class:`~soxs.spatial.RectangleModel` generates photon positions 
on the sky which fill a given rectangle shape, which can be optionally 
rotated through an angle:

.. code-block:: python

    from soxs import RectangleModel
    ra0 = 30.0 # center RA in degrees
    dec0 = 45.0 # center Dec in degrees
    width = 20.0 # width of the rectangle in arcseconds
    height = 10.0 # height of the rectangle in arcseconds
    theta = 20.0 # rotation angle of the rectangle in degrees
    fov_src = RectangleModel(ra0, dec0, fov, theta=theta)

Setting either the ``width`` or ``height`` parameter to 0.0 creates a line source.

"Field of View" Sources
-----------------------

The :class:`~soxs.spatial.FillFOVModel` generates photon positions on 
the sky which fill a given field of view:

.. code-block:: python

    from soxs import FillFOVModel
    ra0 = 30.0 # center RA in degrees
    dec0 = 45.0 # center Dec in degrees
    fov = 20.0 # width of the field of view in arcminutes
    fov_src = FillFOVModel(ra0, dec0, fov)

This may be useful for creating background-like sources.

.. _generate-coords:

Generating Event Coordinates from Spatial Models
------------------------------------------------

To generate coordinates from any :class:`~soxs.spatial.SpatialModel`, the method
:meth:`~soxs.spatial.SpatialModel.generate_coords` is provided. This method takes
the number of events you wish to generate as a required parameter, and a pseudo
random number generator as an optional parameter. It returns two unitful arrays of
RA and Dec coordinates in degrees:

.. code-block:: python

    from soxs import BetaModel
    ra0 = 30.0 # center RA in degrees
    dec0 = 45.0 # center Dec in degrees
    r_c = 20.0 # the core radius in arc seconds
    beta = 2./3. # the beta slope parameter
    beta_src = BetaModel(ra0, dec0, r_c, beta)
    
    # Generate coordinates
    prng = 24 # random seed
    num_events = 1000000 # number of events to generate
    ra, dec = beta_src.generate_coords(num_events, prng=prng)

Normally, :meth:`~soxs.spatial.SpatialModel.generate_coords` will not need to be 
called by the end-user but will be used "under the hood" in the generation of
a :class:`~soxs.simput.PhotonList` as part of a :class:`~soxs.simput.SimputCatalog`.
See :ref:`simput` for more information.
