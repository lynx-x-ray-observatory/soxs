.. _cmd-spatial:

Command Line Scripts for Spatial Models
=======================================

These are scripts that create photon lists for SIMPUT catalogs that can serve
as inputs to the instrument simulator. For details on what's going on under the 
hood, see :ref:`spatial` and :ref:`spectra`.

Each of these scripts accepts a ``specfile`` argument, which is a two-column ASCII
table of energy (keV) and specific flux (:math:`\rm{photons~s^{-1}~cm^{-2}~keV^{-1}}`) 
which can be written using the spectra command-line scripts (see :ref:`cmd-spectra`), 
the Python method :meth:`~soxs.spectrum.Spectrum.write_file` (see :ref:`spectra`), 
or can be created by hand.

``make_point_source``
---------------------

This script creates a SIMPUT photon list of a point source from a spectrum supplied in a
file.

.. code-block:: text

    usage: make_point_source [-h] [--area AREA] [--append] [--clobber]
                             simput_prefix phlist_prefix ra0 dec0 specfile
                             exp_time
    
    Create a SIMPUT photon list of a point source from a spectrum supplied in a
    file.
    
    positional arguments:
      simput_prefix  The prefix of the SIMPUT file to be used as the root of the
                     catalog. If it does not exist, it will be created.
      phlist_prefix  The prefix of the photon list file to be written.
      ra0            The right ascension of the source in degrees.
      dec0           The declination of the source in degrees.
      specfile       The file containing the spectrum to be used.
      exp_time       The exposure time to use, in seconds.
    
    optional arguments:
      -h, --help     show this help message and exit
      --area AREA    The collecting area to use, in cm^2. Default: 30000.0
      --append       If set, append a new source an existing SIMPUT catalog.
      --clobber      Whether or not to clobber an existing file with the same
                     name.
                     
Examples
++++++++

Make a brand-new SIMPUT catalog for a point source photon list, assuming 100 ks exposure. 

.. code-block:: bash

    [~]$ make_point_source pt_src src1 20.0 -32.0 pt_src_spectrum.dat 100000. --clobber

Add a new point source to an existing SIMPUT catalog, assuming 50 ks exposure. 

.. code-block:: bash

    [~]$ make_point_source pt_src src2 19.0 -31.0 pt_src_spectrum.dat 50000. --append --clobber

Specify a different collecting area for the photons. 

.. code-block:: bash

    [~]$ make_point_source pt_src src1 20.0 -32.0 pt_src_spectrum.dat 100000. --area=40000. --clobber

``make_beta_model``
-------------------

This script creates a SIMPUT photon list of a :math:`\beta`-model from a spectrum supplied in a
file. The functional form of the :math:`\beta`-model for a surface brightness profile is:

.. math::

    S(r) = S_0\left[1+\left(\frac{r}{r_c}\right)^2\right]^{(-3\beta+1/2)}

.. code-block:: text

    usage: make_beta_model [-h] [--area AREA] [--theta THETA]
                           [--ellipticity ELLIPTICITY] [--append] [--clobber]
                           simput_prefix phlist_prefix ra0 dec0 r_c beta specfile
                           exp_time
    
    Create a SIMPUT photon list of a beta-model source.
    
    positional arguments:
      simput_prefix         The prefix of the SIMPUT file to be used as the root
                            of the catalog. If it does not exist, it will be
                            created.
      phlist_prefix         The prefix of the photon list file to be written.
      ra0                   The right ascension of the source center in degrees.
      dec0                  The declination of the source center in degrees.
      r_c                   The core radius in arcseconds.
      beta                  The beta parameter.
      specfile              The file containing the spectrum to be used.
      exp_time              The exposure time to use, in seconds.
    
    optional arguments:
      -h, --help            show this help message and exit
      --area AREA           The collecting area to use, in cm^2. Default: 30000.0
      --theta THETA         The angle through which to rotate the beta model in
                            degrees. Only makes sense if ellipticity is added.
                            Default: 0.0
      --ellipticity ELLIPTICITY
                            The ellipticity of the radial profile, expressed as
                            the ratio between the length scales of the x and y
                            coordinates. The value of this parameter will shrink
                            or expand the profile in the direction of the "y"
                            coordinate, so you may need to rotate to get the shape
                            you want. Default: 1.0
      --append              If set, append a new source an existing SIMPUT
                            catalog.
      --clobber             Whether or not to clobber an existing file with the
                            same name.
                     
Examples
++++++++

Make a brand-new SIMPUT catalog for a :math:`\beta`-model photon list, assuming 100 ks exposure. 

.. code-block:: bash

    [~]$ make_beta_model my_srcs beta_src1 20.0 -32.0 10.0 1.0 my_spectrum.dat 100000. --clobber

Add a new :math:`\beta`-model to an existing SIMPUT catalog, assuming 50 ks exposure. 

.. code-block:: bash

    [~]$ make_beta_model my_srcs beta_src2 19.0 -31.0 10.0 1.0 my_spectrum.dat 50000. --append --clobber

Specify a different collecting area for the photons. 

.. code-block:: bash

    [~]$ make_beta_model my_srcs beta_src1 20.0 -32.0 10.0 1.0 my_spectrum.dat 100000. --area=50000. --clobber

Add ellipticity and tilt the model:

.. code-block:: bash

    [~]$ make_beta_model my_srcs beta_src1 20.0 -32.0 10.0 1.0 my_spectrum.dat 100000. --ellipticity=0.5 --theta=45.0 --clobber

``make_annulus_source``
-----------------------

This script creates a SIMPUT photon list of an annulus or disk with constant surface brightness
from a spectrum supplied in a file.

.. code-block:: text

    usage: make_annulus_source [-h] [--theta THETA] [--ellipticity ELLIPTICITY]
                               [--area AREA] [--append] [--clobber]
                               simput_prefix phlist_prefix ra0 dec0 r_in r_out
                               specfile exp_time
    
    Create a SIMPUT photon list of an annulus source with uniform surface
    brightness from a spectrum supplied in a file.
    
    positional arguments:
      simput_prefix         The prefix of the SIMPUT file to be used as the root
                            of the catalog. If it does not exist, it will be
                            created.
      phlist_prefix         The prefix of the photon list file to be written.
      ra0                   The right ascension of the source center in degrees.
      dec0                  The declination of the source center in degrees.
      r_in                  The inner annulus of the source center in arcseconds.
      r_out                 The outer annulus of the source center in arcseconds.
      specfile              The file containing the spectrum to be used.
      exp_time              The exposure time to use, in seconds.
    
    optional arguments:
      -h, --help            show this help message and exit
      --theta THETA         The angle through which to rotate the beta model in
                            degrees. Only makes sense if ellipticity is added.
                            Default: 0.0
      --ellipticity ELLIPTICITY
                            The ellipticity of the radial profile, expressed as
                            the ratio between the length scales of the x and y
                            coordinates. The value of this parameter will shrink
                            or expand the profile in the direction of the "y"
                            coordinate, so you may need to rotate to get the shape
                            you want. Default: 1.0
      --area AREA           The collecting area to use, in cm^2. Default: 30000.0
      --append              If set, append a new source an existing SIMPUT
                            catalog.
      --clobber             Whether or not to clobber an existing file with the
                            same name.

Examples
++++++++

Make a brand-new SIMPUT catalog for an annulus photon list, assuming 100 ks exposure. 

.. code-block:: bash

    [~]$ make_annulus_source my_srcs ann_src1 20.0 -32.0 0.0 30.0 my_spectrum.dat 100000. --clobber

Add a new annulus model to an existing SIMPUT catalog, assuming 50 ks exposure. 

.. code-block:: bash

    [~]$ make_annulus_source my_srcs ann_src2 19.0 -31.0 0.0 30.0 my_spectrum.dat 50000. --append --clobber

Specify a different collecting area for the photons. 

.. code-block:: bash

    [~]$ make_annulus_source my_srcs ann_src1 20.0 -32.0 0.0 30.0 my_spectrum.dat 100000. --area=50000. --clobber

Add ellipticity and tilt the model:

.. code-block:: bash

    [~]$ make_annulus_source my_srcs ann_src1 20.0 -32.0 0.0 30.0 my_spectrum.dat 100000. --ellipticity=2.0 --theta=30.0 --clobber

``make_rectangle_source``
-------------------------

This script creates a SIMPUT photon list of a rectangle shape with constant surface brightness
from a spectrum supplied in a file.

.. code-block:: text

    usage: make_rectangle_source [-h] [--theta THETA] [--area AREA] [--append]
                                 [--clobber]
                                 simput_prefix phlist_prefix ra0 dec0 width height
                                 specfile exp_time

    Create a SIMPUT photon list of a uniformly filled rectangle source from a
    spectrum supplied in a file.

    positional arguments:
      simput_prefix  The prefix of the SIMPUT file to be used as the root of the
                     catalog. If it does not exist, it will be created.
      phlist_prefix  The prefix of the photon list file to be written.
      ra0            The right ascension of the source center in degrees.
      dec0           The declination of the source center in degrees.
      width          The width of the rectangle in arcseconds.
      height         The width of the rectangle in arcseconds.
      specfile       The file containing the spectrum to be used.
      exp_time       The exposure time to use, in seconds.

    optional arguments:
      -h, --help     show this help message and exit
      --theta THETA  The angle through which to rotate the rectangle in degrees.
                     Default: 0.0
      --area AREA    The collecting area to use, in cm^2. Default: 30000.0
      --append       If set, append a new source an existing SIMPUT catalog.
      --clobber      Whether or not to clobber an existing file with the same
                     name.

Examples
++++++++

Make a brand-new SIMPUT catalog for a rectangle photon list, assuming 100 ks exposure.

.. code-block:: bash

    [~]$ make_rectangle_source my_srcs rect_src1 20.0 -32.0 20.0 10.0 my_spectrum.dat 100000. --clobber

Make the same rectangle, but rotate it by 30.0 degrees.

.. code-block:: bash

    [~]$ make_rectangle_source my_srcs rect_src1 20.0 -32.0 20.0 10.0 my_spectrum.dat 100000. --theta=30.0 --clobber

Create a line source with the same width and rotation angle.

.. code-block:: bash

    [~]$ make_rectangle_source my_srcs rect_src1 20.0 -32.0 20.0 0.0 my_spectrum.dat 100000. --theta=30.0 --clobber

Add a new rectangle model to an existing SIMPUT catalog, assuming 50 ks exposure.

.. code-block:: bash

    [~]$ make_rectangle_source my_srcs rect_src2 19.0 -31.0 20.0 10.0 my_spectrum.dat 50000. --append --clobber

Specify a different collecting area for the photons.

.. code-block:: bash

    [~]$ make_rectangle_source my_srcs rect_src1 20.0 -32.0 20.0 10.0 my_spectrum.dat 100000. --area=50000. --clobber

``make_fov_source``
-------------------

This script creates a SIMPUT photon list of a field of view with constant surface brightness
from a spectrum supplied in a file.

.. code-block:: text

    usage: make_fov_source [-h] [--area AREA] [--append] [--clobber]
                           simput_prefix phlist_prefix ra0 dec0 fov specfile
                           exp_time
    
    Create a SIMPUT photon list of a uniformly filled field of view source from a
    spectrum supplied in a file.
    
    positional arguments:
      simput_prefix  The prefix of the SIMPUT file to be used as the root of the
                     catalog. If it does not exist, it will be created.
      phlist_prefix  The prefix of the photon list file to be written.
      ra0            The right ascension of the source center in degrees.
      dec0           The declination of the source center in degrees.
      fov            The field of view on a side in arcminutes.
      specfile       The file containing the spectrum to be used.
      exp_time       The exposure time to use, in seconds.
    
    optional arguments:
      -h, --help     show this help message and exit
      --area AREA    The collecting area to use, in cm^2. Default: 30000.0
      --append       If set, append a new source an existing SIMPUT catalog.
      --clobber      Whether or not to clobber an existing file with the same
                     name.

Examples
++++++++

Make a brand-new SIMPUT catalog for a field-of-view photon list, assuming 100 ks exposure. 

.. code-block:: bash

    [~]$ make_fov_source my_srcs fov_src1 20.0 -32.0 20.0 my_spectrum.dat 100000. --clobber

Add a new field-of-view model to an existing SIMPUT catalog, assuming 50 ks exposure. 

.. code-block:: bash

    [~]$ make_fov_source my_srcs fov_src2 19.0 -31.0 20.0 my_spectrum.dat 50000. --append --clobber

Specify a different collecting area for the photons. 

.. code-block:: bash

    [~]$ make_fov_source my_srcs fov_src1 20.0 -32.0 20.0 my_spectrum.dat 100000. --area=50000. --clobber

``make_phlist_from_ascii``
--------------------------

This script takes a table of photon RA, Dec, and energies from an ASCII-formatted table and writes them
to a new SIMPUT photon list file.

.. code-block:: text

    usage: make_phlist_from_ascii [-h] [--append] [--clobber]
                                  simput_prefix phlist_prefix infile

    Create a SIMPUT photon list from an ASCII table of positions and energies. The file must contain
    the total source flux in erg/s/cm**2 on the first line, commented with #, and must have three
    columns of RA (degrees), Dec (degrees), and energy (keV) for each event.

    Example:

    # 1.194e-15
    30.1  45.5  2.71
    29.67 44.95 0.31
    31.25 45.03 10.01
    29.75 44.44 7.34
    30.05 44.01 12.01
    31.99 45.21 0.05
    ...

    positional arguments:
      simput_prefix  The prefix of the SIMPUT file to be used as the root of the catalog.
                     If it does not exist, it will be created.
      phlist_prefix  The prefix of the photon list file to be written.
      infile         The file containing the flux and positions and energies.

    optional arguments:
      -h, --help     show this help message and exit
      --append       If set, append a new source an existing SIMPUT catalog.
      --clobber      Whether or not to clobber an existing file with the same name.
