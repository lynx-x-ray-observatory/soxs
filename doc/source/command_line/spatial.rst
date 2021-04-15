.. _cmd-spatial:

Command Line Scripts for Spatial Models
=======================================

These are scripts that create photon lists for SIMPUT catalogs that can serve
as inputs to the instrument simulator. For details on what's going on under the 
hood, see :ref:`spatial`, :ref:`spectra`, and :ref:`simput`.

Each of these scripts accepts a ``specfile`` argument, which is a two-column ASCII
table of energy (keV) and specific flux (:math:`\rm{photons~s^{-1}~cm^{-2}~keV^{-1}}`) 
which can be written using the spectra command-line scripts (see :ref:`cmd-spectra`), 
the Python method :meth:`~soxs.spectra.Spectrum.write_file` (see :ref:`spectra`), 
or can be created by hand.

.. _cmd-make-point-source:

``make_point_source``
---------------------

This script creates a SIMPUT source of a point source from a spectrum supplied in a
file.

.. code-block:: text

    usage: make_point_source [-h] [--src_filename SRC_FILENAME] [--append] [--overwrite] filename name ra0 dec0 specfile
    
    Create a SIMPUT source of a point from a spectrum supplied in a file.
    
    positional arguments:
      filename              The filename of the SIMPUT catalog to write or to append to.
      name                  The name of the source in the SIMPUT catalog.
      ra0                   The right ascension of the source in degrees.
      dec0                  The declination of the source in degrees.
      specfile              The file containing the spectrum to be used.
    
    optional arguments:
      -h, --help            show this help message and exit
      --src_filename SRC_FILENAME
                            An optional filename to store the source instead of the SIMPUT catalog file.
      --append              If set, append a new source an existing SIMPUT catalog.
      --overwrite           Overwrite an existing file with the same name.

Examples
++++++++

Make a brand-new SIMPUT catalog for a point source. 

.. code-block:: bash

    [~]$ make_point_source pt_src.simput src1 20.0 -32.0 pt_src_spectrum.dat --overwrite

Add a new point source to an existing SIMPUT catalog. 

.. code-block:: bash

    [~]$ make_point_source pt_src.simput src2 19.0 -31.0 pt_src_spectrum.dat --append --overwrite

Make a brand-new SIMPUT catalog for a point source, but write the source to a different file. 

.. code-block:: bash

    [~]$ make_point_source pt_src.simput src2 19.0 -31.0 pt_src_spectrum.dat --src_filename=pt.fits --overwrite

.. _cmd-make-beta-model-source:

``make_beta_model_source``
--------------------------

This script creates a SIMPUT source of a :math:`\beta`-model from a spectrum supplied in a
file. The functional form of the :math:`\beta`-model for a surface brightness profile is:

.. math::

    S(r) = S_0\left[1+\left(\frac{r}{r_c}\right)^2\right]^{(-3\beta+1/2)}

.. code-block:: text

    usage: make_beta_model_source [-h] [--theta THETA] [--ellipticity ELLIPTICITY] [--src_filename SRC_FILENAME] [--append] [--overwrite]
                                  filename name ra0 dec0 r_c beta specfile image_width nx
    
    Create a SIMPUT source of a beta-model from a spectrum supplied in a file.
    
    positional arguments:
      filename              The filename of the SIMPUT catalog to write or to append to.
      name                  The name of the source in the SIMPUT catalog.
      ra0                   The right ascension of the source center in degrees.
      dec0                  The declination of the source center in degrees.
      r_c                   The core radius in arcseconds.
      beta                  The beta parameter.
      specfile              The file containing the spectrum to be used.
      image_width           The width of the image in arcminutes.
      nx                    The resolution of the image.
    
    optional arguments:
      -h, --help            show this help message and exit
      --theta THETA         The angle through which to rotate the beta model in degrees. Only makes sense if ellipticity is added. Default:
                            0.0
      --ellipticity ELLIPTICITY
                            The ellipticity of the radial profile, expressed as the ratio between the length scales of the x and y
                            coordinates. The value of this parameter will shrink or expand the profile in the direction of the "y"
                            coordinate, so you may need to rotate to get the shape you want. Default: 1.0
      --src_filename SRC_FILENAME
                            An optional filename to store the source instead of the SIMPUT catalog file.
      --append              If set, append a new source an existing SIMPUT catalog.
      --overwrite           Overwrite an existing file with the same name.

Examples
++++++++

Make a brand-new SIMPUT catalog for a :math:`\beta`-model source. 

.. code-block:: bash

    [~]$ make_beta_model_source my_srcs.simput beta_src1 20.0 -32.0 10.0 1.0 my_spectrum.dat 30.0 2000 --overwrite

Add a new :math:`\beta`-model to an existing SIMPUT catalog. 

.. code-block:: bash

    [~]$ make_beta_model_source my_srcs.simput beta_src2 19.0 -31.0 10.0 1.0 my_spectrum.dat 30.0 2000 --append --overwrite

Make a brand-new SIMPUT catalog for a :math:`\beta`-model source, but write the source to a different file. 

.. code-block:: bash

    [~]$ make_beta_model_source my_srcs.simput src2 19.0 -31.0 10.0 1.0 my_spectrum.dat 30.0 2000 --src_filename=beta.fits --overwrite

Add a new :math:`\beta`-model to an existing SIMPUT catalog, but write the source to
a different file.

.. code-block:: bash

    [~]$ make_beta_model_source my_srcs.simput beta_src2 19.0 -31.0 10.0 1.0 my_spectrum.dat 30.0 2000 --append --overwrite --src_filename=pt.fits

Add ellipticity and tilt the model:

.. code-block:: bash

    [~]$ make_beta_model_source my_srcs.simput beta_src1 20.0 -32.0 10.0 1.0 my_spectrum.dat 30.0 2000 --ellipticity=0.5 --theta=45.0 --overwrite

.. _cmd-make-double-beta-model-source:

``make_double_beta_model_source``
---------------------------------

This script creates a SIMPUT source of a double-:math:`\beta`-model from a spectrum 
supplied in a file. The functional form of the double-:math:`\beta`-model for a 
surface brightness profile is:

.. math::

    S(r) = S_{0,1}\left\{\left[1+\left(\frac{r}{r_{c,1}}\right)^2\right]^{(-3\beta_1+1/2)} +
           \frac{S_{0,2}}{S_{0,1}}\left[1+\left(\frac{r}{r_{c,2}}\right)^2\right]^{(-3\beta_2+1/2)}\right\}

.. code-block:: text

    usage: make_double_beta_model_source [-h] [--theta THETA] [--ellipticity ELLIPTICITY] [--src_filename SRC_FILENAME] [--append] [--overwrite]
                                         filename name ra0 dec0 r_c1 beta1 r_c2 beta2 sb_ratio specfile image_width nx
    
    Create a SIMPUT source of a double-beta-model from a spectrum supplied in a file.
    
    positional arguments:
      filename              The filename of the SIMPUT catalog to write or to append to.
      name                  The name of the source in the SIMPUT catalog.
      ra0                   The right ascension of the source center in degrees.
      dec0                  The declination of the source center in degrees.
      r_c1                  The inner core radius in arcseconds.
      beta1                 The inner beta parameter.
      r_c2                  The outer core radius in arcseconds.
      beta2                 The outer beta parameter.
      sb_ratio              The ratio of the outer to the inner SB peak value.
      specfile              The file containing the spectrum to be used.
      image_width           The width of the image in arcminutes.
      nx                    The resolution of the image.
    
    optional arguments:
      -h, --help            show this help message and exit
      --theta THETA         The angle through which to rotate the beta model in degrees. Only makes sense if ellipticity is added. Default:
                            0.0
      --ellipticity ELLIPTICITY
                            The ellipticity of the radial profile, expressed as the ratio between the length scales of the x and y
                            coordinates. The value of this parameter will shrink or expand the profile in the direction of the "y"
                            coordinate, so you may need to rotate to get the shape you want. Default: 1.0
      --src_filename SRC_FILENAME
                            An optional filename to store the source instead of the SIMPUT catalog file.
      --append              If set, append a new source an existing SIMPUT catalog.
      --overwrite           Overwrite an existing file with the same name.

Examples
++++++++

Make a brand-new SIMPUT catalog for a double-:math:`\beta`-model source. 

.. code-block:: bash

    [~]$ make_double_beta_model_source my_srcs.simput beta_src1 20.0 -32.0 10.0 1.0 30.0 0.6666 0.5 my_spectrum.dat 30.0 2000 --overwrite

Add a new double-:math:`\beta`-model to an existing SIMPUT catalog. 

.. code-block:: bash

    [~]$ make_double_beta_model_source my_srcs.simput beta_src2 19.0 -31.0 10.0 1.0 30.0 0.6666 0.5 my_spectrum.dat 30.0 2000 --append --overwrite

Make a brand-new SIMPUT catalog for a double-:math:`\beta`-model, but write the source to a different file. 

.. code-block:: bash

    [~]$ make_double_beta_model_source my_srcs.simput src2 19.0 -31.0 10.0 1.0 30.0 0.6666 0.5 my_spectrum.dat 30.0 2000 --src_filename=pt.fits --overwrite

Add a new double-:math:`\beta`-model to an existing SIMPUT catalog, but write the source to
a different file.

.. code-block:: bash

    [~]$ make_double_beta_model_source my_srcs.simput beta_src2 19.0 -31.0 10.0 1.0 30.0 0.6666 0.5 my_spectrum.dat 30.0 2000 --append --overwrite --src_filename=pt.fits

Add ellipticity and tilt the model:

.. code-block:: bash

    [~]$ make_double_beta_model_source my_srcs.simput beta_src1 20.0 -32.0 10.0 1.0 30.0 0.6666 0.5 my_spectrum.dat 30.0 2000 --ellipticity=0.5 --theta=45.0 --overwrite

.. _cmd-make-annulus-source:

``make_annulus_source``
-----------------------

This script creates a SIMPUT source of an annulus or disk with constant surface brightness
from a spectrum supplied in a file.

.. code-block:: text

    usage: make_annulus_source [-h] [--theta THETA] [--ellipticity ELLIPTICITY] [--src_filename SRC_FILENAME] [--append] [--overwrite]
                               filename name ra0 dec0 r_in r_out specfile image_width nx
    
    Create a SIMPUT source of an annulus with uniform surface brightness from a spectrum supplied in a file.
    
    positional arguments:
      filename              The filename of the SIMPUT catalog to write or to append to.
      name                  The name of the source in the SIMPUT catalog.
      ra0                   The right ascension of the source center in degrees.
      dec0                  The declination of the source center in degrees.
      r_in                  The inner annulus of the source center in arcseconds.
      r_out                 The outer annulus of the source center in arcseconds.
      specfile              The file containing the spectrum to be used.
      image_width           The width of the image in arcminutes.
      nx                    The resolution of the image.
    
    optional arguments:
      -h, --help            show this help message and exit
      --theta THETA         The angle through which to rotate the beta model in degrees. Only makes sense if ellipticity is added. Default:
                            0.0
      --ellipticity ELLIPTICITY
                            The ellipticity of the radial profile, expressed as the ratio between the length scales of the x and y
                            coordinates. The value of this parameter will shrink or expand the profile in the direction of the "y"
                            coordinate, so you may need to rotate to get the shape you want. Default: 1.0
      --src_filename SRC_FILENAME
                            An optional filename to store the source instead of the SIMPUT catalog file.
      --append              If set, append a new source an existing SIMPUT catalog.
      --overwrite           Overwrite an existing file with the same name.

Examples
++++++++

Make a brand-new SIMPUT catalog for an annulus source. 

.. code-block:: bash

    [~]$ make_annulus_source my_srcs.simput ann_src1 20.0 -32.0 0.0 30.0 my_spectrum.dat 30.0 2000 --overwrite

Add a new annulus model to an existing SIMPUT catalog. 

.. code-block:: bash

    [~]$ make_annulus_source my_srcs.simput ann_src2 19.0 -31.0 0.0 30.0 my_spectrum.dat 30.0 2000 --append --overwrite

Add ellipticity and tilt the model:

.. code-block:: bash

    [~]$ make_annulus_source my_srcs.simput ann_src1 20.0 -32.0 0.0 30.0 my_spectrum.dat 30.0 2000 --ellipticity=2.0 --theta=30.0 --overwrite

.. _cmd-make-rectangle-source:

``make_rectangle_source``
-------------------------

This script creates a SIMPUT source of a rectangle shape with constant surface brightness
from a spectrum supplied in a file.

.. code-block:: text

    usage: make_rectangle_source [-h] [--theta THETA] [--src_filename SRC_FILENAME] [--append] [--overwrite]
                                 filename name ra0 dec0 width height specfile image_width nx
    
    Create a SIMPUT source of a uniformly filled rectangle from a spectrum supplied in a file.
    
    positional arguments:
      filename              The filename of the SIMPUT catalog to write or to append to.
      name                  The name of the source in the SIMPUT catalog.
      ra0                   The right ascension of the source center in degrees.
      dec0                  The declination of the source center in degrees.
      width                 The width of the rectangle in arcseconds.
      height                The width of the rectangle in arcseconds.
      specfile              The file containing the spectrum to be used.
      image_width           The width of the image in arcminutes.
      nx                    The resolution of the image.
    
    optional arguments:
      -h, --help            show this help message and exit
      --theta THETA         The angle through which to rotate the rectangle in degrees. Default: 0.0
      --src_filename SRC_FILENAME
                            An optional filename to store the source instead of the SIMPUT catalog file.
      --append              If set, append a new source an existing SIMPUT catalog.
      --overwrite           Overwrite an existing file with the same name.

Examples
++++++++

Make a brand-new SIMPUT catalog for a rectangle source.

.. code-block:: bash

    [~]$ make_rectangle_source my_srcs.simput rect_src1 20.0 -32.0 20.0 10.0 my_spectrum.dat 30.0 2000 --overwrite

Make the same rectangle, but rotate it by 30.0 degrees.

.. code-block:: bash

    [~]$ make_rectangle_source my_srcs.simput rect_src1 20.0 -32.0 20.0 10.0 my_spectrum.dat 30.0 2000 --theta=30.0 --overwrite

Create a line source with the same width and rotation angle.

.. code-block:: bash

    [~]$ make_rectangle_source my_srcs.simput rect_src1 20.0 -32.0 20.0 0.0 my_spectrum.dat 30.0 2000 --theta=30.0 --overwrite

Add a new rectangle model to an existing SIMPUT catalog.

.. code-block:: bash

    [~]$ make_rectangle_source my_srcs.simput rect_src2 19.0 -31.0 20.0 10.0 my_spectrum.dat 30.0 2000 --append --overwrite

.. _cmd-make-fov-source:

``make_fov_source``
-------------------

This script creates a SIMPUT source of a field of view with constant surface brightness
from a spectrum supplied in a file.

.. code-block:: text

    usage: make_fov_source [-h] [--src_filename SRC_FILENAME] [--append] [--overwrite] filename name ra0 dec0 fov specfile image_width nx
    
    Create a SIMPUT source of a uniformly filled field of view from a spectrum supplied in a file.
    
    positional arguments:
      filename              The filename of the SIMPUT catalog to write or to append to.
      name                  The name of the source in the SIMPUT catalog.
      ra0                   The right ascension of the source center in degrees.
      dec0                  The declination of the source center in degrees.
      fov                   The field of view on a side in arcminutes.
      specfile              The file containing the spectrum to be used.
      image_width           The width of the image in arcminutes.
      nx                    The resolution of the image.
    
    optional arguments:
      -h, --help            show this help message and exit
      --src_filename SRC_FILENAME
                            An optional filename to store the source instead of the SIMPUT catalog file.
      --append              If set, append a new source an existing SIMPUT catalog.
      --overwrite           Overwrite an existing file with the same name.

Examples
++++++++

Make a brand-new SIMPUT catalog for a field-of-view source. 

.. code-block:: bash

    [~]$ make_fov_source my_srcs.simput fov_src1 20.0 -32.0 20.0 my_spectrum.dat 30.0 2000 --overwrite

Add a new field-of-view model to an existing SIMPUT catalog. 

.. code-block:: bash

    [~]$ make_fov_source my_srcs.simput fov_src2 19.0 -31.0 20.0 my_spectrum.dat 30.0 2000 --append --overwrite

``make_phlist_from_ascii``
--------------------------

This script takes a table of photon RA, Dec, and energies from an ASCII-formatted table and writes them
to a new SIMPUT catalog with a photon list. 

.. code-block:: text

    usage: make_phlist_from_ascii [-h] [--src_filename SRC_FILENAME] [--append] [--overwrite] name filename infile
    
    Create a SIMPUT source from an ASCII table of positions and energies. The file must contain the total source flux in erg/s/cm**2 on the first line, commented with #, and must have three columns of RA (degrees), Dec (degrees), and energy (keV) for each event.
    
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
      filename              The filename of the SIMPUT catalog to write or to append to.
      name                  The name of the source in the SIMPUT catalog.
      infile                The file containing the flux and positions and energies.
    
    optional arguments:
      -h, --help            show this help message and exit
      --src_filename SRC_FILENAME
                            An optional filename to store the source instead of the SIMPUT catalog file.
      --append              If set, append a new source an existing SIMPUT catalog.
      --overwrite           Overwrite an existing file with the same name.

Examples
++++++++

Read photons from a file and write to a new SIMPUT catalog file. 

.. code-block:: bash

    [~]$ make_phlist_from_ascii my_cat.simput photons events.txt --overwrite

Read photons from a file and write to a new SIMPUT catalog file, but write the photons to a new file.

.. code-block:: bash

    [~]$ make_phlist_from_ascii my_cat.simput photons events.txt --src_filename=photons.fits

Read photons from a file and append to an existing SIMPUT catalog file.

.. code-block:: bash

    [~]$ make_phlist_from_ascii my_cat.simput photons events.txt --append


