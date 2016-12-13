.. _event-tools:

Event File Tools
================

This section documents some helpful tools to take event files produced by the SOXS instrument
simulator and make derivative products from them. 

``write_image``
---------------

:func:`~soxs.events.write_image` bins up events into an image according to the coordinate
system inherent in the event file and writes the image to a FITS file. Images of sky, detector,
or chip coordinates can be written. You can also restrict events within a particular energy range 
to be written to the file.

To write an image in sky coordinates:

.. code-block:: python

    from soxs import write_image
    # Energy bounds are in keV
    write_image("my_evt.fits", "my_sky_img.fits", emin=0.5, emax=7.0)
    
Or in detector coordinates:

.. code-block:: python

    write_image("my_evt.fits", "my_sky_img.fits", coord_type='det', emin=0.5, emax=7.0)

Or in chip coordinates:

.. code-block:: python

    write_image("my_evt.fits", "my_sky_img.fits", coord_type='sky', emin=0.5, emax=7.0)

This image can then be viewed in `ds9 <http://ds9.si.edu>`_ or `APLpy <https://aplpy.github.io>`_.

``write_radial_profile``
------------------------

``write_spectrum``
------------------

:func:`~soxs.events.write_image` bins up events into a spectrum and writes the spectrum
to a FITS file:

.. code-block:: python

    from soxs import write_spectrum
    write_spectrum("my_evt.fits", "my_spec.pha", clobber=True)

This spectrum file can be read and fit with standard X-ray analysis software such as 
`XSPEC <https://heasarc.gsfc.nasa.gov/xanadu/xspec/>`_, `ISIS <http://space.mit.edu/CXC/ISIS/>`_, 
and `Sherpa <http://cxc.harvard.edu/sherpa/>`_. 