.. _instrument:

The Instrument Simulator
========================

Running the Instrument Simulator
--------------------------------

The end product of a mock observation is a "standard" event file which has been 
convolved with a model for the telescope. In SOXS, this is handled by the
instrument simulator. 

:func:`~soxs.instrument.instrument_simulator` reads in a SIMPUT catalog and 
creates a standard event file using the instrument simulator. 
:func:`~soxs.instrument.instrument_simulator` performs the following actions:

1. Uses the effective area curve to determine which events will actually be 
   detected.
2. Projects these events onto the detector plane and perform PSF blurring and 
   dithering of their positions.
3. Add particle/instrumental and astrophysical background events.
4. Convolves the event energies with the response matrix to produce channels.
5. Writes everything to an event file.

All of the photon lists in the SIMPUT catalog will be processed. A typical 
invocation of :func:`~soxs.instrument.instrument_simulator` looks like the 
following:

.. code-block:: python

    from soxs import instrument_simulator
    simput_file = "snr_simput.fits" # SIMPUT file to be read
    out_file = "evt_mucal.fits" # event file to be written
    exp_time = 30000. # The exposure time in seconds
    instrument = "mucal" # short name for instrument to be used
    sky_center = [30., 45.] # RA, Dec of pointing in degrees
    instrument_simulator(simput_file, out_file, exp_time, instrument, 
                         sky_center, clobber=True)
 
The ``clobber`` argument allows an existing file to be overwritten.

.. _instrument-arg:

The ``instrument`` Argument
+++++++++++++++++++++++++++

SOXS currently supports instrument configurations for *Lynx* and *Athena* "out 
of the box". Any of these can be specified with the ``instrument`` argument:

Lynx
~~~~

For *Lynx*, there are currently two base instruments, ``"hdxi"`` for the 
High-Definition X-ray Imager, and ``"mucal"`` for the microcalorimeter. There 
are also variations on these instruments which use different mirror parameters. 
The different variations on mirror parameters available are:

* :math:`d` = 3 m, :math:`f` = 10 m (default case)
* :math:`d` = 3 m, :math:`f` = 15 m
* :math:`d` = 3 m, :math:`f` = 20 m
* :math:`d` = 6 m, :math:`f` = 20 m

Where :math:`d` is the diameter of the outermost mirror shell, and :math:`f` is
the focal length. To use a different case other than the default, append it to 
the instrument string in a ``dxf`` pattern, e.g. ``"hdxi_3x20"``, 
``"mucal_6x20"``.

Athena
~~~~~~

For simulating *Athena* observations, two instrument specifications are 
available, for the WFI (Wide-Field Imager) and the X-IFU (X-ray Integral Field 
Unit). For both of these specifications, a 12-meter focal length is assumed, 
along with a 5-arcsecond Gaussian PSF, and observations are not dithered. For 
more information about the specification of the *Athena* instruments assumed 
here, consult 
`the Athena simulation tools web portal <http://www.the-athena-x-ray-observatory.eu/resources/simulation-tools.html>`_.

Chandra
~~~~~~~

For simulating *Chandra* observations, two instrument specifications are 
available, both for the ACIS-I instrument. These specifications are almost 
identical with a 10-meter focal length, 0.5-arcsecond Gaussian PSF, dithering, 
0.492-arcsecond pixels, and roughly 16.9 arcminute field of view. However, The 
two separate specifications, ``"acisi_cy0"`` and ``"acisi_cy18"``, use the 
instrumental responses from shortly after launch ("Cycle 0") and from more 
recently ("Cycle 18"), respectively. The main effect is that the effective area 
at low energies for ``"acisi_cy18"`` is much lower due to the buildup of 
contamination on the ACIS optical blocking filters compared to the 
``"acisi_cy0"`` responses.

.. _other-mods:

Other Modifications
+++++++++++++++++++

You can also change other aspects of the observation with 
:func:`~soxs.instrument.instrument_simulator`. For example, you can change the
shape and size of the dither pattern. The default dither pattern is a square of
width 16.0 arcseconds on a side. You can change it to be a circle dither pattern
or turn off dithering entirely:

.. code-block:: python

    # this invocation makes the dither shape a circle and 
    # sets the radius to 8 arcsec
    instrument_simulator(simput_file, out_file, exp_time, instrument, 
                         sky_center, clobber=True, dither_shape="circle", 
                         dither_size=8.0)
    
.. code-block:: python

    # this invocation turns off dithering entirely
    instrument_simulator(simput_file, out_file, exp_time, instrument, 
                         sky_center, clobber=True, dither_shape=None) 

You can also specify a non-zero roll angle:

.. code-block:: python

    # adds a roll of 45.0 degrees
    instrument_simulator(simput_file, out_file, exp_time, instrument, 
                         sky_center, clobber=True, roll_angle=45.0) 

.. note:: 

    Dithering will only be enabled if the instrument specification allows for 
    it. For *Lynx*, dithering is on by default, but for *Athena* it is off. 

The astrophysical and instrumental backgrounds can be turned on and off using 
the ``astro_bkgnd`` and ``instr_bkgnd`` arguments:

.. code-block:: python

    # decreases the particle background intensity by half
    instrument_simulator(simput_file, out_file, exp_time, instrument, 
                         sky_center, clobber=True, instr_bkgnd=False,
                         astro_bkgnd=True) 

Finally, to simulate an observation of backgrounds only without a source, simply 
run :func:`~soxs.instrument.instrument_simulator` with ``None`` as the argument 
for the ``simput_file``:

.. code-block:: python

    # simulates backgrounds only
    instrument_simulator(None, "bkg_evt.fits", exp_time, instrument, 
                         sky_center, clobber=True)

.. _instrument-registry:

Creating New Instrument Specifications
--------------------------------------

SOXS provides the ability to customize the models of the different components of
the instrument being simulated. This is provided by the use of the instrument 
registry and JSON files which contain prescriptions for different instrument 
configurations.

The Instrument Registry
+++++++++++++++++++++++

The instrument registry is simply a Python dictionary containing various 
instrument specifications. You can see the contents of the instrument registry 
by calling :func:`~soxs.instrument.show_instrument_registry`:

.. code-block:: python

    import soxs
    soxs.show_instrument_registry()

gives (showing only a subset for brevity):

.. code-block:: pycon

    Instrument: hdxi
        num_pixels: 4096
        fov: 5.0
        bkgnd: acisi
        psf: ['gaussian', 0.5]
        name: hdxi_3x10
        arf: xrs_hdxi_3x10.arf
        rmf: xrs_hdxi.rmf
        focal_length: 10.0
        dither: True
    Instrument: mucal
        num_pixels: 300
        fov: 5.0
        bkgnd: mucal
        psf: ['gaussian', 0.5]
        name: mucal_3x10
        arf: xrs_mucal_3x10.arf
        rmf: xrs_mucal.rmf
        focal_length: 10.0
        dither: True
    Instrument: mucal_3x15
        num_pixels: 300
        fov: 5.0
        bkgnd: mucal
        psf: ['gaussian', 0.5]
        name: mucal_3x15
        arf: xrs_mucal_3x15.arf
        rmf: xrs_mucal.rmf
        focal_length: 15.0
        dither: True
    Instrument: hdxi_3x15
        num_pixels: 4096
        fov: 20.0
        bkgnd: acisi
        psf: ['gaussian', 0.5]
        name: hdxi_3x15
        arf: xrs_hdxi_3x15.arf
        rmf: xrs_hdxi.rmf
        focal_length: 15.0
        dither: True
    Instrument: hdxi_3x10
        num_pixels: 4096
        fov: 20.0
        bkgnd: acisi
        psf: ['gaussian', 0.5]
        name: hdxi_3x10
        arf: xrs_hdxi_3x10.arf
        rmf: xrs_hdxi.rmf
        focal_length: 10.0
        dither: True
    ...

The various parts of each instrument specification are:

* ``"name"``: The name of the instrument specification. 
* ``"arf"``: The file containing the ARF.
* ``"num_pixels"``: The number of resolution elements on a side of the field of 
  view.
* ``"bkgnd"``: The name of the instrumental background to use, stored in the 
  background registry (see :ref:`background` for more details). This can also be
  set to ``None`` for no particle background.
* ``"psf"``: The PSF specification to use. At time of writing, the only one 
  available is that of a Gaussian PSF, with a single parameter, the HPD of the 
  PSF. This is specified using a Python list, e.g. ``["gaussian", 0.5]``. This 
  can also be set to ``None`` for no PSF.
* ``"rmf"``: The file containing the RMF.
* ``"fov"``: The field of view in arcminutes. 
* ``"focal_length"``: The focal length of the telescope in meters.
* ``"dither"``: Whether or not the instrument dithers by default. 

As SOXS matures, this list of specifications will likely expand, and the number 
of options for some of them (e.g., the PSF) will also expand.

Making Custom Instruments
+++++++++++++++++++++++++

To make a custom instrument, you can take an existing instrument specification 
and modify it, giving it a new name, or write a new specification to a 
`JSON <http://www.json.org>`_ file and read it in. To make a new specification 
from a dictionary, construct the dictionary and feed it to 
:func:`~soxs.instrument.add_instrument_to_registry`. For example, if you wanted 
to take the default calorimeter specification and change the plate scale, you 
would do it this way, using :func:`~soxs.instrument.get_instrument_from_registry`
to get the specification so that you can alter it:

.. code-block:: python

    from soxs import get_instrument_from_registry, add_instrument_to_registry
    new_mucal = get_instrument_from_registry("mucal")
    new_mucal["name"] = "mucal_high_res" # Must change the name, otherwise an error will be thrown
    new_mucal["num_pixels"] = 12000 # Results in an ambitiously smaller plate scale, 0.1 arcsec per pixel
    name = add_instrument_to_registry(new_mucal)
    
You can also store an instrument specification in a JSON file and import it:

.. code-block:: python

    name = add_instrument_to_registry("my_mucal.json")
    
You can download an example instrument specification JSON file 
`here <../example_mucal_spec.json>`_. 

You can also take an existing instrument specification and write it to a JSON 
file for editing using :func:`~soxs.instrument.write_instrument_json`:

.. code-block:: python

    from soxs import write_instrument_json
    # Using the "new_mucal" from above
    write_instrument_json("mucal_high_res", "mucal_high_res.json")

.. warning::

    Since JSON files use Javascript-style notation instead of Python's, there 
    are two differences one must note when creating JSON-based instrument 
    specifications:
    1. Python's ``None`` will convert to ``null``, and vice-versa.
    2. ``True`` and ``False`` are capitalized in Python, in JSON they are 
       lowercase.