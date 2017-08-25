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
   dithering of their positions (if dithering is enabled for that particular
   instrument).
3. Add background events.
4. Convolves the event energies with the response matrix to produce channels.
5. Writes everything to an event file.

All of the photon lists in the SIMPUT catalog will be processed. A typical 
invocation of :func:`~soxs.instrument.instrument_simulator` looks like the 
following:

.. code-block:: python

    from soxs import instrument_simulator
    simput_file = "snr_simput.fits" # SIMPUT file to be read
    out_file = "evt_mucal.fits" # event file to be written
    exp_time = (30.0, "ks") # The exposure time
    instrument = "mucal" # short name for instrument to be used
    sky_center = [30., 45.] # RA, Dec of pointing in degrees
    instrument_simulator(simput_file, out_file, exp_time, instrument, 
                         sky_center, overwrite=True)
 
The ``overwrite`` argument allows an existing file to be overwritten.

.. _coords:

Coordinate Systems in SOXS
--------------------------

SOXS event files produced by the instrument simulator have two coordinate systems:
the (X, Y) "sky" coordinate system and the (DETX, DETY) "detector" coordinate system.

For a given instrument specification, the detector space is defined by the field of 
view parameter ``fov``, which is in arcminutes, and is divided into ``num_pixels``
pixels on a side. The field of view is shown in the schematic diagram in Figure 1 
as the dashed red square. The center of the field of view has detector coordinates 
0,0, as can be seen in Figure 1. 

The sky coordinate system is defined to be twice the size of the ``fov`` parameter, with
twice as many pixels. The center of the sky coordinate system is given by pixel 
coordinates ``0.5*(2*num_pixels+1),0.5*(2*num_+pixels+1)``. The sky coordinate system is
also shown in Figure 1. In event files and images, standard world coordinate system 
(WCS) keywords are used to translate between sky coordinates and RA and Dec. 

.. figure:: ../images/det_schematic.png
    :width: 700px

    Figure 1: Schematic showing the layout of sky and detector coordinate systems, 
    as well as multiple chips, for an example instrument similar to *Chandra*/ACIS-I. 
    A roll angle of 45 degrees has been specified. 

If the ``roll_angle`` parameter of the instrument simulation is 0, the sky and detector
coordinate systems will be aligned, but otherwise they will not. Figure 1 shows the 
orientation of the detector in the sky coordinate system for a roll angle of 45 degrees. 
For observations which have dither, the sky coordinates and the detector coordinates
will not have a one-to-one mapping, but will change as a function of time. 

Finally, Figure 1 also shows that multiple chips can be specified. In this case, only
events which fall within the chip regions are detected. For more information on how
multiple chips can be specified for a particlular instrument, see :ref:`chips`.

.. warning::

    At the present time, the coordinate systems specified in SOXS do not correspond 
    directly to those systems in event files produced by actual X-ray observatories.
    This is particularly true of detector coordinates. The conventions chosen by 
    SOXS are mainly for convenience. 

.. _instrument-arg:

The ``instrument`` Argument
+++++++++++++++++++++++++++

SOXS currently supports instrument configurations for *Lynx*, *Athena*, *Chandra*, 
and *Hitomi* "out of the box". Any of these can be specified with the ``instrument`` 
argument:

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
along with a 5-arcsecond Gaussian PSF, and observations are not dithered. The
WFI detector consists of four chips laid out in a 2x2 shape with a field of view
of approximately 40 arcminutes, and the X-IFU detector has a single hexagonal 
shape with an approximate diameter of 5 arcminutes. For more information about 
the specification of the *Athena* instruments assumed here, consult 
`the Athena simulation tools web portal <http://www.the-athena-x-ray-observatory.eu/resources/simulation-tools.html>`_.

Chandra
~~~~~~~

For simulating *Chandra* observations, a number of instrument specifications are 
available. All specifications assume a 10-meter focal length, 0.5-arcsecond Gaussian 
PSF, dithering, and 0.492-arcsecond pixels.
 
ACIS-I
######

The two ACIS-I specifications have a square field of view of roughly 20 
arcminutes, laid out in four chips 8 arcminutes on a side arranged 2x2. However,
The two separate specifications, ``"acisi_cy0"`` and ``"acisi_cy19"``, use the 
instrumental responses from shortly after launch ("Cycle 0") and from more 
recently ("Cycle 19"), respectively. The main effect is that the effective area
at low energies for ``"acisi_cy19"`` is much lower due to the buildup of 
contamination on the ACIS optical blocking filters compared to the 
``"acisi_cy0"`` responses.

ACIS-S
######

The two ACIS-S specifications have 6 chips 8 arcminutes on a side in a single row.
As in the ACIS-I case, the two specifications are for Cycle 0 ``"aciss_cy0"``, and
Cycle 19, ``"aciss_cy19"``. 

Gratings
########


Hitomi
~~~~~~

A single instrument specification is available for *Hitomi*, for the Soft X-ray
Spectrometer (SXS) instrument. It has a 5.6-meter focal length, a 1.2-arcminute
Gaussian PSF, no dithering, a 3-arcminute field of view, and 0.5-arcminute pixels.
The ARF and RMF for this specification were generated by Eric Miller (MIT). The 
RMF was produced for a single pixel response using the HEASOFT FTOOL ``sxsrmf``. 
The ARF was produced using ``aharfgen`` for the full field of view response to 
a point source observed on axis, assuming an empty filter position and the gate 
valve open. HEASOFT v6.20 FTOOLS were used, along with Hitomi CALDB v5 (release 
date 2016-12-23).

.. _bkgnds:

Backgrounds
+++++++++++

The instrument simulator simulates background events as well as the source
events provided by the user. There are three background components: the 
Galactic foreground, a background comprised of discrete point sources, and the 
instrumental/particle background. Complete information about these components 
can be found in :ref:`background`, but here the keyword arguments pertaining to
backgrounds for :func:`~soxs.instrument.instrument_simulator` will be detailed. 

The various background components can be turned on and off using 
the ``ptsrc_bkgnd``, ``instr_bkgnd``, and ``foreground`` arguments. They are all
on by default, but can be turned on or off individually:

.. code-block:: python

    # turns off the astrophysical background but leaves in the instrumental
    instrument_simulator(simput_file, out_file, exp_time, instrument, 
                         sky_center, overwrite=True, instr_bkgnd=False,
                         foreground=True) # ptsrc_bkgnd True by default

For long exposures, backgrounds may take a long time to generate. For this
reason, SOXS provides a way to add a background stored in a previously
generated event file to the simulation of a source, via the ``bkgnd_file``
argument:

.. code-block:: python

    # loads the background from a file
    instrument_simulator(simput_file, out_file, exp_time, instrument, 
                         sky_center, overwrite=True, bkgnd_file="my_bkgnd.fits") 

In this case the values of ``instr_bkgnd``, ``ptsrc_bkgnd``, and ``foreground``
are ignored regardless of their value. The required background event file can be
generated using :func:`~soxs.instrument.make_background_file`, and is documented
at :ref:`make-bkgnd`. The background event file must be for the same instrument 
as the one that is being simulated for the source and must have an exposure time
at least as long as the source exposure. 

.. _other-mods:

Other Modifications
+++++++++++++++++++

You can also change other aspects of the observation with 
:func:`~soxs.instrument.instrument_simulator`. For example, you can change the
size and period of the Lissajous dither pattern, for instruments which have 
dithering enabled. The default dither pattern has amplitudes of 8.0 arcseconds 
in the DETX and DETY directions, and a period of 1000.0 seconds in the DETX 
direction and a period of 707.0 seconds in the DETY direction. You can change
these numbers by supplying a list of parameters to the ``dither_params`` argument:

.. code-block:: python

    import soxs
    # The order of dither_params is [x_amp, y_amp, x_period, y_period]
    # the units of the amplitudes are in arcseconds and the periods are in
    # seconds
    dither_params = [8.0, 16.0, 1000.0, 2121.0]
    soxs.instrument_simulator(simput_file, out_file, exp_time, instrument, 
                              sky_center, overwrite=True, 
                              dither_params=dither_params)
    
To turn dithering off entirely for instruments that enable it, use the 
``no_dither`` argument:

.. code-block:: python

    import soxs
    soxs.instrument_simulator(simput_file, out_file, exp_time, instrument, 
                              sky_center, overwrite=True, 
                              no_dither=True)

.. note:: 

    Dithering will only be enabled if the instrument specification allows for 
    it. For example, for *Lynx*, dithering is on by default, but for *Athena* 
    it is off. 

.. _simulate-spectrum:

Simulating Spectra Only 
-----------------------

If you would like to use an instrument specification and a 
:class:`~soxs.spectra.Spectrum` object to generate a spectrum file only (without
including spatial effects), SOXS provides a function 
:func:`~soxs.instrument.simulate_spectrum` which can take an unconvolved
spectrum and generate a convolved one from it. This is similar to what the XSPEC
command "fakeit" does. 

.. code-block:: python

    spec = soxs.Spectrum.from_file("lots_of_lines.dat")
    instrument = "mucal"
    out_file = "lots_of_lines.pha"
    simulate_spectrum(spec, instrument, exp_time, out_file, overwrite=True)

This spectrum file then can be read in and analyzed by standard software such as
XSPEC, Sherpa, ISIS, etc.

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
        name: hdxi_3x10
        arf: xrs_hdxi_3x10.arf
        rmf: xrs_hdxi.rmf
        bkgnd: acisi
        fov: 20.0
        num_pixels: 4096
        aimpt_coords: [0.0, 0.0]
        chips: None
        focal_length: 10.0
        dither: True
        psf: ['gaussian', 0.5]
    Instrument: athena_xifu
        name: athena_xifu
        arf: athena_xifu_1469_onaxis_pitch249um_v20160401.arf
        rmf: athena_xifu_rmf_v20160401.rmf
        bkgnd: athena_xifu
        fov: 5.991992621478149
        num_pixels: 84
        aimpt_coords: [0.0, 0.0]
        chips: [['Polygon', 
                 [-33, 0, 33, 33, 0, -33], 
                 [20, 38, 20, -20, -38, -20]]]
        focal_length: 12.0
        dither: False
        psf: ['gaussian', 5.0]
    Instrument: acisi_cy19
        name: acisi_cy19
        arf: acisi_aimpt_cy19.arf
        rmf: acisi_aimpt_cy19.rmf
        bkgnd: acisi
        fov: 20.008
        num_pixels: 2440
        aimpt_coords: [86.0, 57.0]
        chips: [['Box', -523, -523, 1024, 1024], 
                ['Box', 523, -523, 1024, 1024], 
                ['Box', -523, 523, 1024, 1024], 
                ['Box', 523, 523, 1024, 1024]]
        psf: ['gaussian', 0.5]
        focal_length: 10.0
        dither: True
    Instrument: hitomi_sxs
        name: hitomi_sxs
        arf: hitomi_sxs_ptsrc.arf
        rmf: hitomi_sxs.rmf
        bkgnd: hitomi_sxs
        num_pixels: 6
        fov: 3.06450576
        aimpt_coords: [0.0, 0.0]
        chips: None
        focal_length: 5.6
        dither: False
        psf: ['gaussian', 72.0]
    ...

The various parts of each instrument specification are:

* ``"name"``: The name of the instrument specification. 
* ``"arf"``: The file containing the ARF.
* ``"rmf"``: The file containing the RMF.
* ``"fov"``: The field of view in arcminutes. This may represent a single chip
  or an area within which chips are embedded.
* ``"num_pixels"``: The number of resolution elements on a side of the field of 
  view.
* ``"chips"``: The specification for multiple chips, if desired. For more details
  on how to specify chips, see :ref:`chips`. 
* ``"bkgnd"``: The name of the instrumental background to use, stored in the 
  background registry (see :ref:`background` for more details). This can also be
  set to ``None`` for no particle background.
* ``"psf"``: The PSF specification to use. At time of writing, the only one 
  available is that of a Gaussian PSF, with a single parameter, the HPD of the 
  PSF. This is specified using a Python list, e.g. ``["gaussian", 0.5]``. This 
  can also be set to ``None`` for no PSF.
* ``"focal_length"``: The focal length of the telescope in meters.
* ``"dither"``: Whether or not the instrument dithers by default. 

As SOXS matures, this list of specifications will likely expand, and the number 
of options for some of them (e.g., the PSF) will also expand.

.. _custom-instruments:

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
    2. ``True`` and ``False`` are capitalized in Python, in JSON they are lowercase.

.. _chips:

Defining Instruments with Multiple Chips
++++++++++++++++++++++++++++++++++++++++

If the ``"chips"`` entry in the instrument specification is ``None``, then there
will only be one chip which covers the entire field of view. However, it is also 
possible to specify multiple chips with essentially arbitary shapes. In this case, 
the ``"chips"`` entry needs to be a list containing a set of lists, one for each
chip, that specifies a region expression parseable by the 
`pyregion <https://pyregion.readthedocs.io>`_ package. 

Three options are currently recognized by SOXS for chip shapes:

* Rectangle shapes, which use the ``Box`` region. The four arguments are ``xc``
  (center in the x-coordinate), ``yc`` (center in the y-coordinate), ``width``,
  and ``height``.
* Circle shapes, which use the ``Circle`` region. The three arguments are ``xc``
  (center in the x-coordinate), ``yc`` (center in the y-coordinate), and ``radius``.
* Generic polygon shapes, which use the ``Polygon`` region. The two arguments are
  ``x`` and ``y``, which are lists of x and y coordinates for each point of the
  polygon. 

To create a chip, simply supply a list starting with the name of the region 
type and followed by the arguments in order. All coordinates and distances are
in detector coordinates. For example, a ``Box`` region at detector coordinates
(0,0) with a width of 100 pixels and a height of 200 pixels would be specified
as ``["Box", 0.0, 0.0, 100, 200]``. 

For example, the *Chandra* ACIS-I instrument configurations have a list of four 
``Box`` regions to specify the four I-array square-shaped chips:

.. code-block:: python

    instrument_registry["acisi_cy18"] = {"name": "acisi_cy18",
                                         "arf": "acisi_aimpt_cy18.arf",
                                         "rmf": "acisi_aimpt_cy18.rmf",
                                         "bkgnd": "acisi",
                                         "fov": 20.008,
                                         "num_pixels": 2440,
                                         "aimpt_coords": [86.0, 57.0],
                                         "chips": [["Box", -523, -523, 1024, 1024],
                                                   ["Box", 523, -523, 1024, 1024],
                                                   ["Box", -523, 523, 1024, 1024],
                                                   ["Box", 523, 523, 1024, 1024]],
                                         "psf": ["gaussian", 0.5],
                                         "focal_length": 10.0,
                                         "dither": True}

whereas the *Athena* XIFU instrument configuration uses a ``Polygon`` region:

.. code-block:: python

    instrument_registry["athena_xifu"] = {"name": "athena_xifu",
                                          "arf": "athena_xifu_1469_onaxis_pitch249um_v20160401.arf",
                                          "rmf": "athena_xifu_rmf_v20160401.rmf",
                                          "bkgnd": "athena_xifu",
                                          "fov": 5.991992621478149,
                                          "num_pixels": 84,
                                          "aimpt_coords": [0.0, 0.0],
                                          "chips": [["Polygon", 
                                                     [-33, 0, 33, 33, 0, -33],
                                                     [20, 38, 20, -20, -38, -20]]],
                                          "focal_length": 12.0,
                                          "dither": False,
                                          "psf": ["gaussian", 5.0]}
