.. _background:

Simulating Background in SOXS
=============================

SOXS contains support for including models of astrophysical and 
instrumental/particle backgrounds (or foregrounds) in a simulation of an 
observation. 

Backgrounds are turned on in the instrument simulator by default. Either the
astrophysical background or the instrumental background can be turned off
entirely in the call to :func:`~soxs.instrument.instrument_simulator` by setting
``astro_bkgnd`` and/or ``instr_bkgnd`` to ``True`` or ``False``:

.. code-block:: python

    # turn off the astrophysical background
    instrument_simulator(simput_file, out_file, exp_time, instrument, 
                         sky_center, clobber=True, astro_bkgnd=False)

.. code-block:: python

    # turn off the instrumental background
    instrument_simulator(simput_file, out_file, exp_time, instrument, 
                         sky_center, clobber=True, instr_bkgnd=False)

.. code-block:: python

    # turn off both backgrounds
    instrument_simulator(simput_file, out_file, exp_time, instrument, 
                         sky_center, clobber=True, astro_bkgnd=False,
                         instr_bkgnd=False)

Astrophysical Background
------------------------

The astrophysical background in SOXS is comprised of three components: a 
galactic foreground, a point-source background, and a cosmological background.
The astrophysical background is not tied to any particular instrument 
specification.

Galactic Foreground Model
+++++++++++++++++++++++++

The galactic foreground component is modeled as a sum of two thermal models, 
``apec+apec``, with parameters:

``apec`` **Model 1**

* ``kT``: :math:`\rm{0.2~keV}`
* ``abund``: :math:`\rm{1.0~Z_\odot}`
* ``redshift``: :math:`0.0`
* ``norm``: :math:`\rm{6.82251 \times 10^{-7}~10^{-14}\frac{\int{n_en_HdV}}{4{\pi}D_A(1+z)^2}}`
 
``apec`` **Model 2**

* ``kT``: :math:`\rm{0.099~keV}`
* ``abund``: :math:`\rm{1.0~Z_\odot}`
* ``redshift``: :math:`0.0`
* ``norm``: :math:`\rm{1.12328 \times 10^{-6}~10^{-14}\frac{\int{n_en_HdV}}{4{\pi}D_A(1+z)^2}}`

For more details on how this model was derived see 
`Hickox & Markevitch 2007 <http://adsabs.harvard.edu/abs/2007ApJ...661L.117H>`_.
This background is diffuse and uniformly fills the entire field of view of the
instrument you choose to simulate. 

Point Source Background Model
+++++++++++++++++++++++++++++

Cosmological Background Model
+++++++++++++++++++++++++++++

Creating SIMPUT Catalogs for Point Source and Cosmological Background Events
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Instrumental Background
-----------------------

Each instrument specification in the SOXS instrument registry has a default 
instrumental/particle background given by its ``"bkgnd"`` entry, which simply 
refers to the entry in the background registry:

.. code-block:: python

    from soxs import get_instrument_from_registry
    hdxi = get_instrument_from_registry("hdxi")
    print(hdxi)
 
.. code-block:: pycon

    {'arf': 'xrs_hdxi_3x10.arf',
     'bkgnd': 'acisi',
     'dither': True,
     'focal_length': 10.0,
     'fov': 20.0,
     'name': 'hdxi_3x10',
     'num_pixels': 4096,
     'psf': ['gaussian', 0.5],
     'rmf': 'xrs_hdxi.rmf'}

The ``"bkgnd"`` entry can also be set to ``None``, which corresponds to no 
particle background. To change the particle background, one would need to 
define a new instrument specification with a different background. 

Default Instrumental Backgrounds
++++++++++++++++++++++++++++++++

Lynx
~~~~

The default instrumental background in SOXS for the *Lynx* HDXI model is the 
*Chandra*/ACIS-I particle background, named ``"acisi"``, and the default 
instrumental background for the *Lynx* microcalorimeter is based on a 
model developed for the *Athena* calorimeter 
(`see here for details <http://adsabs.harvard.edu/abs/2014A%26A...569A..54L>`_), 
named ``"mucal"``.

Athena
~~~~~~

The default instrumental backgrounds in SOXS for the *Athena* WFI and 
X-IFU are based on the specifications that can be found at 
`the Athena simulation tools web portal <http://www.the-athena-x-ray-observatory.eu/resources/simulation-tools.html>`_.

Chandra
~~~~~~~

The default instrumental background in SOXS for the *Chandra* ACIS-I models is 
the *Chandra*/ACIS-I particle background, named ``"acisi"``.

Adding Your Own Instrumental Backgrounds to SOXS
++++++++++++++++++++++++++++++++++++++++++++++++

You can add your own instrumental background to the SOXS background registry. 
What you need is an ASCII table file with two columns, one with the bin energy 
in keV and the background in that bin in units of 
:math:`\rm{photons~s^{-1}~cm^{-2}~arcmin^{-2}~keV^{-1}}`. The binning must be 
linear and the bins must be equally spaced. Then you can supply it to SOXS using 
:func:`~soxs.background.instrumental.add_instrumental_background`, along with 
the name you want to give it: 

.. code-block:: python

    import soxs
    soxs.add_instrumental_background("my_particle_bkg", "my_bkg.dat")

Then you will need to create a new instrument specification (this example shows
how to clone an existing one and change the background, but one could also 
create one from scratch):

.. code-block:: python

    from soxs import get_instrument_from_registry, add_instrument_to_registry
    # Create a new instrument with that background
    new_hdxi = get_instrument_from_registry("hdxi")
    new_hdxi["name"] = "hdxi_new_bkg" # Must change the name, otherwise an error will be thrown
    new_hdxi["bkgnd"] = "my_particle_bkg"
    name = add_instrument_to_registry(new_hdxi)

Using a Background From an Event File
-------------------------------------

Creating a new background every time SOXS is run may be time-consuming for 
long exposures. SOXS provides a way to generate background events for a
particular instrument, save them to a standard event file, and then use this
file as input to :func:`~soxs.instrument.instrument_simulator`. The
:func:`~soxs.instrument.make_background_file` allows for this:

