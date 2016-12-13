.. _background:

Simulating Background in SOXS
=============================

SOXS contains support for including models of astrophysical and instrumental/particle backgrounds (or foregrounds)
in a simulation of an observation. Backgrounds are turned on in the instrument simulator by default. 

Backgrounds in SOXS, both astrophysical and particle, are stored in a background registry, which is a 
Python dictionary. Each background is a ``BackgroundSpectrum`` object, which simply stores the background
spectrum itself (in units of :math:`\rm{photons~s^{-1}~cm^{-2}~arcmin^{-2}~keV^{-1}}`) as well as the background
type. You can see what backgrounds are registered by calling :func:`~soxs.background.show_background_registry`:

.. code-block:: python
    
    from soxs import show_background_registry
    show_background_registry()

which prints something like:

.. code-block:: pycon

    Background: hm_cxb
        Type: astrophysical
        Total Flux (0.1 keV - 10.0 keV): 1.100008388272249e-05 ph / (arcmin2 cm2 s)
    Background: mucal
        Type: instrumental
        Total Flux (0.1 keV - 10.0 keV): 0.0029680200000000035 ph / (arcmin2 cm2 s)
    Background: acisi
        Type: instrumental
        Total Flux (0.1 keV - 10.0 keV): 0.004312855625210535 ph / (arcmin2 cm2 s)
    ...
    
Applying instrumental and astrophysical backgrounds are handled somewhat differently. Each instrument 
specification in the SOXS instrument registry has a default instrumental/particle background given by its ``"bkgnd"``
entry, which simply refers to the entry in the background registry:

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

The ``"bkgnd"`` entry can also be set to ``None``, which corresponds to no particle background. To change 
the particle background, one would need to define a new instrument specification with a different background. 

The default instrumental background in SOXS for the *Lynx* HDXI is the *Chandra*/ACIS-I particle 
background, named ``"acisi"``, and the default instrumental background for the *Lynx* microcalorimeter 
is based on a model developed for the *Athena* calorimeter 
(`see here for details <http://adsabs.harvard.edu/abs/2014A%26A...569A..54L>`_), named ``"mucal"``.

The default instrumental backgrounds in SOXS for the *Athena* WFI and X-IFU are based on the specifications
that can be found at `the Athena simulation tools web portal <http://www.the-athena-x-ray-observatory.eu/resources/simulation-tools.html>`_.

The astrophysical background is not tied to a particular instrument specification. It can be turned off
entirely in the call to :func:`~soxs.instrument.instrument_simulator` by setting ``astro_bkgnd=False``:

.. code-block:: python

    # turn off the astrophysical background
    instrument_simulator(simput_file, out_file, exp_time, instrument, 
                         sky_center, clobber=True, astro_bkgnd=False)

The default astrophysical background in SOXS is from 
`Hickox & Markevitch 2007 <http://adsabs.harvard.edu/abs/2007ApJ...661L.117H>`_, named ``"hm_cxb"``, and
is modeled as a sum of two thermal models, ``apec+apec``, with parameters:

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

Adding Your Own Backgrounds to SOXS
-----------------------------------

You can add your own instrumental background to the SOXS background registry. What you need is an 
ASCII table file with two columns, one with the bin energy in keV and the background in that bin in 
units of :math:`\rm{photons~s^{-1}~cm^{-2}~arcmin^{-2}~keV^{-1}}`. The binning must be linear and 
the bins must be equally spaced. Then you can supply it to SOXS using 
:func:`~soxs.background.add_background_to_registry`, along with the name you want to give it and
the background type, which in this case is ``"instrumental"``:

.. code-block:: python

    import soxs
    soxs.add_background_to_registry("my_particle_bkg", "my_bkg.dat", "instrumental")

Then you will need to create a new instrument specification:

.. code-block:: python

    from soxs import get_instrument_from_registry, add_instrument_to_registry
    # Create a new instrument with that background
    new_hdxi = get_instrument_from_registry("hdxi")
    new_hdxi["name"] = "hdxi_new_bkg" # Must change the name, otherwise an error will be thrown
    new_hdxi["bkgnd"] = "my_particle_bkg"
    name = add_instrument_to_registry(new_hdxi)

