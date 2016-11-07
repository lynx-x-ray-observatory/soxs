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

which prints:

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

Applying instrumental and astrophysical backgrounds are handled somewhat differently. Each instrument 
specification has a default instrumental/particle background given in its entry in the SOXS instrument 
registry, which simply refers to the entry in the background registry. To change the instrumental background,
one would need to define a new instrument specification with a different background. The default instrumental
background in SOXS for the HDXI is the *Chandra*/ACIS-I particle background, named ``"acisi"``, and the default
instrumental background for the calorimeter is from a model developed for the *Athena* calorimeter 
(`see here for details <http://adsabs.harvard.edu/abs/2014A%26A...569A..54L>`_), named ``"mucal"``.

Astrophysical backgrounds are not tied to a particular instrument specification and can be specified in the
call to :func:`~soxs.instrument.instrument_simulator`. To specify a particular astrophysical background,
simply supply the name of the background in the background registry, or supply ``None`` to turn it off
entirely:

.. code-block:: python

    # uses the default astrophysical background
    instrument_simulator(simput_file, out_file, exp_time, instrument, 
                         sky_center, clobber=True, astro_bkgnd="hm_cxb")
                          
.. code-block:: python

    # turns off the astrophysical background entirely
    instrument_simulator(simput_file, out_file, exp_time, instrument, 
                         sky_center, clobber=True, astro_bkgnd=None)

The default astrophysical background in SOXS is from 
`Hickox & Markevitch 2007 <http://adsabs.harvard.edu/abs/2007ApJ...661L.117H>`_, named ``"hm_cxb"``.

Adding Your Own Backgrounds to SOXS
-----------------------------------

You can add your own background to the SOXS background registry. What you need is an ASCII table
file with two columns, one with the bin energy in keV and the background in that bin in units of 
:math:`\rm{photons~s^{-1}~cm^{-2}~arcmin^{-2}~keV^{-1}}`. The binning must be linear and the bins 
must be equally spaced. Then you can supply it to SOXS using 
:func:`~soxs.background.add_background_to_registry`, along with the name you want to give it and
the background type, ``"instrumental"``, or ``"astrophysical"``:

.. code-block:: python

    import soxs
    soxs.add_background_to_registry("zuhone_bkg", "my_bkg.dat", "astrophysical")

Then, in the case of an astrophysical background, it can be used with the instrument simulator 
in the same way:

.. code-block:: python

    # uses your astrophysical background
    instrument_simulator(simput_file, out_file, exp_time, instrument, 
                         sky_center, clobber=True, astro_bkgnd="zuhone_bkg")

In the case of an instrmental background, you will need to create a new instrument specification:

.. code-block:: python

    from soxs import get_instrument_from_registry, add_instrument_to_registry, \
        add_background_to_registry
    # First add the background
    add_background_to_registry("my_particle_bkg", "my_pbkg.dat", "instrumental")
    # Then create a new instrument with that background
    new_hdxi = get_instrument_from_registry("hdxi")
    new_hdxi["name"] = "hdxi_new_bkg" # Must change the name, otherwise an error will be thrown
    new_hdxi["bkgnd"] = "my_particle_bkg"
    name = add_instrument_to_registry(new_hdxi)

