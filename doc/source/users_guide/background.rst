.. _background:

Background Models in SOXS
=========================

SOXS simulates background for every observation. The background in SOXS is
comprised of three components: a uniform galactic foreground, a point-source 
background, and an instrumental/particle background. The former two components
are not tied to any particular instrument specification, whereas the latter 
depends on the instrument being simulated. We will describe each of these
background components in turn. 

Galactic Foreground Model
-------------------------

The galactic foreground component is modeled as a sum of two thermal models, 
``apec+apec``, with parameters:

``apec`` **Model 1**

* ``kT``: :math:`\rm{0.2~keV}`
* ``abund``: :math:`\rm{1.0~Z_\odot}`
* ``redshift``: :math:`0.0`
* ``norm``: :math:`\rm{6.82251 \times 10^{-7}~10^{-14}\frac{\int{n_en_HdV}}{4{\pi}D_A^2(1+z)^2}}`
 
``apec`` **Model 2**

* ``kT``: :math:`\rm{0.099~keV}`
* ``abund``: :math:`\rm{1.0~Z_\odot}`
* ``redshift``: :math:`0.0`
* ``norm``: :math:`\rm{1.12328 \times 10^{-6}~10^{-14}\frac{\int{n_en_HdV}}{4{\pi}D_A^2(1+z)^2}}`

This fit is based on the results of 
`Hickox & Markevitch 2007 <http://adsabs.harvard.edu/abs/2007ApJ...661L.117H>`_.
This background is diffuse and uniformly fills the entire field of view of the
instrument you choose to simulate. 

.. _ptsrc-bkgnd:

Point Source Background Model
-----------------------------

Another astrophysical component of the background in SOXS comes from resolved
point sources. The emission of these sources is assumed to originate from 
cosmologically distant AGN and galaxies. The fluxes for these sources are drawn
from :math:`\rm{log}~N-\rm{log}~S` distributions taken from
`Lehmer et al. (2012) <http://adsabs.harvard.edu/abs/2012ApJ...752...46L>`_'s
study of the *Chandra* Deep Field South. The point sources have fluxes in the 
0.5-2 keV band in the :math:`7.63 \times 10^{-22} - 1.0 \times 10^{-13}~\rm{erg}~\rm{s}^{-1}~\rm{cm}^{-2}` 
range.

Each point source is given a power-law spectrum. Galaxies are assumed to have a
spectral index of :math:`\alpha = 2.0`. The spectral indices of AGN sources are
drawn from a fit to the spectral index distribution of sources given in 
Figure 13a of `Hickox & Markevitch 2006 <http://adsabs.harvard.edu/abs/2006ApJ...645...95H>`_. 
Sources are absorbed by foreground Galactic neutral hydrogen assuming a neutral 
hydrogen column of :math:`n_H = 5 \times 10^{20}~\rm{cm}^{-2}`. The position of 
each point source is uniformly randomly distributed within the field of view. 

Though a point-source population is automatically created as a background 
component when an observation is simulated, one can also create a SIMPUT catalog
of point sources using the same machinery, with the ability more finely control
the input parameters of the model. For more information, see :ref:`point-source-catalog`.

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

Turning Background Components On and Off
----------------------------------------

All components of the background are turned on in the instrument simulator by
default. The various components of the background can be turned on or off 
entirely in the call to :func:`~soxs.instrument.instrument_simulator` by setting
the parameters ``ptsrc_bkgnd``, ``foreground``, and/or ``instr_bkgnd`` to 
``True`` or ``False``:

.. code-block:: python

    # turn off the astrophysical foreground
    soxs.instrument_simulator(simput_file, out_file, exp_time, instrument, 
                              sky_center, overwrite=True, foreground=False)

.. code-block:: python

    # turn off the instrumental background
    soxs.instrument_simulator(simput_file, out_file, exp_time, instrument, 
                              sky_center, overwrite=True, instr_bkgnd=False)

.. code-block:: python

    # turn off all backgrounds
    soxs.instrument_simulator(simput_file, out_file, exp_time, instrument, 
                              sky_center, overwrite=True, ptsrc_bkgnd=False,
                              instr_bkgnd=False, foreground=False)

.. _make-bkgnd:

Using a Background From an Event File
-------------------------------------

Creating a new background every time SOXS is run may be time-consuming for 
long exposures. SOXS provides a way to generate background events for a
particular instrument, save them to a standard event file, and then use this
file as input to :func:`~soxs.instrument.instrument_simulator`. The
:func:`~soxs.instrument.make_background_file` allows for this:

.. code-block:: python

    out_file = 'bkgnd_evt.fits'
    exp_time = 1000000.0 # seconds
    instrument = "hdxi"
    sky_center = [24., 12.] # degrees
    soxs.make_background_file(out_file, exp_time, instrument, sky_center, 
                              overwrite=True, foreground=True, instr_bkgnd=False,
                              ptsrc_bkgnd=True)

As can be noted from this example, :func:`~soxs.instrument.make_background_file`
allows one to turn any of the three background components on or off using the
boolean arguments ``foreground``, ``instr_bkgnd``, or ``ptsrc_bkgnd``. 

One can also take finer control of the point-source contribution to the background
by supplying an ASCII table of point-source properties generated by 
:func:`~soxs.background.make_point_sources_file` or :func:`~soxs.background.make_point_source_list`
using the ``input_sources`` keyword argument:

.. code-block:: python

    fov = 20.0 # arcmin
    out_file = 'bkgnd_evt.fits'
    exp_time = 1000000.0 # seconds
    instrument = "hdxi"
    sky_center = [24., 12.] # degrees
    soxs.make_point_source_list('my_ptsrc.dat', exp_time, fov, sky_center)
    soxs.make_background_file(out_file, exp_time, instrument, sky_center, 
                              overwrite=True, input_sources='my_ptsrc.dat')

See :ref:`point-source-list` for more information on this feature. 

:func:`~soxs.instrument.instrument_simulator` can use this background file when
it is supplied with the ``bkgnd_file`` argument, provided that the same
instrument is used and the exposure time of the source observation is not longer
than the exposure time that the background was run with:

.. code-block:: python

    exp_time = 500000.0 # seconds
    instrument = "hdxi"
    simput_file = "beta_model_simput.fits"
    out_file = "evt.fits"
    sky_center = [30., 45.]
    soxs.instrument_simulator(simput_file, out_file, exp_time, instrument, 
                              sky_center, overwrite=True, bkgnd_file="bkgnd_evt.fits") 

Note that the pointing of the background event file does not to be the same as
the source pointing--the background events will be re-projected to match the
pointing of the source. 