.. _simput:

Working with SIMPUT Files
=========================

The default storage format for unconvolved events in SOXS is SIMPUT, which is 
fast becoming a standard for making mock X-ray observations. There are two 
functions to read and write SIMPUT files in SOXS.

Writing SIMPUT Files
--------------------

If you have created a set of simulated events which you wish to convolve with 
the instrument simulator or with some other tool, you can write them to a SIMPUT 
file using :func:`~soxs.simput.write_photon_list`. This will produce two files: 
the SIMPUT filecontaining the parameters for the source, and a photon list file 
linked to the SIMPUT file which contains the actual event energies and 
positions. For example, say we have created a :class:`~soxs.spectra.Spectrum` 
object named ``spec`` which we've generated energies from, and that we're going 
to assign these energies to a point source. We can then create the SIMPUT file 
and photon list file like this:

.. code-block:: python

    from soxs import write_photon_list, PointSource
    
    exp_time = (500.0, "ks")
    area = (3.0, "m**2")
    energies = spec.generate_energies(exp_time, area)
    num_events = len(energies)
    pt_src = PointSource(30.0, 45.0)
    
    write_photon_list("point_source", "source1", energies.flux,
                      pt_src.ra, pt_src.dec, energies, overwrite=True)
                         
The ``energies`` returned by :meth:`~soxs.spectra.Spectrum.generate_energies` 
is an augmented NumPy array with unit information and the value of the flux 
for that set of energies. This flux needs to be passed to 
:func:`~soxs.simput.write_photon_list` as the third argument.

Alternatively, you may already have a SIMPUT file associated with a photon 
list file, but want to add another source to the same SIMPUT catalog. You can
accomplish this by making the same call to 
:func:`~soxs.simput.write_photon_list` but setting ``append=True``:

.. code-block:: python

    write_photon_list("point_source", "source2" energies2.flux,
                      pt2.ra, pt2.dec, energies2, append=True) 

SOXS will give each photon list source a name in the catalog, determined by the
scheme ``"soxs_src_n"`` where ``n`` is the n-th source in the file, but you can 
supply an alternative name for the source in the call to 
:func:`~soxs.simput.write_photon_list` using the ``src_name`` keyword argument: 

.. code-block:: python

    write_photon_list("point_source", "source2" energies2.flux,
                      pt2.ra, pt2.dec, energies2, append=True, 
                      src_name="my_point_source") 

Reading SIMPUT Files
--------------------

A SIMPUT catalog can be read using :func:`~soxs.simput.read_simput_catalog`:

.. code-block:: python

    from soxs import read_simput_catalog
    events, parameters = read_simput_catalog("point_source_simput.fits")
    
It returns two arguments, ``events`` and ``parameters``. ``events`` is a list of 
Python dictionaries, one for each source in the file. Each dictionary contains 
NumPy arrays for the positions and energies of the events. For example, for a 
catalog which only has one source they would look like this:

.. code-block:: python

    print(events)
    
.. code-block:: pycon

    [{'dec': array([ 44.98377818,  44.99404092,  44.99444754, ...,  45.00548515,
             45.0052105 ,  45.00658426]),
      'energy': array([ 5.11127663,  0.58575863,  2.00386882, ...,  1.09081411,
             1.31414783,  2.21034932], dtype=float32),
      'ra': array([ 30.2032835 ,  29.95447951,  29.95380409, ...,  30.04756871,
             30.04568841,  30.04643141])}]

.. code-block:: python

    print(parameters)
    
.. code-block:: pycon

    {'emax': array([ 10.99995703]), 
     'flux': array([  1.12239243e-11]), 
     'emin': array([ 0.12598762])}

Energies are in keV, flux is in :math:`{\rm erg~s^{-1}~cm^{-2}}`, and sky 
coordinates are in degrees. :func:`~soxs.simput.read_simput_catalog` is used by
the instrument simulator to read sources from a SIMPUT catalog. 