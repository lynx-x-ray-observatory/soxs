.. _simput:

Working with SIMPUT Files
=========================

Writing SIMPUT Files
--------------------

If you have created a set of simulated events which you wish to convolve with the instrument
simulator or with some other tool, you can write them to a SIMPUT file using
:func:`~xrs_tools.simput.write_simput_catalog`. This will produce two files: the SIMPUT file
containing the parameters for the source, and a photon list file linked to the SIMPUT file which
contains the actual event energies and positions. For example, say we have created a 
:class:`~xrs_tools.spectra.Spectrum` object named ``spec`` which we've generated energies 
from, and that we're going to assign these energies to a point source. We can then create 
the SIMPUT file and photon list file like this:

.. code-block:: python

    from xrs_tools import write_simput_catalog
    
    exp_time = 100000. # in seconds
    area = 30000. # in cm^2
    energies = spec.generate_energies(exp_time, area)
    num_events = len(energies)
    ra = 30.0*np.ones(num_events) # RA in degrees
    dec = 45.0*np.ones(num_events) # Dec in degrees
    
    write_simput_catalog("point_source", "source1", 
                         exp_time, area, ra, dec, energies, 
                         clobber=True) 
                         
We have to give :func:`~xrs_tools.simput.write_simput_catalog` the exposure time and area because
we need to compute a flux for the source.

Alternatively, you may already have a SIMPUT file associated with a photon list file, but want to 
add another source to the same SIMPUT catalog. You can accomplish this by making the same call to
:func:`~xrs_tools.simput.write_simput_catalog` but setting ``append=True``:

.. code-block:: python

    write_simput_catalog("point_source", "source2", 
                         exp_time, area, ra, dec, energies, 
                         append=True) 

Reading SIMPUT Files
--------------------
