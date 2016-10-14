.. _simput:

Reading and Writing SIMPUT Files
================================

Reading SIMPUT Files
--------------------

Writing SIMPUT Files
--------------------

If you have created a set of simulated events which you wish to convolve with the instrument
simulator or with some other tool, you can write them to a SIMPUT file using
:function:`~xrs_tools.simput.write_simput_phlist`. This will produce two files: the SIMPUT file
containing the parameters for the source, and a photon list file linked to the SIMPUT file which
contains the actual event energies and positions.


