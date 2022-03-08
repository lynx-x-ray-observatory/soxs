.. _config:

SOXS Configuration File
=======================

In SOXS, configuration elements can be set in the configuration file. 
On most systems, this is placed in the ``XDG_CONFIG_HOME`` environment 
variable, which is ``$HOME/.config`` for most systems. The SOXS configuration 
file is therefore ``XDG_CONFIG_HOME/soxs/soxs.cfg``.

In versions of SOXS post v3.1.0, these are the options available for 
customization in the configuration file:

.. code-block:: text

    [soxs]
    soxs_data_dir = /does/not/exist # The path to instrument files and APEC tables
    abund_table = angr # The abundance table to use for APEC thermal spectra
    apec_vers = 3.0.9 # The default version of APEC to use
    bkgnd_nH = 0.018 # neutral hydrogen column for backgrounds, units of 1e22 cm**-2
    bkgnd_absorb_model = wabs # absorption model, currently either wabs or tbabs
    frgnd_spec_model = default # foreground spectrum model, currently either default or halosat

If ``soxs_data_dir`` is not set in the configuration file, or is
set to an invalid directory, a default directory will be chosen:

.. code-block:: pycon

    soxs : [WARNING  ] 2021-04-14 22:05:49,790 Setting 'soxs_data_dir' to /Users/jzuhone/Library/Caches/soxs for this session. Please update your configuration if you want it somewhere else.

The configuration can be changed in the file, or it can be changed from within
a Python script or notebook itself, using the :func:`~soxs.utils.set_soxs_config`
function:

.. code-block::

    import soxs

    soxs.set_soxs_config("apec_vers", "2.0.1")

    # this will now use APEC version 2.0.1 by default
    agen = soxs.ApecGenerator(0.1, 10.0, 10000)

.. note::

    Changes to the ``"abund_table"``, ``"apec_vers"``, ``"bkgnd_nH"``,
    ``"bkgnd_absorb_model"``, or ``"frgnd_spec_model"`` config options using
    :func:`~soxs.utils.set_soxs_config` will trigger a re-creation of the 
    astrophysical foreground model. 

.. _mission-config:

Mission-Specific Configurations
-------------------------------

Though not normally the case, there may be cases in which it will make sense
to apply a mission-specific configuration, which will change a number of
configuration values at once. For this, SOXS provides the handy function
:func:`~soxs.utils.set_mission_config`. Currently, the only option is ``"lem"``
for the `LEM Probe Concept <https://lem.physics.wisc.edu>`_:

.. code-block::

    import soxs
    soxs.set_mission_config("lem")
