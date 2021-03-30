.. _config:

SOXS Configuration File
=======================

In SOXS version 2.0 and later, configuration elements can be set in the 
configuration file. On most systems, this is placed in the ``XDG_CONFIG_HOME``
environment variable, which is ``$HOME/.config`` for most systems. The
SOXS configuration file is therefore ``XDG_CONFIG_HOME/soxs/soxs.cfg``.

In versions of SOXS post v3.0.0, these are the options available for 
customization in the configuration file:

.. code-block:: text

    [soxs]
    soxs_data_dir = /does/not/exist # The path to instrument files and APEC tables
    abund_table = angr # The abundance table to use for APEC thermal spectra
