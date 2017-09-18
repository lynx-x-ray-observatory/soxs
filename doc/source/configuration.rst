.. _config:

SOXS Configuration File
=======================

In SOXS version 2.0 and later, configuration elements can be set in the 
configuration file. On most systems, this is placed in the ``XDG_CONFIG_HOME``
environment variable, which is ``$HOME/.config`` for most systems. The
SOXS configuration file is therefore ``XDG_CONFIG_HOME/soxs/soxs.cfg``.

In version 2.0, these are the options available for customization in the
configuration file:

.. code-block:: text

    [soxs]
    response_path = /does/not/exist # The path to the ARF and RMF files
    abund_table = angr # The abundance table to use for APEC thermal spectra
