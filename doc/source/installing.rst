.. _installing:

Installation
============

The packages which are part of ``xrs_tools`` may be installed separately (the hard
way) or together using a script (the easy way).

The Easy Way: Installation script
---------------------------------

To download all of the packages (with some options of what to choose to install),
download the all-in-one install script and run it from the command line:

.. code-block:: bash

    bash install_script.sh

This will minimally install an Anaconda Python stack with all of the packages
included, with the exception of the following options, which may be tuned by
editing them at the top of the script:

* ``DEST_DIR``: The directory into which the packages will be installed. The default
  is ``${HOME}/xrs_tools``.
* ``PY_VERSION``: The Python version to install, whether 2.7 (set to 2) or 3.5 (set
  to 3). Default 2.
* ``INST_SIMX``: If set to 0, SIMX won't be installed, if set to 1 (default), it will.
* ``INST_PYXSIM``: If set to 0 (default), pyXSIM won't be installed, if set to 1, it will.
* ``INST_GIT``: If set to 0, git won't be installed, if set to 1 (default), it will.
  If you already have git, you should set this to 0.

The Hard Way: Custom install
----------------------------

If you'd prefer to install the various tools on your own, they may be installed
separately. For example, the Python tools may be installed into any Python
installation you choose. For Anaconda Python, use the standard ``conda install``:

.. code-block:: bash

    conda install -c xrs_tools mksrcs

or if you are using any other Python stack, use pip:

.. code-block:: bash

    pip install mksrcs

SIMX may be downloaded and installed separately from its own
`website <http://hea-www.cfa.harvard.edu/simx>`_.
