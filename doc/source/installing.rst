.. _installing:

Installation
============

SOXS and its dependencies are installed as a standard Python package. You may use ``pip``:

.. code-block:: bash

    pip install soxs

or, if you use `Anaconda Python <https://www.continuum.io/anaconda-overview>`_, you may 
install SOXS using ``conda``:

.. code-block:: bash

    conda install -c jzuhone soxs
  
These methods install both the Python interface and the command-line scripts. 

Of course, you can always clone the source from `GitHub <http://github.com/XRStools/soxs>`_
and install it manually:

.. code-block:: bash
    
    git clone http://github.com/XRStools/soxs
    cd soxs
    python setup.py install
    
or run ``python setup.py develop`` if you want to make changes to the code and see them 
reflected without recompiling. 