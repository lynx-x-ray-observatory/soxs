.. _installing:

Installation
============

SOXS and its dependencies are installed as a standard Python package, and it is 
compatible with Python 3.8 and higher. You may use ``pip`` to install it (if 
you do not have pip, check that your executable is not named ``pip3``, otherwise 
visit https://pip.pypa.io/ to download it):

.. code-block:: bash

    pip install soxs

If the Python distribution is not "owned" by you on your machine you might have
to call ``sudo pip install soxs``. If you need to upgrade from a previous 
version of SOXS, issue ``[sudo] pip install -U soxs`` from the command line. 

If you use `Anaconda Python <https://www.continuum.io/anaconda-overview>`_, you
may install SOXS using ``conda``:

.. code-block:: bash

    conda install -c jzuhone -c astropy soxs
  
Note both the ``jzuhone`` and ``astropy`` channels are required. These methods 
install both the Python interface and the command-line scripts. 

.. warning::

    Currently, there is no Anaconda package for 
    `regions <https://astropy-regions.readthedocs.io/>`_ on Python 3.10, which 
    is a SOXS dependency. It must be installed via pip. 

Of course, you can always clone the source from 
`GitHub <https://github.com/lynx-x-ray-observatory/soxs>`_ and install it 
manually:

.. code-block:: bash
    
    git clone https://github.com/lynx-x-ray-observatory/soxs
    cd soxs
    python setup.py install
    
or run ``python setup.py develop`` instead if you want to make changes to the 
code and see them reflected without recompiling (though if you make updates to 
the command-line scripts you will have to run ``python setup.py develop`` 
again). 

SOXS Dependencies
=================

SOXS has the following Python dependencies:

* `NumPy <https://numpy.org>`_
* `AstroPy <https://www.astropy.org>`_
* `SciPy <https://www.scipy.org>`_
* `h5py <https://www.h5py.org>`_
* `tqdm <https://github.com/noamraph/tqdm>`_
* `regions <https://astropy-regions.readthedocs.io/>`_
* `pooch <https://www.fatiando.org/pooch>`_

Using any installation method, these dependencies should automatically install 
(if you do not already have them) provided you are connected to the internet.

Optional Packages
=================

There are also a number of optional packages which may be used with SOXS to
enhance its capabilities. 

* To make mock X-ray observations from 3D hydrodynamics models and other models
  of astrophysical sources, use the 
  `pyXSIM <http://hea-www.cfa.harvard.edu/~jzuhone/pyxsim>`_ package with SOXS.
* To make mosaics of SOXS observations (see :ref:`mosaic`), install the 
  `reproject <https://reproject.readthedocs.io>`_ package.