.. soxs documentation master file, created by
   sphinx-quickstart on Sat Aug 13 15:23:51 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

SOXS: Simulated Observations of X-ray Sources
=============================================

.. image:: images/SOXS_Wordmark.png
   :width: 70%
   :align: center

SOXS is a software suite which creates simulated X-ray observations of
astrophysical sources. The goal of SOXS is to provide a comprehensive set 
of tools to design source models and convolve them with simulated models 
of X-ray observatories. In particular, SOXS is the primary simulation tool
for simulations of *Lynx* observations.  

There are two main entry points to SOXS: a Python interface, and a command-line
interface. The former has more power and flexibility, but the latter is often
simpler to use. Both of these entry points are extensively documented here with
examples. Though you will find details on usage and example runs of the command
line scripts in the :ref:`command-line`, it is recommended to also look over
the corresponding documentation in the :ref:`python` for details about what is 
going on under the hood. 

Why Another Mock X-ray Observation Package?
-------------------------------------------

.. raw:: html

   <figure style="display: table; float: right; margin: 0 0 20px 20px;">
   <img src="_images/soxs_showcase.png" width="600" style="float: right;"/>
   </figure>

There are already a number of successful efforts to create mock X-ray 
observations, and it is a sensible question as to why it would not have made
more sense to adapt one of the existing projects for *Lynx* simulations. There
are three basic reasons why we have chosen to develop a new package:

First, SOXS is explicitly geared towards simulating observations of X-ray
sources with *Lynx* (formerly "X-ray Surveyor"). This tool is designed to 
provide support for developing a science case for the *Lynx* mission concept 
and help drive the design of the instruments. Having a standalone package 
(that still plays nice with others) that serves as a "one-stop shop" for *Lynx*
simulations simplfies this task. 

Second, SOXS is being developed in Python, reflecting the growing popularity 
of the use of Python in astronomy and astrophysics. Though there are a number
of command-line scripts provided in SOXS, which can carry out the most important
tasks, the Python interface to SOXS is more powerful and flexible. Working with 
the Python interface also provides a direct connection to other Python packages
for science in general and astronomy in particular, including 
`NumPy <http://www.numpy.org>`_, `SciPy <http://www.scipy.org>`_, 
`AstroPy <http://www.astropy.org>`_, `yt <http://yt-project.org/>`_, and 
`pyXSIM <http://hea-www.cfa.harvard.edu/~jzuhone/pyXSIM>`_. 

Thirdly, SOXS is being developed out in the open, on 
`GitHub <http://github.com/lynx-x-ray-observatory/soxs>`_, 
to encourage contributions in the form of bugfixes and enhancements, and to make
contributing as simple as forking the code and submitting a pull request for 
review.

License
-------

SOXS is released under a `BSD 3-clause license <https://opensource.org/licenses/BSD-3-Clause>`_.

Current Version
---------------

The current stable version is 3.4.0. See the :ref:`changelog` for details on changes 
from previous versions.

Documentation Contents
----------------------

.. toctree::
   :maxdepth: 2

   installing
   overview
   configuration
   responses
   users_guide/index
   command_line/index
   cookbook/index
   getting_help
   api/index

Contributors
------------

* John ZuHone (CfA)
* Scott Randall (CfA)
* Felipe Andrade-Santos (CfA)
* Herve Bourdin (CfA)
* Alexey Vikhlinin (CfA)
* Reinout Van Weeren (CfA)
* Akos Bogdan (CfA)
* Grant Tremblay (CfA)

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

