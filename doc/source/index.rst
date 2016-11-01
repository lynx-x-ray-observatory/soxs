.. sox documentation master file, created by
   sphinx-quickstart on Sat Aug 13 15:23:51 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

SOXS: Simulated Observations with X-ray Surveyor
================================================

.. raw:: html

   <figure style="display: table; float: right; margin: 0 0 20px 20px;">
   <img src="galaxy_compare.png" width="600" style="float: right;"/>
   <figcaption style="display: table-caption; caption-side: bottom;">
   Simulated X-ray events from a spiral galaxy from the Illustris simulation,
   assuming <em>Chandra</em>/ACIS-I on the left and X-ray Surveyor/HDXI on the 
   right. Both exposures are 50 ks. 
   </figcaption>
   </figure>

SOXS is a software suite which creates simulated X-ray observations of
astrophysical sources with the mission concept X-ray Surveyor. The goal of 
SOXS is to provide a comprehensive set of tools to design source models and
convolve them with simulated models of the X-ray Surveyor instrument. 

There are two main entry points to SOXS: a command-line interface, and a
Python interface. The former is simpler to use, but the latter has more power
and flexibility. Both of these entry points are extensively documented here with
examples. Though you will find details on usage and example runs of the command
line scripts in the :ref:`command-line`, it is recommended to also look over
the corresponding documentation in the :ref:`python` for details about what is 
going on under the hood. 

Why Another Mock X-ray Observation Package?
-------------------------------------------

There are already a number of successful efforts to create mock X-ray 
observations, and it is a sensible question as to why it would not have made
more sense to adapt one of the existing projects for X-ray Surveyor simulations.
There are three basic reasons why we have chosen to develop a new package:

Firstly, SOXS is explicitly geared towards simulating observations of X-ray
sources with X-ray Surveyor. This tool is designed to provide support for developing
a science case for the X-ray Surveyor mission and help drive the design of the 
instruments. Having a standalone package (that still plays nice with others) that
serves as a "one-stop shop" for X-ray Surveyor simulations simplfies this task. 

Secondly, SOXS is being developed in Python, reflecting the growing popularity 
of the use of Python in astronomy and astrophysics. Though there are a number
of command-line scripts provided in SOXS, which can carry out the most important
tasks, the Python interface to SOXS is more powerful and flexible. Working with the
Python interface also provides a direct connection to other Python packages for 
science in general and astronomy in particular, including `NumPy <http://www.numpy.org>`_, 
`SciPy <http://www.scipy.org>`_, `AstroPy <http://www.astropy.org>`_, `yt <http://yt-project.org/>`_,
and `pyXSIM <http://hea-www.cfa.harvard.edu/~jzuhone/pyXSIM>`_. 

Thirdly, SOXS is being developed out in the open, on `GitHub <http://github.com/XRStools/soxs>`_, 
to encourage contributions in the form of bugfixes and enhancements, and to make
contributing as simple as forking the code and submitting a pull request for review.

License
-------

SOXS is released under a `BSD 3-clause license <https://opensource.org/licenses/BSD-3-Clause>`_.

Current Version
---------------

The current stable version is 0.2.0. See the :ref:`changelog` for details on changes from previous
versions.

Documentation Contents
----------------------

.. toctree::
   :maxdepth: 2

   installing
   overview
   responses
   command_line/index
   python/index
   cookbook/index
   getting_help
   api/index


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

