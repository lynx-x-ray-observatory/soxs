.. _charge-exchange:

Charge Exchange Spectra
=======================

The phenomenon of charge exchange is a process in which a neutral atom or molecule
collides with an ion, resulting in the transfer of an electron from the neutral atom
to the ion. The recombined ion is then left in an excited state, which can lead to
the emission of X-rays as the ion returns to its ground state. Charge exchange lines
may be expected in astrophysical situations where neutral atoms or molecules are in the
presence of highly ionized plasma, and there is a significant relative velocity between
the two.

Charge exchange spectra are typically characterized by a series of emission lines
(and continuum 2-photon emission from 2s-1s transitions). To implement charge
exchange spectra in SOXS, we use the ``acx2`` package, which is an optional dependency
and is not installed automatically when you install SOXS.

Installation
------------

.. |pyatomdb| replace:: install the ``pyatomdb`` package
.. _pyatomdb: https://atomdb.readthedocs.io/en/master/installation.html

.. |acx2| replace:: install the ``acx2`` package
.. _acx2: https://acx2.readthedocs.io/en/latest/#installation

To use the charge exchange model in SOXS, you must first |pyatomdb|_, which is a
Python package for atomic data and modeling of X-ray spectra. Following that, you
can |acx2|_.

Usage
-----

There are two ways to use the charge exchange model in SOXS.
