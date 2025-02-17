

![SOXS: Simulated Observations of X-ray Sources](doc/source/images/SOXS_Wordmark.png)

[![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)
[![Build and Test](https://github.com/lynx-x-ray-observatory/soxs/actions/workflows/build-test.yml/badge.svg)](https://github.com/lynx-x-ray-observatory/soxs/actions/workflows/build-test.yml/)
[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/lynx-x-ray-observatory/soxs/main.svg)](https://results.pre-commit.ci/latest/github/lynx-x-ray-observatory/soxs/main)

# What is SOXS?

SOXS is a software suite which can create simulated X-ray observations of
astrophysical sources with almost *any* existing or planned X-ray observatory. The goal of
SOXS is to provide a comprehensive set of tools to design source models and
convolve them with simulated models of X-ray instruments. This package was originally
developed to support the [*Lynx X-ray Observatory*](https://www.lynxobservatory.org)
mission concept, but has evolved to support other missions as well.

There are two main entry points to SOXS: a command-line interface, and a
Python interface. The former is simpler to use, but the latter has more power
and flexibility. Both of these entry points are extensively documented with
examples.

# Installing SOXS

SOXS and its dependencies are installed as a standard Python package, and it is compatible
with Python 3.11 or later. You may use `pip` to install it (if you do not have pip, check
that your executable is not named `pip3`, otherwise visit https://pip.pypa.io/ to download
it):

```
[~]$ pip install soxs
```

If the Python distribution is not "owned" by you on your machine you might have to call
`sudo pip install soxs`. If you need to upgrade from a previous version of SOXS, issue
`[sudo] pip install -U soxs` from the command line.

If you use [Anaconda Python](https://www.continuum.io/anaconda-overview), you may
install SOXS using `conda`:

```
[~]$ conda install -c conda-forge soxs
```

These methods install both the Python interface and the command-line scripts.

Alternatively, to install into your Python distribution from [source](https://github.com/lynx-x-ray-observatory/soxs):

```
[~]$ pip install .
```

# Getting Help

There are a number of ways to get help with SOXS.

## Documentation

Documentation for SOXS lives at https://hea-www.cfa.harvard.edu/soxs.

## Mailing List

There's a [Google Group](https://groups.google.com/g/soxs-pyxsim-sims) to get help on SOXS and
discuss related matters.

## GitHub Issues Page

If you have a specific code issue that seems like a bug or have a feature or enhancement request,
the best place to note it is on the [GitHub issues page](https://github.com/lynx-x-ray-observatory/soxs/issues)
so that we can keep track of it.
