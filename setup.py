#!/usr/bin/env python
from setuptools import setup
from setuptools.extension import Extension
import numpy as np
import glob

scripts = glob.glob("scripts/*")

cython_extensions = [
    Extension("xrs_tools.cutils",
              sources=["xrs_tools/cutils.pyx"],
              language="c", libraries=["m"],
              include_dirs=[np.get_include()])]

setup(name='xrs_tools',
      packages=['xrs_tools'],
      use_scm_version=False,
      version="0.1",
      description='Tools for X-Ray Surveyor simulations',
      author='John ZuHone',
      author_email='jzuhone@gmail.com',
      url='http://github.com/XRStools/xrs_tools',
      setup_requires=["numpy","cython>=0.24","setuptools_scm"],
      install_requires=["six","numpy","astropy"],
      scripts=scripts,
      classifiers=[
          'Intended Audience :: Science/Research',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3.5',
          'Topic :: Scientific/Engineering :: Visualization',
      ],
      ext_modules=cython_extensions,
      )
