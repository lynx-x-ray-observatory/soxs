#!/usr/bin/env python
from setuptools import setup
from setuptools.extension import Extension
import numpy as np

VERSION = '0.1.dev0'

cython_extensions = [
    Extension("pyxsim.cutils",
              sources=["pyxsim/cutils.pyx"],
              language="c", libraries=["m"],
              include_dirs=[np.get_include()])]

setup(name='xrs_tools',
      packages=['xrs_tools'],
      version=VERSION,
      description='Tools for X-Ray Surveyor simulations',
      author='John ZuHone',
      author_email='jzuhone@gmail.com',
      url='http://github.com/XRStools/xrs_tools',
      setup_requires=["numpy","cython>=0.24"],
      install_requires=["six","numpy","astropy"],
      classifiers=[
          'Intended Audience :: Science/Research',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3.5',
          'Topic :: Scientific/Engineering :: Visualization',
      ],
      ext_modules=cython_extensions,
      )
