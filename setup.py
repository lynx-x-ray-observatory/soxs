#!/usr/bin/env python
from setuptools import setup
from setuptools.extension import Extension
import numpy as np
import glob

scripts = glob.glob("scripts/*")

cython_extensions = [
    Extension("sox.cutils",
              sources=["sox/cutils.pyx"],
              language="c", libraries=["m"],
              include_dirs=[np.get_include()])]

setup(name='sox',
      packages=['sox'],
      version="0.1",
      description='Simulated Observations with X-ray Surveyor',
      author='John ZuHone',
      author_email='jzuhone@gmail.com',
      url='http://github.com/XRStools/sox',
      setup_requires=["numpy","cython>=0.24"],
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
