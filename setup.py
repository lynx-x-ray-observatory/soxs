#!/usr/bin/env python
from setuptools import setup
from setuptools.extension import Extension
import numpy as np
import glob

scripts = glob.glob("scripts/*")

cython_extensions = [
    Extension("soxs.cutils",
              sources=["soxs/cutils.pyx"],
              language="c", libraries=["m"],
              include_dirs=[np.get_include()])]

setup(name='soxs',
      packages=['soxs'],
      version="0.5.1",
      description='Simulated Observations with X-ray Surveyor',
      author='John ZuHone',
      author_email='john.zuhone@cfa.harvard.edu',
      url='http://github.com/XRStools/soxs',
      setup_requires=["numpy","cython>=0.24"],
      install_requires=["six","numpy","astropy","tqdm"],
      include_package_data=True,
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
