#!/usr/bin/env python
from setuptools import setup

setup(name='xrs_tools',
      packages=['xrs_tools'],
      version='0.1',
      description='Tools for X-Ray Surveyor simulations',
      author='John ZuHone',
      author_email='jzuhone@gmail.com',
      url='http://github.com/XRStools/xrs_tools',
      setup_requires=["numpy"],
      install_requires=["six","numpy","astropy"],
      classifiers=[
          'Intended Audience :: Science/Research',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3.5',
          'Topic :: Scientific/Engineering :: Visualization',
      ],
      )
