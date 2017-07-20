#!/usr/bin/env python
from setuptools import setup, find_packages
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
      packages=find_packages(),
      version="1.3.0",
      description='Simulated Observations with X-ray Surveyor',
      author='John ZuHone',
      author_email='john.zuhone@cfa.harvard.edu',
      url='https://github.com/XRStools/soxs/',
      setup_requires=["numpy","cython>=0.24"],
      install_requires=["six","numpy","astropy>=1.3","tqdm",
                        "h5py","scipy","pyyaml","pyregion"],
      include_package_data=True,
      scripts=scripts,
      classifiers=[
          'Intended Audience :: Science/Research',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3.5',
          'Programming Language :: Python :: 3.6',
          'Topic :: Scientific/Engineering :: Visualization',
      ],
      ext_modules=cython_extensions,
      )
