#!/usr/bin/env python
from setuptools import setup, find_packages
from setuptools.extension import Extension
import numpy as np
import glob
import os


if os.name == "nt":
    std_libs = []
else:
    std_libs = ["m"]

scripts = glob.glob("scripts/*")

cython_extensions = [
    Extension("soxs.lib.broaden_lines",
              ["soxs/lib/broaden_lines.pyx"],
              language="c", libraries=std_libs,
              include_dirs=[np.get_include()]),
]

setup(name='soxs',
      packages=find_packages(),
      description='Simulated Observations of X-ray Sources',
      author='John ZuHone',
      author_email='john.zuhone@cfa.harvard.edu',
      url='https://github.com/lynx-x-ray-observatory/soxs/',
      install_requires=["numpy", "astropy>=3.0", "tqdm", "pooch", 
                        "h5py", "scipy", "pyyaml", "regions", "appdirs"],
      include_package_data=True,
      scripts=scripts,
      classifiers=[
          'Intended Audience :: Science/Research',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 3',
          'Topic :: Scientific/Engineering :: Visualization',
      ],
      ext_modules=cython_extensions,
      )
