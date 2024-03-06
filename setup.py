#!/usr/bin/env python
import glob
import os

import numpy as np
from setuptools import find_packages, setup
from setuptools.extension import Extension

if os.name == "nt":
    std_libs = []
else:
    std_libs = ["m"]

scripts = glob.glob("scripts/*")

cython_extensions = [
    Extension(
        "soxs.lib.broaden_lines",
        ["soxs/lib/broaden_lines.pyx"],
        language="c",
        libraries=std_libs,
        include_dirs=[np.get_include()],
    ),
    Extension(
        "soxs.lib.psf_cdf",
        ["soxs/lib/psf_cdf.pyx"],
        language="c",
        libraries=std_libs,
        include_dirs=[np.get_include()],
    ),
]

setup(
    packages=find_packages(),
    url="https://github.com/lynx-x-ray-observatory/soxs/",
    include_package_data=True,
    scripts=scripts,
    ext_modules=cython_extensions,
)
