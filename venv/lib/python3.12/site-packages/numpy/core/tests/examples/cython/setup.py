"""
Provide python-space access to the functions exposed in numpy/__init__.pxd
for testing.
"""

import os
from distutils.core import setup

import numpy as np
from Cython.Build import cythonize
from setuptools.extension import Extension

macros = [("NPY_NO_DEPRECATED_API", 0)]

checks = Extension(
    "checks",
    sources=[os.path.join(".", "checks.pyx")],
    include_dirs=[np.get_include()],
    define_macros=macros,
)

extensions = [checks]

setup(ext_modules=cythonize(extensions))
