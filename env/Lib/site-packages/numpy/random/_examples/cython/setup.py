#!/usr/bin/env python3
"""
Build the Cython demonstrations of low-level access to NumPy random

Usage: python setup.py build_ext -i
"""

import numpy as np
from distutils.core import setup
from Cython.Build import cythonize
from setuptools.extension import Extension
from os.path import join, dirname

path = dirname(__file__)
defs = [('NPY_NO_DEPRECATED_API', 0)]

extending = Extension("extending",
                      sources=[join(path, 'extending.pyx')],
                      include_dirs=[
                            np.get_include(),
                            join(path, '..', '..')
                        ],
                      define_macros=defs,
                      )
distributions = Extension("extending_distributions",
                          sources=[join(path, 'extending_distributions.pyx')],
                          include_dirs=[np.get_include()],
                          define_macros=defs,
                         )

extensions = [extending, distributions]

setup(
    ext_modules=cythonize(extensions)
)
