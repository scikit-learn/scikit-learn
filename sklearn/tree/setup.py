import os

import numpy
from numpy.distutils.misc_util import Configuration

from sklearn._build_utils import add_cython_extension


def configuration(parent_package="", top_path=None):
    config = Configuration("tree", parent_package, top_path)
    libraries = []
    if os.name == 'posix':
        libraries.append('m')
    add_cython_extension(top_path,
                         config,
                         "_tree",
                         sources=["_tree.c"],
                         include_dirs=[numpy.get_include()],
                         libraries=libraries,
                         extra_compile_args=["-O3"])
    add_cython_extension(top_path,
                         config,
                         "_splitter",
                         sources=["_splitter.c"],
                         include_dirs=[numpy.get_include()],
                         libraries=libraries,
                         extra_compile_args=["-O3"])
    add_cython_extension(top_path,
                         config,
                         "_criterion",
                         sources=["_criterion.c"],
                         include_dirs=[numpy.get_include()],
                         libraries=libraries,
                         extra_compile_args=["-O3"])
    add_cython_extension(top_path,
                         config,
                         "_utils",
                         sources=["_utils.c"],
                         include_dirs=[numpy.get_include()],
                         libraries=libraries,
                         extra_compile_args=["-O3"])

    config.add_subpackage("tests")

    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(**configuration().todict())
