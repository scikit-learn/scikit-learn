import os

import numpy
from numpy.distutils.misc_util import Configuration


def configuration(parent_package="", top_path=None):
    config = Configuration("manifold", parent_package, top_path)
    config.add_extension("_binary_search",
                         sources=["_binary_search.c"],
                         include_dirs=[numpy.get_include()],
                         extra_compile_args=["-O3"])

    config.add_subpackage("tests")

    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(**configuration().todict())
