# License: BSD 3 clause
import os

import numpy


def configuration(parent_package="", top_path=None):
    from numpy.distutils.misc_util import Configuration

    libraries = []
    if os.name == "posix":
        libraries.append("m")

    config = Configuration("_hdbscan", parent_package, top_path)

    # HDBSCAN subpackage
    config.add_subpackage("tests")
    config.add_extension(
        "_linkage",
        sources=["_linkage.pyx"],
        include_dirs=[numpy.get_include()],
        libraries=libraries,
    )
    config.add_extension(
        "_reachability",
        sources=["_reachability.pyx"],
        include_dirs=[numpy.get_include()],
        libraries=libraries,
    )
    config.add_extension(
        "_tree",
        sources=["_tree.pyx"],
        include_dirs=[numpy.get_include()],
        libraries=libraries,
    )
    return config


if __name__ == "__main__":
    from numpy.distutils.core import setup

    setup(**configuration(top_path="").todict())
