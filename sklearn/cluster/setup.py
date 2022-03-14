# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
# License: BSD 3 clause
import os

import numpy


def configuration(parent_package="", top_path=None):
    from numpy.distutils.misc_util import Configuration

    libraries = []
    if os.name == "posix":
        libraries.append("m")

    config = Configuration("cluster", parent_package, top_path)

    config.add_extension(
        "_dbscan_inner",
        sources=["_dbscan_inner.pyx"],
        include_dirs=[numpy.get_include()],
        language="c++",
    )

    config.add_extension(
        "_hierarchical_fast",
        sources=["_hierarchical_fast.pyx"],
        language="c++",
        include_dirs=[numpy.get_include()],
        libraries=libraries,
    )

    config.add_extension(
        "_k_means_common",
        sources=["_k_means_common.pyx"],
        include_dirs=[numpy.get_include()],
        libraries=libraries,
    )

    config.add_extension(
        "_k_means_lloyd",
        sources=["_k_means_lloyd.pyx"],
        include_dirs=[numpy.get_include()],
        libraries=libraries,
    )

    config.add_extension(
        "_k_means_elkan",
        sources=["_k_means_elkan.pyx"],
        include_dirs=[numpy.get_include()],
        libraries=libraries,
    )

    config.add_extension(
        "_k_means_minibatch",
        sources=["_k_means_minibatch.pyx"],
        include_dirs=[numpy.get_include()],
        libraries=libraries,
    )

    config.add_subpackage("tests")

    # HDBSCAN subpackage
    config.add_subpackage("_hdbscan.tests")
    config.add_extension(
        "_hdbscan._hdbscan_boruvka",
        sources=["_hdbscan/_hdbscan_boruvka.pyx"],
        include_dirs=[numpy.get_include(), "_hdbscan"],
        libraries=libraries,
    )
    config.add_extension(
        "_hdbscan._hdbscan_linkage",
        sources=["_hdbscan/_hdbscan_linkage.pyx"],
        include_dirs=[numpy.get_include()],
        libraries=libraries,
    )
    config.add_extension(
        "_hdbscan._hdbscan_reachability",
        sources=["_hdbscan/_hdbscan_reachability.pyx"],
        include_dirs=[numpy.get_include()],
        libraries=libraries,
    )
    config.add_extension(
        "_hdbscan._hdbscan_tree",
        sources=["_hdbscan/_hdbscan_tree.pyx"],
        include_dirs=[numpy.get_include()],
        libraries=libraries,
    )
    config.add_extension(
        "_hdbscan._prediction_utils",
        sources=["_hdbscan/_prediction_utils.pyx"],
        include_dirs=[numpy.get_include()],
        libraries=libraries,
    )
    return config


if __name__ == "__main__":
    from numpy.distutils.core import setup

    setup(**configuration(top_path="").todict())
