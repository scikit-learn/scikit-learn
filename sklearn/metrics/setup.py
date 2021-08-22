import os

from numpy.distutils.misc_util import Configuration


def configuration(parent_package="", top_path=None):
    config = Configuration("metrics", parent_package, top_path)

    libraries = []
    if os.name == "posix":
        libraries.append("m")

    config.add_subpackage("_plot")
    config.add_subpackage("_plot.tests")
    config.add_subpackage("cluster")

    config.add_extension(
        "_pairwise_fast", sources=["_pairwise_fast.pyx"], libraries=libraries
    )

    config.add_subpackage("tests")

    return config


if __name__ == "__main__":
    from numpy.distutils.core import setup

    setup(**configuration().todict())
