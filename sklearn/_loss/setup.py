import numpy
from numpy.distutils.misc_util import Configuration


def configuration(parent_package="", top_path=None):
    config = Configuration("_loss", parent_package, top_path)
    config.add_extension(
        "_loss",
        sources=["_loss.pyx"],
        include_dirs=[numpy.get_include()],
        define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_13_API_VERSION")],
    )
    return config


if __name__ == "__main__":
    from numpy.distutils.core import setup

    setup(**configuration().todict())
