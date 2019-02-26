import numpy
from numpy.distutils.misc_util import Configuration


def configuration(parent_package="", top_path=None):
    config = Configuration("_fast_gradient_boosting", parent_package, top_path)

    config.add_extension("_gradient_boosting",
                         sources=["_gradient_boosting.pyx"],
                         include_dirs=[numpy.get_include()])

    config.add_extension("histogram",
                         sources=["histogram.pyx"],
                         include_dirs=[numpy.get_include()])

    config.add_extension("splitting",
                         sources=["splitting.pyx"],
                         include_dirs=[numpy.get_include()])

    config.add_extension("_binning",
                         sources=["_binning.pyx"],
                         include_dirs=[numpy.get_include()])

    config.add_extension("_predictor",
                         sources=["_predictor.pyx"],
                         include_dirs=[numpy.get_include()])

    config.add_extension("_loss",
                         sources=["_loss.pyx"],
                         include_dirs=[numpy.get_include()])

    config.add_extension("types",
                         sources=["types.pyx"],
                         include_dirs=[numpy.get_include()])

    config.add_extension("utils",
                         sources=["utils.pyx"],
                         include_dirs=[numpy.get_include()])

    config.add_subpackage("tests")

    return config


if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(**configuration().todict())
