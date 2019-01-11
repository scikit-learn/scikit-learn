import numpy
from numpy.distutils.misc_util import Configuration


def configuration(parent_package="", top_path=None):
    config = Configuration("ensemble", parent_package, top_path)
    config.add_extension("_gradient_boosting",
                         sources=["_gradient_boosting.pyx"],
                         include_dirs=[numpy.get_include()])

    # config.add_extension("gbm._gradient_boosting",
    #                      sources=["gbm/_gradient_boosting.pyx"],
    #                      include_dirs=[numpy.get_include()],
    #                      extra_compile_args=['-fopenmp'],
    #                      extra_link_args=['-fopenmp'])

    # config.add_extension("gbm.histogram",
    #                      sources=["gbm/histogram.pyx"],
    #                      include_dirs=[numpy.get_include()])

    # config.add_extension("gbm.splitting",
    #                      sources=["gbm/splitting.pyx"],
    #                      include_dirs=[numpy.get_include()])

    # config.add_extension("gbm.binning",
    #                      sources=["gbm/binning.pyx"],
    #                      include_dirs=[numpy.get_include()],
    #                      extra_compile_args=['-fopenmp'],
    #                      extra_link_args=['-fopenmp'])

    # config.add_extension("gbm.predictor",
    #                      sources=["gbm/predictor.pyx"],
    #                      include_dirs=[numpy.get_include()])

    # config.add_extension("gbm.loss",
    #                      sources=["gbm/loss.pyx"],
    #                      include_dirs=[numpy.get_include()],
    #                      extra_compile_args=['-fopenmp'],
    #                      extra_link_args=['-fopenmp'])

    # config.add_extension("gbm.playground",
    #                      sources=["gbm/playground.pyx"],
    #                      include_dirs=[numpy.get_include()])

    config.add_subpackage("tests")
    # config.add_data_files("gbm/histogram.pxd")

    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(**configuration().todict())
