import numpy
from numpy.distutils.misc_util import Configuration


def configuration(parent_package="", top_path=None):
    config = Configuration("gbm", parent_package, top_path)

    config.add_extension("_gradient_boosting",
                         sources=["_gradient_boosting.pyx"],
                         include_dirs=[numpy.get_include()],
                         extra_compile_args=['-fopenmp'],
                         extra_link_args=['-fopenmp'])

    config.add_extension("histogram",
                         sources=["histogram.pyx"],
                         include_dirs=[numpy.get_include()])

    config.add_extension("splitting",
                         sources=["splitting.pyx"],
                         include_dirs=[numpy.get_include()],
                         extra_compile_args=['-fopenmp'],
                         extra_link_args=['-fopenmp'])

    config.add_extension("binning",
                         sources=["binning.pyx"],
                         include_dirs=[numpy.get_include()],
                         extra_compile_args=['-fopenmp'],
                         extra_link_args=['-fopenmp'])

    config.add_extension("predictor",
                         sources=["predictor.pyx"],
                         include_dirs=[numpy.get_include()],
                         extra_compile_args=['-fopenmp'],
                         extra_link_args=['-fopenmp'])

    config.add_extension("loss",
                         sources=["loss.pyx"],
                         include_dirs=[numpy.get_include()],
                         extra_compile_args=['-fopenmp'],
                         extra_link_args=['-fopenmp'])

    config.add_extension("types",
                         sources=["types.pyx"],
                         include_dirs=[numpy.get_include()])

    config.add_extension("playground",
                         sources=["playground.pyx"],
                         include_dirs=[numpy.get_include()],
                         extra_compile_args=['-fopenmp'],
                         extra_link_args=['-fopenmp'])

    config.add_subpackage("tests")
    # config.add_data_files("histogram.pxd")

    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(**configuration().todict())
