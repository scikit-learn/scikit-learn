import numpy
from numpy.distutils.misc_util import Configuration


def configuration(parent_package="", top_path=None):
    config = Configuration("earth", parent_package, top_path)
    config.add_extension("_basis",
                         sources=["_basis.pyx"],
                         include_dirs=[numpy.get_include()])
    config.add_extension("_forward",
                         sources=["_forward.pyx"],
                         include_dirs=[numpy.get_include()])
    config.add_extension("_pruning",
                         sources=["_pruning.pyx"],
                         include_dirs=[numpy.get_include()])
    config.add_extension("_record",
                         sources=["_record.pyx"],
                         include_dirs=[numpy.get_include()])
    config.add_extension("_util",
                         sources=["_util.pyx"],
                         include_dirs=[numpy.get_include()])
    config.add_extension("_types",
                         sources=["_types.pyx"],
                         include_dirs=[numpy.get_include()])
    config.add_subpackage("tests")

    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(**configuration().todict())
