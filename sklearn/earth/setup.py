import numpy
from numpy.distutils.misc_util import Configuration

def configuration(parent_package="", top_path=None):
    config = Configuration("earth", parent_package, top_path)
    config.add_extension("_basis",
                         sources=["_basis.c"],
                         include_dirs=[numpy.get_include()])
    config.add_extension("_forward",
                         sources=["_forward.c"],
                         include_dirs=[numpy.get_include()])
    config.add_extension("_pruning",
                         sources=["_pruning.c"],
                         include_dirs=[numpy.get_include()])
    config.add_extension("_record",
                         sources=["_record.c"],
                         include_dirs=[numpy.get_include()])
    config.add_extension("_util",
                         sources=["_util.c"],
                         include_dirs=[numpy.get_include()])
    config.add_subpackage("tests")

    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(**configuration().todict())
