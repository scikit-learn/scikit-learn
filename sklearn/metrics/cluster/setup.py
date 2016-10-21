import os

import numpy
from numpy.distutils.misc_util import Configuration

from sklearn._build_utils import add_cython_extension


def configuration(parent_package="", top_path=None):
    config = Configuration("metrics/cluster", parent_package, top_path)
    libraries = []
    if os.name == 'posix':
        libraries.append('m')
    add_cython_extension(top_path,
                         config,
                         "expected_mutual_info_fast",
                         sources=["expected_mutual_info_fast.c"],
                         include_dirs=[numpy.get_include()],
                         libraries=libraries)

    config.add_subpackage("tests")

    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(**configuration().todict())
