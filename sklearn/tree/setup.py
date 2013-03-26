import os

import numpy
from numpy.distutils.misc_util import Configuration, dot_join
from numpy.distutils.core import Extension
from Cython.Build import cythonize


def configuration(parent_package="", top_path=None):
    config = Configuration("tree", parent_package, top_path)
    libraries = []
    if os.name == 'posix':
        libraries.append('m')

    extensions = [
        Extension(
            dot_join(config.name, '_tree'),
            sources=config.paths(['_tree.pyx']),
            include_dirs=[numpy.get_include()],
            libraries=libraries)
    ]
    extensions = cythonize(extensions)

    config.ext_modules.extend(extensions)

    config.add_subpackage("tests")

    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(**configuration().todict())
