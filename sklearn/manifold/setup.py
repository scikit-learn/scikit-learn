import os
from os.path import join

import numpy
from numpy.distutils.misc_util import Configuration
from sklearn._build_utils import get_blas_info


def configuration(parent_package="", top_path=None):
    config = Configuration("manifold", parent_package, top_path)
    libraries = []
    if os.name == 'posix':
        libraries.append('m')
    config.add_extension("_utils",
                         sources=["_utils.pyx"],
                         include_dirs=[numpy.get_include()],
                         libraries=libraries,
                         extra_compile_args=["-O3"])
    cblas_libs, blas_info = get_blas_info()
    eca = blas_info.pop('extra_compile_args', [])
    eca.append("-O4")
    config.add_extension("_barnes_hut_tsne",
                         libraries=cblas_libs,
                         sources=["_barnes_hut_tsne.pyx"],
                         include_dirs=[join('..', 'src', 'cblas'),
                                       numpy.get_include(),
                                       blas_info.pop('include_dirs', [])],
                         extra_compile_args=eca, **blas_info)

    config.add_subpackage('tests')

    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(**configuration().todict())
