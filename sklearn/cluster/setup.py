# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
# License: BSD Style.
import os
from os.path import join

import numpy


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info

    blas_info = get_info('blas_opt', 0)
    if (not blas_info) or (
        ('NO_ATLAS_INFO', 1) in blas_info.get('define_macros', [])):
        cblas_libs = ['cblas']
        blas_info.pop('libraries', None)
    else:
        cblas_libs = blas_info.pop('libraries', [])
    libraries = []
    if os.name == 'posix':
        cblas_libs.append('m')
        libraries.append('m')

    config = Configuration('cluster', parent_package, top_path)
    config.add_extension('_hierarchical',
                         sources=['_hierarchical.c'],
                         include_dirs=[numpy.get_include()],
                         libraries=libraries)

    config.add_extension(
        '_k_means',
        libraries=cblas_libs,
        sources=['_k_means.c'],
        include_dirs=[join('..', 'src', 'cblas'),
                      numpy.get_include(),
                      blas_info.pop('include_dirs', [])],
        extra_compile_args=blas_info.pop('extra_compile_args', []),
        **blas_info
    )
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
