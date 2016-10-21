import os
from os.path import join

import numpy

from sklearn._build_utils import get_blas_info
from sklearn._build_utils import add_cython_extension


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration('linear_model', parent_package, top_path)

    cblas_libs, blas_info = get_blas_info()

    if os.name == 'posix':
        cblas_libs.append('m')

    add_cython_extension(top_path,
                         config,
                         'cd_fast', sources=['cd_fast.c'],
                         libraries=cblas_libs,
                         include_dirs=[join('..', 'src', 'cblas'),
                                       numpy.get_include(),
                                       blas_info.pop('include_dirs', [])],
                         extra_compile_args=blas_info.pop('extra_compile_args',
                                                          []), **blas_info)

    add_cython_extension(top_path,
                         config,
                         'sgd_fast',
                         sources=['sgd_fast.c'],
                         include_dirs=[join('..', 'src', 'cblas'),
                                       numpy.get_include(),
                                       blas_info.pop('include_dirs', [])],
                         libraries=cblas_libs,
                         extra_compile_args=blas_info.pop('extra_compile_args',
                                                          []),
                         **blas_info)

    add_cython_extension(top_path,
                         config,
                         'sag_fast',
                         sources=['sag_fast.c'],
                         include_dirs=numpy.get_include())

    # add other directories
    config.add_subpackage('tests')

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
