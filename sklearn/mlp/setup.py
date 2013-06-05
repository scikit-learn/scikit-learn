import os
from os.path import join
import numpy

from sklearn._build_utils import get_blas_info


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    #from numpy.distutils.system_info import get_info
    config = Configuration('mlp', parent_package, top_path)

    cblas_libs, blas_info = get_blas_info()

    if os.name == 'posix':
        cblas_libs.append('m')

    config.add_extension('mlp_fast',
         sources=['mlp_fast.c'],
         libraries=cblas_libs,
         include_dirs=[join('..', 'src', 'cblas'),
             numpy.get_include(),
             blas_info.pop('include_dirs', [])],
         extra_compile_args=blas_info.pop('extra_compile_args', []),
         **blas_info)

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
