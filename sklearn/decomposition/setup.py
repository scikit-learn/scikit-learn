import os
from os.path import join

import numpy
from numpy.distutils.misc_util import Configuration

from sklearn._build_utils import get_blas_info


def configuration(parent_package="", top_path=None):
    config = Configuration("decomposition", parent_package, top_path)

    cblas_libs, blas_info = get_blas_info()

    libraries = []
    if os.name == 'posix':
        libraries.append('m')

    config.add_extension("_online_lda",
                         sources=["_online_lda.pyx"],
                         include_dirs=[numpy.get_include()],
                         libraries=libraries)

    config.add_extension('cdnmf_fast',
                         sources=['cdnmf_fast.pyx'],
                         include_dirs=[numpy.get_include()],
                         libraries=libraries)

    config.add_extension('_update_dict_fast',
                         sources=['_update_dict_fast.pyx'],
                         include_dirs=[join('..', 'src', 'cblas'),
                                       numpy.get_include(),
                                       blas_info.pop('include_dirs', [])],
                         libraries=libraries + cblas_libs,
                         extra_compile_args=blas_info.pop('extra_compile_args',
                                                          []), **blas_info)

    config.add_subpackage("tests")

    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(**configuration().todict())
