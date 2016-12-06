# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
# License: BSD 3 clause
import os
from os.path import join

import numpy

from sklearn._build_utils import get_blas_info


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration

    cblas_libs, blas_info = get_blas_info()

    libraries = []
    if os.name == 'posix':
        cblas_libs.append('m')
        libraries.append('m')

    config = Configuration('cluster', parent_package, top_path)
    config.add_extension('_dbscan_inner',
                         sources=['_dbscan_inner.pyx'],
                         include_dirs=[numpy.get_include()],
                         language="c++")
    config.add_extension('_optics_inner',
                         sources=['_optics_inner.pyx'],
                         include_dirs=[numpy.get_include()],
                         libraries=libraries)

    config.add_extension('_hierarchical',
                         sources=['_hierarchical.pyx'],
                         language="c++",
                         include_dirs=[numpy.get_include()],
                         libraries=libraries)
    config.add_extension('_k_means_elkan',
                         sources=['_k_means_elkan.pyx'],
                         include_dirs=[numpy.get_include()],
                         libraries=libraries)

    config.add_extension('_k_means',
                         libraries=cblas_libs,
                         sources=['_k_means.pyx'],
                         include_dirs=[join('..', 'src', 'cblas'),
                                       numpy.get_include(),
                                       blas_info.pop('include_dirs', [])],
                         extra_compile_args=blas_info.pop(
                             'extra_compile_args', []),
                         **blas_info
                         )

    config.add_subpackage('tests')

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
