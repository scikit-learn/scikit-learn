import os
from os.path import join
import numpy

from sklearn._build_utils import get_blas_info


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration('svm', parent_package, top_path)

    config.add_subpackage('tests')

    # Section LibSVM

    # we compile both libsvm and libsvm_sparse
    config.add_library('libsvm-skl',
                       sources=[join('src', 'libsvm', 'libsvm_template.cpp')],
                       depends=[join('src', 'libsvm', 'svm.cpp'),
                                join('src', 'libsvm', 'svm.h')],
                       # Force C++ linking in case gcc is picked up instead
                       # of g++ under windows with some versions of MinGW
                       extra_link_args=['-lstdc++'],
                       )

    libsvm_sources = ['libsvm.c']
    libsvm_depends = [join('src', 'libsvm', 'libsvm_helper.c'),
                      join('src', 'libsvm', 'libsvm_template.cpp'),
                      join('src', 'libsvm', 'svm.cpp'),
                      join('src', 'libsvm', 'svm.h')]

    config.add_extension('libsvm',
                         sources=libsvm_sources,
                         include_dirs=[numpy.get_include(),
                                       join('src', 'libsvm')],
                         libraries=['libsvm-skl'],
                         depends=libsvm_depends,
                         )

    ### liblinear module
    cblas_libs, blas_info = get_blas_info()
    if os.name == 'posix':
        cblas_libs.append('m')

    liblinear_sources = ['liblinear.c',
                         join('src', 'liblinear', '*.cpp')]

    liblinear_depends = [join('src', 'liblinear', '*.h'),
                         join('src', 'liblinear', 'liblinear_helper.c')]

    config.add_extension('liblinear',
                         sources=liblinear_sources,
                         libraries=cblas_libs,
                         include_dirs=[join('..', 'src', 'cblas'),
                                       numpy.get_include(),
                                       blas_info.pop('include_dirs', [])],
                         extra_compile_args=blas_info.pop('extra_compile_args',
                                                          []),
                         depends=liblinear_depends,
                         # extra_compile_args=['-O0 -fno-inline'],
                         ** blas_info)

    ## end liblinear module

    # this should go *after* libsvm-skl
    libsvm_sparse_sources = ['libsvm_sparse.c']
    config.add_extension('libsvm_sparse', libraries=['libsvm-skl'],
                         sources=libsvm_sparse_sources,
                         include_dirs=[numpy.get_include(),
                                       join("src", "libsvm")],
                         depends=[join("src", "libsvm", "svm.h"),
                                  join("src", "libsvm",
                                       "libsvm_sparse_helper.c")])

    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
