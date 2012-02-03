from os.path import join
import numpy

import warnings


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info, BlasNotFoundError

    config = Configuration('svm', parent_package, top_path)

    config.add_subpackage('tests')

    # Section LibSVM

    # we compile both libsvm and lisvm_sparse
    config.add_library('libsvm-skl',
                       sources=[join('src', 'libsvm', 'libsvm_template.cpp')],
                       depends=[join('src', 'libsvm', 'svm.cpp'),
                                join('src', 'libsvm', 'svm.h')],
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
    blas_sources = [join('src', 'blas', 'daxpy.c'),
                    join('src', 'blas', 'ddot.c'),
                    join('src', 'blas', 'dnrm2.c'),
                    join('src', 'blas', 'dscal.c')]

    liblinear_sources = ['liblinear.c',
                         join('src', 'liblinear', '*.cpp')]

    liblinear_depends = [join('src', 'liblinear', '*.h'),
                         join('src', 'liblinear', 'liblinear_helper.c')]

    # we try to link agains system-wide blas
    blas_info = get_info('blas_opt', 0)

    if not blas_info:
        config.add_library('blas', blas_sources)
        warnings.warn(BlasNotFoundError.__doc__)

    config.add_extension('liblinear',
                         sources=liblinear_sources,
                         libraries=blas_info.pop('libraries', ['blas']),
                         include_dirs=['src',
                                       numpy.get_include(),
                                       blas_info.pop('include_dirs', [])],
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

    config.add_subpackage('sparse')

    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
