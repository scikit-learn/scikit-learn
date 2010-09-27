from os.path import join
import warnings
import numpy
import sys
if sys.version_info[0] < 3:
    from ConfigParser import ConfigParser
else:
    from configparser import ConfigParser

def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info, get_standard_file, \
        BlasNotFoundError
    config = Configuration('learn', parent_package, top_path)

    site_cfg  = ConfigParser()
    site_cfg.read(get_standard_file('site.cfg'))

    config.add_subpackage('datasets')
    config.add_subpackage('features')
    config.add_subpackage('features/tests')
    config.add_subpackage('cluster')
    config.add_subpackage('cluster/tests')
    config.add_subpackage('covariance')
    config.add_subpackage('covariance/tests')
    config.add_subpackage('feature_selection')
    config.add_subpackage('feature_selection/tests')
    config.add_subpackage('sparse')
    config.add_subpackage('sparse/tests/')
    config.add_subpackage('utils')
    config.add_subpackage('utils/tests')
    config.add_subpackage('externals')

    # Section LibSVM
    libsvm_includes = [numpy.get_include()]
    libsvm_libraries = []
    libsvm_library_dirs = []
    libsvm_sources = [join('src', 'libsvm', '_libsvm.c')]

    # we try to link against system-wide libsvm
    if site_cfg.has_section('libsvm'):
        libsvm_includes.append(site_cfg.get('libsvm', 'include_dirs'))
        libsvm_libraries.append(site_cfg.get('libsvm', 'libraries'))
        libsvm_library_dirs.append(site_cfg.get('libsvm', 'library_dirs'))
    else:
        # if not specified, we build our own libsvm
        libsvm_sources.append(join('src', 'libsvm', 'svm.cpp'))

    config.add_extension('_libsvm',
                         sources=libsvm_sources,
                         include_dirs=libsvm_includes,
                         libraries=libsvm_libraries,
                         library_dirs=libsvm_library_dirs,
                         depends=[join('src', 'libsvm', 'svm.h'),
                                  join('src', 'libsvm', 'libsvm_helper.c')],
                         # add this for gdb debug
                         extra_compile_args=['-O0 -fno-inline']
                                  )

    ### liblinear module
    blas_sources = [join('src', 'blas', 'daxpy.c'),
                    join('src', 'blas', 'ddot.c'),
                    join('src', 'blas', 'dnrm2.c'),
                    join('src', 'blas', 'dscal.c')]

    liblinear_sources = [join('src', 'liblinear', '_liblinear.c'),
                         join('src', 'liblinear', '*.cpp')]

    # we try to link agains system-wide blas
    blas_info = get_info('blas_opt', 0)

    if not blas_info:
        config.add_library('blas', blas_sources)
        warnings.warn(BlasNotFoundError.__doc__)

    config.add_extension('_liblinear',
                         sources=liblinear_sources,
                         libraries = blas_info.pop('libraries', ['blas']),
                         include_dirs=['src',
                                       numpy.get_include(),
                                       blas_info.pop('include_dirs', [])],
                         depends=[join('src', 'liblinear', '*.h')],
                         **blas_info)

    ## end liblinear module

    # some libs needs cblas, fortran-compiled BLAS will not be sufficient
    blas_info = get_info('blas_opt', 0)
    if (not blas_info) or (
        ('NO_ATLAS_INFO', 1) in blas_info.get('define_macros', [])) :
        config.add_library('cblas',
                           sources=[
                               join('src', 'cblas', '*.c'),
                               ]
                           )
        cblas_libs = ['cblas']
        blas_info.pop('libraries', None)
    else:
        cblas_libs = blas_info.pop('libraries', [])

    config.add_extension('ball_tree',
                         sources=[join('src', 'BallTree.cpp')],
                         include_dirs=[numpy.get_include()]
                         )


    # this has to be build *after* cblas
    config.add_subpackage('glm')
    config.add_subpackage('utils')

    # add the test directory
    config.add_subpackage('tests')

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
