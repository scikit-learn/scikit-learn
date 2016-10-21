"""
Utilities useful during the build.
"""
# author: Andy Mueller, Gael Varoquaux
# license: BSD

from __future__ import division, print_function, absolute_import

import os

DEFAULT_ROOT = 'sklearn'

from numpy.distutils.system_info import get_info


def get_blas_info():
    def atlas_not_found(blas_info_):
        def_macros = blas_info.get('define_macros', [])
        for x in def_macros:
            if x[0] == "NO_ATLAS_INFO":
                # if x[1] != 1 we should have lapack
                # how do we do that now?
                return True
            if x[0] == "ATLAS_INFO":
                if "None" in x[1]:
                    # this one turned up on FreeBSD
                    return True
        return False

    blas_info = get_info('blas_opt', 0)
    if (not blas_info) or atlas_not_found(blas_info):
        cblas_libs = ['cblas']
        blas_info.pop('libraries', None)
    else:
        cblas_libs = blas_info.pop('libraries', [])

    return cblas_libs, blas_info


def get_cython_source(filename):
    is_c_filename = filename.endswith('.c')
    is_cpp_filename = filename.endswith('.cpp')

    # files in src are .c and .cpp files that are not cython-generated
    if 'src/' in filename and (is_c_filename or is_cpp_filename):
        return filename
    elif is_c_filename:
        filename = filename[:-1]
    elif is_cpp_filename:
        filename = filename[:-3]
    else:
        raise ValueError('Only .c and .cpp files are supported. '
                         'Got {0!r} instead'.format(filename))
    return filename + 'pyx'


def add_cython_extension(top_path, config, name, sources, **kwargs):
    is_dev_version = not os.path.exists(os.path.join(top_path, 'PKG-INFO'))

    if is_dev_version:
        try:
            from Cython.Build import cythonize
        except ImportError:
            raise ValueError('Please install cython in order '
                             'to build a scikit-learn development version')

        sources = [get_cython_source(filename) for filename in sources]

    config.add_extension(name, sources, **kwargs)

    if is_dev_version:
        config.ext_modules = cythonize(config.ext_modules)
