"""
Utilities useful during the build.
"""
# author: Andy Mueller, Gael Varoquaux
# license: BSD

from __future__ import division, print_function, absolute_import

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
