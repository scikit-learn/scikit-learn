from __future__ import division, absolute_import, print_function

import warnings

import numpy as np
import numpy.distutils.system_info

from numpy.distutils.system_info import (system_info,
                                         numpy_info,
                                         NotFoundError,
                                         BlasNotFoundError,
                                         LapackNotFoundError,
                                         AtlasNotFoundError,
                                         LapackSrcNotFoundError,
                                         BlasSrcNotFoundError,
                                         dict_append,
                                         get_info as old_get_info)

from scipy._lib._version import NumpyVersion


if NumpyVersion(np.__version__) >= "1.15.0.dev":
    # For new enough numpy.distutils, the ACCELERATE=None environment
    # variable in the top-level setup.py is enough, so no need to
    # customize BLAS detection.
    get_info = old_get_info
else:
    # For numpy < 1.15.0, we need overrides.

    def get_info(name, notfound_action=0):
        # Special case our custom *_opt_info
        cls = {'lapack_opt': lapack_opt_info,
               'blas_opt': blas_opt_info}.get(name.lower())
        if cls is None:
            return old_get_info(name, notfound_action)
        return cls().get_info(notfound_action)

    #
    # The following is copypaste from numpy.distutils.system_info, with
    # OSX Accelerate-related parts removed.
    #

    class lapack_opt_info(system_info):

        notfounderror = LapackNotFoundError

        def calc_info(self):

            lapack_mkl_info = get_info('lapack_mkl')
            if lapack_mkl_info:
                self.set_info(**lapack_mkl_info)
                return

            openblas_info = get_info('openblas_lapack')
            if openblas_info:
                self.set_info(**openblas_info)
                return

            openblas_info = get_info('openblas_clapack')
            if openblas_info:
                self.set_info(**openblas_info)
                return

            atlas_info = get_info('atlas_3_10_threads')
            if not atlas_info:
                atlas_info = get_info('atlas_3_10')
            if not atlas_info:
                atlas_info = get_info('atlas_threads')
            if not atlas_info:
                atlas_info = get_info('atlas')

            need_lapack = 0
            need_blas = 0
            info = {}
            if atlas_info:
                l = atlas_info.get('define_macros', [])
                if ('ATLAS_WITH_LAPACK_ATLAS', None) in l \
                       or ('ATLAS_WITHOUT_LAPACK', None) in l:
                    need_lapack = 1
                info = atlas_info

            else:
                warnings.warn(AtlasNotFoundError.__doc__, stacklevel=2)
                need_blas = 1
                need_lapack = 1
                dict_append(info, define_macros=[('NO_ATLAS_INFO', 1)])

            if need_lapack:
                lapack_info = get_info('lapack')
                #lapack_info = {} ## uncomment for testing
                if lapack_info:
                    dict_append(info, **lapack_info)
                else:
                    warnings.warn(LapackNotFoundError.__doc__, stacklevel=2)
                    lapack_src_info = get_info('lapack_src')
                    if not lapack_src_info:
                        warnings.warn(LapackSrcNotFoundError.__doc__, stacklevel=2)
                        return
                    dict_append(info, libraries=[('flapack_src', lapack_src_info)])

            if need_blas:
                blas_info = get_info('blas')
                if blas_info:
                    dict_append(info, **blas_info)
                else:
                    warnings.warn(BlasNotFoundError.__doc__, stacklevel=2)
                    blas_src_info = get_info('blas_src')
                    if not blas_src_info:
                        warnings.warn(BlasSrcNotFoundError.__doc__, stacklevel=2)
                        return
                    dict_append(info, libraries=[('fblas_src', blas_src_info)])

            self.set_info(**info)
            return

    class blas_opt_info(system_info):

        notfounderror = BlasNotFoundError

        def calc_info(self):

            blas_mkl_info = get_info('blas_mkl')
            if blas_mkl_info:
                self.set_info(**blas_mkl_info)
                return

            blis_info = get_info('blis')
            if blis_info:
                self.set_info(**blis_info)
                return

            openblas_info = get_info('openblas')
            if openblas_info:
                self.set_info(**openblas_info)
                return

            atlas_info = get_info('atlas_3_10_blas_threads')
            if not atlas_info:
                atlas_info = get_info('atlas_3_10_blas')
            if not atlas_info:
                atlas_info = get_info('atlas_blas_threads')
            if not atlas_info:
                atlas_info = get_info('atlas_blas')

            need_blas = 0
            info = {}
            if atlas_info:
                info = atlas_info
            else:
                warnings.warn(AtlasNotFoundError.__doc__, stacklevel=2)
                need_blas = 1
                dict_append(info, define_macros=[('NO_ATLAS_INFO', 1)])

            if need_blas:
                blas_info = get_info('blas')
                if blas_info:
                    dict_append(info, **blas_info)
                else:
                    warnings.warn(BlasNotFoundError.__doc__, stacklevel=2)
                    blas_src_info = get_info('blas_src')
                    if not blas_src_info:
                        warnings.warn(BlasSrcNotFoundError.__doc__, stacklevel=2)
                        return
                    dict_append(info, libraries=[('fblas_src', blas_src_info)])

            self.set_info(**info)
            return
