cimport numpy as cnp

from ._base cimport BaseDistanceReducer64
from ._gemm_term_computer cimport GEMMTermComputer64

from ...utils._typedefs cimport DTYPE_t

cnp.import_array()

cdef class PairwiseDistances64(BaseDistanceReducer64):
    """64bit implementation of PairwiseDistances."""

    cdef:
        DTYPE_t[:, ::1] pairwise_distances_matrix


cdef class EuclideanPairwiseDistances64(PairwiseDistances64):
    """EuclideanDistance-specialized 64bit implementation for PairwiseDistances."""
    cdef:
        GEMMTermComputer64 gemm_term_computer
        const DTYPE_t[::1] X_norm_squared
        const DTYPE_t[::1] Y_norm_squared

        bint use_squared_distances
