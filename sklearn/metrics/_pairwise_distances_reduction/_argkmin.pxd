
cimport numpy as cnp
from ...utils._typedefs cimport ITYPE_t, DTYPE_t

cnp.import_array()

from ._base cimport PairwiseDistancesReduction64
from ._gemm_term_computer cimport GEMMTermComputer64

cdef class PairwiseDistancesArgKmin64(PairwiseDistancesReduction64):
    """64bit implementation of PairwiseDistancesArgKmin."""

    cdef:
        ITYPE_t k

        ITYPE_t[:, ::1] argkmin_indices
        DTYPE_t[:, ::1] argkmin_distances

        # Used as array of pointers to private datastructures used in threads.
        DTYPE_t ** heaps_r_distances_chunks
        ITYPE_t ** heaps_indices_chunks


cdef class FastEuclideanPairwiseDistancesArgKmin64(PairwiseDistancesArgKmin64):
    """EuclideanDistance-specialized 64bit implementation for PairwiseDistancesArgKmin."""
    cdef:
        GEMMTermComputer64 gemm_term_computer
        const DTYPE_t[::1] X_norm_squared
        const DTYPE_t[::1] Y_norm_squared

        bint use_squared_distances

from ._base cimport PairwiseDistancesReduction32
from ._gemm_term_computer cimport GEMMTermComputer32

cdef class PairwiseDistancesArgKmin32(PairwiseDistancesReduction32):
    """32bit implementation of PairwiseDistancesArgKmin."""

    cdef:
        ITYPE_t k

        ITYPE_t[:, ::1] argkmin_indices
        DTYPE_t[:, ::1] argkmin_distances

        # Used as array of pointers to private datastructures used in threads.
        DTYPE_t ** heaps_r_distances_chunks
        ITYPE_t ** heaps_indices_chunks


cdef class FastEuclideanPairwiseDistancesArgKmin32(PairwiseDistancesArgKmin32):
    """EuclideanDistance-specialized 32bit implementation for PairwiseDistancesArgKmin."""
    cdef:
        GEMMTermComputer32 gemm_term_computer
        const DTYPE_t[::1] X_norm_squared
        const DTYPE_t[::1] Y_norm_squared

        bint use_squared_distances
