cimport numpy as cnp

from ._base cimport (
    PairwiseDistancesReduction64,
)
from ._gemm_term_computer cimport GEMMTermComputer64

cnp.import_array()

cdef class PairwiseDistancesArgKmin64(PairwiseDistancesReduction64):
    """64bit implementation of PairwiseDistancesArgKmin."""

    cdef:
        cnp.intp_t k

        cnp.intp_t[:, ::1] argkmin_indices
        cnp.float64_t[:, ::1] argkmin_distances

        # Used as array of pointers to private datastructures used in threads.
        cnp.float64_t ** heaps_r_distances_chunks
        cnp.intp_t ** heaps_indices_chunks


cdef class FastEuclideanPairwiseDistancesArgKmin64(PairwiseDistancesArgKmin64):
    """EuclideanDistance-specialized 64bit implementation for PairwiseDistancesArgKmin."""
    cdef:
        GEMMTermComputer64 gemm_term_computer
        const cnp.float64_t[::1] X_norm_squared
        const cnp.float64_t[::1] Y_norm_squared

        bint use_squared_distances
