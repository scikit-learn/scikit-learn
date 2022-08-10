
cimport numpy as cnp
from ...utils._typedefs cimport ITYPE_t, DTYPE_t

cnp.import_array()

from ._base cimport ComputationTemplate64
from ._gemm_term_computer cimport GEMMTermComputer64

cdef class ArgKmin64(ComputationTemplate64):
    """64bit implementation of ComputationTemplate64 for the `ArgKmin` reduction."""

    cdef:
        ITYPE_t k

        ITYPE_t[:, ::1] argkmin_indices
        DTYPE_t[:, ::1] argkmin_distances

        # Used as array of pointers to private datastructures used in threads.
        DTYPE_t ** heaps_r_distances_chunks
        ITYPE_t ** heaps_indices_chunks


cdef class EuclideanArgKmin64(ArgKmin64):
    """EuclideanDistance-specialized 64bit implementation of ArgKmin64."""
    cdef:
        GEMMTermComputer64 gemm_term_computer
        const DTYPE_t[::1] X_norm_squared
        const DTYPE_t[::1] Y_norm_squared

        bint use_squared_distances

from ._base cimport ComputationTemplate32
from ._gemm_term_computer cimport GEMMTermComputer32

cdef class ArgKmin32(ComputationTemplate32):
    """32bit implementation of ComputationTemplate32 for the `ArgKmin` reduction."""

    cdef:
        ITYPE_t k

        ITYPE_t[:, ::1] argkmin_indices
        DTYPE_t[:, ::1] argkmin_distances

        # Used as array of pointers to private datastructures used in threads.
        DTYPE_t ** heaps_r_distances_chunks
        ITYPE_t ** heaps_indices_chunks


cdef class EuclideanArgKmin32(ArgKmin32):
    """EuclideanDistance-specialized 32bit implementation of ArgKmin32."""
    cdef:
        GEMMTermComputer32 gemm_term_computer
        const DTYPE_t[::1] X_norm_squared
        const DTYPE_t[::1] Y_norm_squared

        bint use_squared_distances
