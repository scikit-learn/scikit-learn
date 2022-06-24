from ...utils._typedefs cimport DTYPE_t, ITYPE_t
from libcpp.vector cimport vector


cdef class GEMMTermComputer64:
    cdef:
        const DTYPE_t[:, ::1] X
        const DTYPE_t[:, ::1] Y

        ITYPE_t effective_n_threads
        ITYPE_t chunks_n_threads
        ITYPE_t dist_middle_terms_chunks_size
        ITYPE_t n_features
        ITYPE_t chunk_size

        # Buffers for the `-2 * X_c @ Y_c.T` term computed via GEMM
        vector[vector[DTYPE_t]] dist_middle_terms_chunks

    cdef void _parallel_on_X_pre_compute_and_reduce_distances_on_chunks(
        self,
        ITYPE_t X_start,
        ITYPE_t X_end,
        ITYPE_t Y_start,
        ITYPE_t Y_end,
        ITYPE_t thread_num,
    ) nogil

    cdef void _parallel_on_X_parallel_init(self, ITYPE_t thread_num) nogil

    cdef void _parallel_on_X_init_chunk(
        self,
        ITYPE_t thread_num,
        ITYPE_t X_start,
        ITYPE_t X_end,
    ) nogil

    cdef void _parallel_on_Y_init(self) nogil

    cdef void _parallel_on_Y_parallel_init(
        self,
        ITYPE_t thread_num,
        ITYPE_t X_start,
        ITYPE_t X_end,
    ) nogil

    cdef void _parallel_on_Y_pre_compute_and_reduce_distances_on_chunks(
        self,
        ITYPE_t X_start,
        ITYPE_t X_end,
        ITYPE_t Y_start,
        ITYPE_t Y_end,
        ITYPE_t thread_num
    ) nogil

    cdef DTYPE_t * _compute_distances_on_chunks(
        self,
        ITYPE_t X_start,
        ITYPE_t X_end,
        ITYPE_t Y_start,
        ITYPE_t Y_end,
        ITYPE_t thread_num,
    ) nogil
