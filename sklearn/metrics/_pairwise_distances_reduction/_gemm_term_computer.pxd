cimport numpy as cnp

from ...utils._typedefs cimport DTYPE_t, ITYPE_t, SPARSE_INDEX_TYPE_t
from libcpp.vector cimport vector


cdef void _middle_term_sparse_sparse_64(
    const DTYPE_t[:] X_data,
    const SPARSE_INDEX_TYPE_t[:] X_indices,
    const SPARSE_INDEX_TYPE_t[:] X_indptr,
    ITYPE_t X_start,
    ITYPE_t X_end,
    const DTYPE_t[:] Y_data,
    const SPARSE_INDEX_TYPE_t[:] Y_indices,
    const SPARSE_INDEX_TYPE_t[:] Y_indptr,
    ITYPE_t Y_start,
    ITYPE_t Y_end,
    DTYPE_t * D,
) nogil


cdef class MiddleTermComputer64:
    cdef:
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

    cdef DTYPE_t * _compute_dist_middle_terms(
        self,
        ITYPE_t X_start,
        ITYPE_t X_end,
        ITYPE_t Y_start,
        ITYPE_t Y_end,
        ITYPE_t thread_num,
    ) nogil


cdef class DenseDenseMiddleTermComputer64(MiddleTermComputer64):
    cdef:
        const DTYPE_t[:, ::1] X
        const DTYPE_t[:, ::1] Y


    cdef void _parallel_on_X_pre_compute_and_reduce_distances_on_chunks(
        self,
        ITYPE_t X_start,
        ITYPE_t X_end,
        ITYPE_t Y_start,
        ITYPE_t Y_end,
        ITYPE_t thread_num,
    ) nogil

    cdef void _parallel_on_X_init_chunk(
        self,
        ITYPE_t thread_num,
        ITYPE_t X_start,
        ITYPE_t X_end,
    ) nogil

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

    cdef DTYPE_t * _compute_dist_middle_terms(
        self,
        ITYPE_t X_start,
        ITYPE_t X_end,
        ITYPE_t Y_start,
        ITYPE_t Y_end,
        ITYPE_t thread_num,
    ) nogil


cdef class SparseSparseMiddleTermComputer64(MiddleTermComputer64):
    cdef:
        const DTYPE_t[:] X_data
        const SPARSE_INDEX_TYPE_t[:] X_indices
        const SPARSE_INDEX_TYPE_t[:] X_indptr

        const DTYPE_t[:] Y_data
        const SPARSE_INDEX_TYPE_t[:] Y_indices
        const SPARSE_INDEX_TYPE_t[:] Y_indptr

    cdef void _parallel_on_X_pre_compute_and_reduce_distances_on_chunks(
        self,
        ITYPE_t X_start,
        ITYPE_t X_end,
        ITYPE_t Y_start,
        ITYPE_t Y_end,
        ITYPE_t thread_num
    ) nogil

    cdef void _parallel_on_Y_pre_compute_and_reduce_distances_on_chunks(
        self,
        ITYPE_t X_start,
        ITYPE_t X_end,
        ITYPE_t Y_start,
        ITYPE_t Y_end,
        ITYPE_t thread_num
    ) nogil

    cdef DTYPE_t * _compute_dist_middle_terms(
        self,
        ITYPE_t X_start,
        ITYPE_t X_end,
        ITYPE_t Y_start,
        ITYPE_t Y_end,
        ITYPE_t thread_num,
    ) nogil


cdef void _middle_term_sparse_sparse_32(
    const cnp.float32_t[:] X_data,
    const SPARSE_INDEX_TYPE_t[:] X_indices,
    const SPARSE_INDEX_TYPE_t[:] X_indptr,
    ITYPE_t X_start,
    ITYPE_t X_end,
    const cnp.float32_t[:] Y_data,
    const SPARSE_INDEX_TYPE_t[:] Y_indices,
    const SPARSE_INDEX_TYPE_t[:] Y_indptr,
    ITYPE_t Y_start,
    ITYPE_t Y_end,
    DTYPE_t * D,
) nogil


cdef class MiddleTermComputer32:
    cdef:
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

    cdef DTYPE_t * _compute_dist_middle_terms(
        self,
        ITYPE_t X_start,
        ITYPE_t X_end,
        ITYPE_t Y_start,
        ITYPE_t Y_end,
        ITYPE_t thread_num,
    ) nogil


cdef class DenseDenseMiddleTermComputer32(MiddleTermComputer32):
    cdef:
        const cnp.float32_t[:, ::1] X
        const cnp.float32_t[:, ::1] Y

        # Buffers for upcasting chunks of X and Y from 32bit to 64bit
        vector[vector[DTYPE_t]] X_c_upcast
        vector[vector[DTYPE_t]] Y_c_upcast

    cdef void _parallel_on_X_pre_compute_and_reduce_distances_on_chunks(
        self,
        ITYPE_t X_start,
        ITYPE_t X_end,
        ITYPE_t Y_start,
        ITYPE_t Y_end,
        ITYPE_t thread_num,
    ) nogil

    cdef void _parallel_on_X_init_chunk(
        self,
        ITYPE_t thread_num,
        ITYPE_t X_start,
        ITYPE_t X_end,
    ) nogil

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

    cdef DTYPE_t * _compute_dist_middle_terms(
        self,
        ITYPE_t X_start,
        ITYPE_t X_end,
        ITYPE_t Y_start,
        ITYPE_t Y_end,
        ITYPE_t thread_num,
    ) nogil


cdef class SparseSparseMiddleTermComputer32(MiddleTermComputer32):
    cdef:
        const cnp.float32_t[:] X_data
        const SPARSE_INDEX_TYPE_t[:] X_indices
        const SPARSE_INDEX_TYPE_t[:] X_indptr

        const cnp.float32_t[:] Y_data
        const SPARSE_INDEX_TYPE_t[:] Y_indices
        const SPARSE_INDEX_TYPE_t[:] Y_indptr

    cdef void _parallel_on_X_pre_compute_and_reduce_distances_on_chunks(
        self,
        ITYPE_t X_start,
        ITYPE_t X_end,
        ITYPE_t Y_start,
        ITYPE_t Y_end,
        ITYPE_t thread_num
    ) nogil

    cdef void _parallel_on_Y_pre_compute_and_reduce_distances_on_chunks(
        self,
        ITYPE_t X_start,
        ITYPE_t X_end,
        ITYPE_t Y_start,
        ITYPE_t Y_end,
        ITYPE_t thread_num
    ) nogil

    cdef DTYPE_t * _compute_dist_middle_terms(
        self,
        ITYPE_t X_start,
        ITYPE_t X_end,
        ITYPE_t Y_start,
        ITYPE_t Y_end,
        ITYPE_t thread_num,
    ) nogil
