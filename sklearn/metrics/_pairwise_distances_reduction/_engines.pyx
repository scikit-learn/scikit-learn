cimport numpy as cnp

from libcpp.vector cimport vector

from ...utils._cython_blas cimport (
  BLAS_Order,
  BLAS_Trans,
  NoTrans,
  RowMajor,
  Trans,
  _gemm,
)
from ...utils._typedefs cimport DTYPE_t, ITYPE_t, SPARSE_INDEX_TYPE_t

# TODO: change for `libcpp.algorithm.fill` once Cython 3 is used
# Introduction in Cython:
#
# https://github.com/cython/cython/blob/05059e2a9b89bf6738a7750b905057e5b1e3fe2e/Cython/Includes/libcpp/algorithm.pxd#L50 #noqa
cdef extern from "<algorithm>" namespace "std" nogil:
    void fill[Iter, T](Iter first, Iter last, const T& value) except + #noqa

import numpy as np
from scipy.sparse import issparse, csr_matrix
from ...utils._typedefs import DTYPE, SPARSE_INDEX_TYPE
from ...utils import check_array

cdef class BaseEngine:
    def __init__(self):
        return

    cdef void _parallel_on_X_parallel_init(
        self,
        ITYPE_t thread_num,
    ) nogil:
        return

    cdef void _parallel_on_X_init_chunk(
        self,
        ITYPE_t thread_num,
        ITYPE_t X_start,
        ITYPE_t X_end,
    ) nogil:
        return

    cdef void _parallel_on_X_pre_compute_and_reduce_distances_on_chunks(
        self,
        ITYPE_t X_start,
        ITYPE_t X_end,
        ITYPE_t Y_start,
        ITYPE_t Y_end,
        ITYPE_t thread_num,
    ) nogil:
        return

    cdef void _parallel_on_X_prange_iter_finalize(
        self,
        ITYPE_t thread_num,
        ITYPE_t X_start,
        ITYPE_t X_end,
    ) nogil:
        return

    cdef void _parallel_on_X_parallel_finalize(
        self,
        ITYPE_t thread_num
    ) nogil:
        return

    cdef void _parallel_on_Y_init(
        self,
    ) nogil:
        return

    cdef void _parallel_on_Y_parallel_init(
        self,
        ITYPE_t thread_num,
        ITYPE_t X_start,
        ITYPE_t X_end,
    ) nogil:
        return

    cdef void _parallel_on_Y_pre_compute_and_reduce_distances_on_chunks(
        self,
        ITYPE_t X_start,
        ITYPE_t X_end,
        ITYPE_t Y_start,
        ITYPE_t Y_end,
        ITYPE_t thread_num,
    ) nogil:
        return

    cdef void _parallel_on_Y_synchronize(
        self,
        ITYPE_t X_start,
        ITYPE_t X_end,
    ) nogil:
        return

    cdef void _compute_and_reduce_distances_on_chunks(
        self,
        ITYPE_t X_start,
        ITYPE_t X_end,
        ITYPE_t Y_start,
        ITYPE_t Y_end,
        ITYPE_t thread_num,
    ) nogil:
        return

# TODO: If possible optimize this routine to efficiently treat cases where
# `n_samples_X << n_samples_Y` met in practise when X_test consists of a
# few samples, and thus when there's a single chunk of X whose number of
# samples is less that the default chunk size.

# TODO: compare this routine with the similar ones in SciPy, especially
# `csr_matmat` which might implement a better algorithm.
# See: https://github.com/scipy/scipy/blob/e58292e066ba2cb2f3d1e0563ca9314ff1f4f311/scipy/sparse/sparsetools/csr.h#L603-L669  # noqa
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
) nogil:
    # This routine assumes that D points to the first element of a
    # zeroed buffer of length at least equal to n_X × n_Y, conceptually
    # representing a 2-d C-ordered array.
    cdef:
        ITYPE_t i, j, k
        ITYPE_t n_X = X_end - X_start
        ITYPE_t n_Y = Y_end - Y_start
        ITYPE_t X_i_col_idx, X_i_ptr, Y_j_col_idx, Y_j_ptr

    for i in range(n_X):
        for X_i_ptr in range(X_indptr[X_start+i], X_indptr[X_start+i+1]):
            X_i_col_idx = X_indices[X_i_ptr]
            for j in range(n_Y):
                k = i * n_Y + j
                for Y_j_ptr in range(Y_indptr[Y_start+j], Y_indptr[Y_start+j+1]):
                    Y_j_col_idx = Y_indices[Y_j_ptr]
                    if X_i_col_idx == Y_j_col_idx:
                        D[k] += -2 * X_data[X_i_ptr] * Y_data[Y_j_ptr]


from ._base cimport _sqeuclidean_row_norms64

cdef class EuclideanEngine64(BaseEngine):
    """Helper class to compute a Euclidean distance matrix in chunks.

    This is an abstract base class that is further specialized depending
    on the type of data (dense or sparse).

    `EuclideanDistance` subclasses relies on the squared Euclidean
    distances between chunks of vectors X_c and Y_c using the
    following decomposition for the (i,j) pair :


         ||X_c_i - Y_c_j||² = ||X_c_i||² - 2 X_c_i.Y_c_j^T + ||Y_c_j||²


    This helper class is in charge of wrapping the common logic to compute
    the middle term, i.e. `- 2 X_c_i.Y_c_j^T`.
    """

    @classmethod
    def get_for(
        cls,
        X,
        Y,
        pda,
    ) -> EuclideanEngine64:
        """Return the DatasetsPair implementation for the given arguments.

        Parameters
        ----------
        X : ndarray or CSR sparse matrix of shape (n_samples_X, n_features)
            Input data.
            If provided as a ndarray, it must be C-contiguous.

        Y : ndarray or CSR sparse matrix of shape (n_samples_Y, n_features)
            Input data.
            If provided as a ndarray, it must be C-contiguous.

        Returns
        -------
        engine: EuclideanEngine64
            The suited EuclideanEngine64 implementation.
        """
        X_is_sparse = issparse(X)
        Y_is_sparse = issparse(Y)
        dist_middle_terms_chunks_size = pda.Y_n_samples_chunk * pda.X_n_samples_chunk
        if not X_is_sparse and not Y_is_sparse:
            return DenseDenseEuclideanEngine64(
                X,
                Y,
                effective_n_threads=pda.effective_n_threads,
                chunks_n_threads=pda.chunks_n_threads,
                dist_middle_terms_chunks_size=dist_middle_terms_chunks_size,
                chunk_size=pda.chunk_size,
                metric_kwargs=pda.metric_kwargs,
            )
        if X_is_sparse and Y_is_sparse:
            return SparseSparseEuclideanEngine64(
                X,
                Y,
                effective_n_threads=pda.effective_n_threads,
                chunks_n_threads=pda.chunks_n_threads,
                dist_middle_terms_chunks_size=dist_middle_terms_chunks_size,
                chunk_size=pda.chunk_size,
                metric_kwargs=pda.metric_kwargs,
            )

        raise NotImplementedError(
            "X and Y must be both CSR sparse matrices or both numpy arrays."
        )


    @classmethod
    def unpack_csr_matrix(cls, X: csr_matrix):
        """Ensure that the CSR matrix is indexed with SPARSE_INDEX_TYPE."""
        X_data = np.asarray(X.data, dtype=DTYPE)
        X_indices = np.asarray(X.indices, dtype=SPARSE_INDEX_TYPE)
        X_indptr = np.asarray(X.indptr, dtype=SPARSE_INDEX_TYPE)
        return X_data, X_indices, X_indptr

    def __init__(
        self,
        X,
        Y,
        ITYPE_t effective_n_threads,
        ITYPE_t chunks_n_threads,
        ITYPE_t dist_middle_terms_chunks_size,
        ITYPE_t chunk_size,
        dict metric_kwargs=None,
    ):
        self.effective_n_threads = effective_n_threads
        self.chunks_n_threads = chunks_n_threads
        self.dist_middle_terms_chunks_size = dist_middle_terms_chunks_size
        self.n_features = X.shape[1]
        self.chunk_size = chunk_size

        self.dist_middle_terms_chunks = vector[vector[DTYPE_t]](self.effective_n_threads)

        if metric_kwargs is not None and "Y_norm_squared" in metric_kwargs:
            self.Y_norm_squared = check_array(
                metric_kwargs.pop("Y_norm_squared"),
                ensure_2d=False,
                input_name="Y_norm_squared",
                dtype=np.float64,
            )
        else:
            self.Y_norm_squared = _sqeuclidean_row_norms64(
                Y,
                self.effective_n_threads,
            )

        if metric_kwargs is not None and "X_norm_squared" in metric_kwargs:
            self.X_norm_squared = check_array(
                metric_kwargs.pop("X_norm_squared"),
                ensure_2d=False,
                input_name="X_norm_squared",
                dtype=np.float64,
            )
        else:
            # Do not recompute norms if datasets are identical.
            self.X_norm_squared = (
                self.Y_norm_squared if X is Y else
                _sqeuclidean_row_norms64(
                    X,
                    self.effective_n_threads,
                )
            )

    cdef void _parallel_on_X_parallel_init(self, ITYPE_t thread_num) nogil:
        self.dist_middle_terms_chunks[thread_num].resize(self.dist_middle_terms_chunks_size)

    cdef void _parallel_on_X_init_chunk(
        self,
        ITYPE_t thread_num,
        ITYPE_t X_start,
        ITYPE_t X_end,
    ) nogil:
        return

    cdef void _parallel_on_Y_init(self) nogil:
        for thread_num in range(self.chunks_n_threads):
            self.dist_middle_terms_chunks[thread_num].resize(
                self.dist_middle_terms_chunks_size
            )

    cdef DTYPE_t * _compute_dist_middle_terms(
        self,
        ITYPE_t X_start,
        ITYPE_t X_end,
        ITYPE_t Y_start,
        ITYPE_t Y_end,
        ITYPE_t thread_num,
    ) nogil:
        return NULL

    cdef DTYPE_t _compute_pair_distance(
        self,
        ITYPE_t i, # Index of X sample
        ITYPE_t j, # Index of Y sample
        ITYPE_t X_start, # Index offset
        ITYPE_t Y_start, # Index offset
        DTYPE_t * dist_middle_terms, # Array of pre-computeted middle terms
    ) nogil:

        cdef ITYPE_t n_Y = len(self.Y_norm_squared)
        # Index of middle term
        cdef ITYPE_t k = n_Y * i + j
        cdef DTYPE_t val = (
            self.X_norm_squared[i + X_start] +
            dist_middle_terms[i * n_Y + j] +
            self.Y_norm_squared[j + Y_start]
        )
        # Catastrophic cancellation might cause -0. to be present,
        # e.g. when computing d(x_i, y_i) when X is Y.
        return max(0., val)


cdef class DenseDenseEuclideanEngine64(EuclideanEngine64):
    """Computes the middle term of the Euclidean distance between two chunked dense matrices
    X_c and Y_c.

                        dist_middle_terms = - 2 X_c_i.Y_c_j^T

    This class use the BLAS gemm routine to perform the dot product of each chunks
    of the distance matrix with improved arithmetic intensity and vector instruction (SIMD).
    """

    def __init__(
        self,
        const DTYPE_t[:, ::1] X,
        const DTYPE_t[:, ::1] Y,
        ITYPE_t effective_n_threads,
        ITYPE_t chunks_n_threads,
        ITYPE_t dist_middle_terms_chunks_size,
        ITYPE_t n_features,
        ITYPE_t chunk_size,
        dict metric_kwargs=None,
    ):
        super().__init__(
            X, Y,
            effective_n_threads,
            chunks_n_threads,
            dist_middle_terms_chunks_size,
            n_features,
            chunk_size,
            metric_kwargs=None,
        )
        self.X = X
        self.Y = Y

    cdef void _parallel_on_X_pre_compute_and_reduce_distances_on_chunks(
        self,
        ITYPE_t X_start,
        ITYPE_t X_end,
        ITYPE_t Y_start,
        ITYPE_t Y_end,
        ITYPE_t thread_num,
    ) nogil:
        return

    cdef void _parallel_on_X_init_chunk(
        self,
        ITYPE_t thread_num,
        ITYPE_t X_start,
        ITYPE_t X_end,
    ) nogil:
        return

    cdef void _parallel_on_Y_parallel_init(
        self,
        ITYPE_t thread_num,
        ITYPE_t X_start,
        ITYPE_t X_end,
    ) nogil:
        return

    cdef void _parallel_on_Y_pre_compute_and_reduce_distances_on_chunks(
        self,
        ITYPE_t X_start,
        ITYPE_t X_end,
        ITYPE_t Y_start,
        ITYPE_t Y_end,
        ITYPE_t thread_num
    ) nogil:
        return

    cdef DTYPE_t * _compute_dist_middle_terms(
        self,
        ITYPE_t X_start,
        ITYPE_t X_end,
        ITYPE_t Y_start,
        ITYPE_t Y_end,
        ITYPE_t thread_num,
    ) nogil:
        cdef:
            DTYPE_t *dist_middle_terms = self.dist_middle_terms_chunks[thread_num].data()

            # Careful: LDA, LDB and LDC are given for F-ordered arrays
            # in BLAS documentations, for instance:
            # https://www.netlib.org/lapack/explore-html/db/dc9/group__single__blas__level3_gafe51bacb54592ff5de056acabd83c260.html #noqa
            #
            # Here, we use their counterpart values to work with C-ordered arrays.
            BLAS_Order order = RowMajor
            BLAS_Trans ta = NoTrans
            BLAS_Trans tb = Trans
            ITYPE_t m = X_end - X_start
            ITYPE_t n = Y_end - Y_start
            ITYPE_t K = self.n_features
            DTYPE_t alpha = - 2.
            # Casting for A and B to remove the const is needed because APIs exposed via
            # scipy.linalg.cython_blas aren't reflecting the arguments' const qualifier.
            # See: https://github.com/scipy/scipy/issues/14262
            DTYPE_t * A = <DTYPE_t *> &self.X[X_start, 0]
            DTYPE_t * B = <DTYPE_t *> &self.Y[Y_start, 0]
            ITYPE_t lda = self.n_features
            ITYPE_t ldb = self.n_features
            DTYPE_t beta = 0.
            ITYPE_t ldc = Y_end - Y_start

        # dist_middle_terms = `-2 * X[X_start:X_end] @ Y[Y_start:Y_end].T`
        _gemm(order, ta, tb, m, n, K, alpha, A, lda, B, ldb, beta, dist_middle_terms, ldc)

        return dist_middle_terms


cdef class SparseSparseEuclideanEngine64(EuclideanEngine64):
    """Middle term of the Euclidean distance between two chunked CSR matrices.

    The result is return as a contiguous array.

            dist_middle_terms = - 2 X_c_i.Y_c_j^T

    The logic of the computation is wrapped in the routine _middle_term_sparse_sparse_64.
    This routine iterates over the data, indices and indptr arrays of the sparse matrices without
    densifying them.
    """

    def __init__(
        self,
        X,
        Y,
        ITYPE_t effective_n_threads,
        ITYPE_t chunks_n_threads,
        ITYPE_t dist_middle_terms_chunks_size,
        ITYPE_t n_features,
        ITYPE_t chunk_size,
    ):
        super().__init__(
            X, Y,
            effective_n_threads,
            chunks_n_threads,
            dist_middle_terms_chunks_size,
            n_features,
            chunk_size,
            metric_kwargs=None,
        )
        self.X_data, self.X_indices, self.X_indptr = self.unpack_csr_matrix(X)
        self.Y_data, self.Y_indices, self.Y_indptr = self.unpack_csr_matrix(Y)

    cdef void _parallel_on_X_pre_compute_and_reduce_distances_on_chunks(
        self,
        ITYPE_t X_start,
        ITYPE_t X_end,
        ITYPE_t Y_start,
        ITYPE_t Y_end,
        ITYPE_t thread_num,
    ) nogil:
        # Flush the thread dist_middle_terms_chunks to 0.0
        fill(
            self.dist_middle_terms_chunks[thread_num].begin(),
            self.dist_middle_terms_chunks[thread_num].end(),
            0.0,
        )

    cdef void _parallel_on_Y_pre_compute_and_reduce_distances_on_chunks(
        self,
        ITYPE_t X_start,
        ITYPE_t X_end,
        ITYPE_t Y_start,
        ITYPE_t Y_end,
        ITYPE_t thread_num,
    ) nogil:
        # Flush the thread dist_middle_terms_chunks to 0.0
        fill(
            self.dist_middle_terms_chunks[thread_num].begin(),
            self.dist_middle_terms_chunks[thread_num].end(),
            0.0,
        )

    cdef DTYPE_t * _compute_dist_middle_terms(
        self,
        ITYPE_t X_start,
        ITYPE_t X_end,
        ITYPE_t Y_start,
        ITYPE_t Y_end,
        ITYPE_t thread_num,
    ) nogil:
        cdef:
            DTYPE_t *dist_middle_terms = (
                self.dist_middle_terms_chunks[thread_num].data()
            )

        _middle_term_sparse_sparse_64(
            self.X_data,
            self.X_indices,
            self.X_indptr,
            X_start,
            X_end,
            self.Y_data,
            self.Y_indices,
            self.Y_indptr,
            Y_start,
            Y_end,
            dist_middle_terms,
        )

        return dist_middle_terms

from ._base cimport _sqeuclidean_row_norms32

cdef class EuclideanEngine32(BaseEngine):
    """Helper class to compute a Euclidean distance matrix in chunks.

    This is an abstract base class that is further specialized depending
    on the type of data (dense or sparse).

    `EuclideanDistance` subclasses relies on the squared Euclidean
    distances between chunks of vectors X_c and Y_c using the
    following decomposition for the (i,j) pair :


         ||X_c_i - Y_c_j||² = ||X_c_i||² - 2 X_c_i.Y_c_j^T + ||Y_c_j||²


    This helper class is in charge of wrapping the common logic to compute
    the middle term, i.e. `- 2 X_c_i.Y_c_j^T`.
    """

    @classmethod
    def get_for(
        cls,
        X,
        Y,
        pda,
    ) -> EuclideanEngine32:
        """Return the DatasetsPair implementation for the given arguments.

        Parameters
        ----------
        X : ndarray or CSR sparse matrix of shape (n_samples_X, n_features)
            Input data.
            If provided as a ndarray, it must be C-contiguous.

        Y : ndarray or CSR sparse matrix of shape (n_samples_Y, n_features)
            Input data.
            If provided as a ndarray, it must be C-contiguous.

        Returns
        -------
        engine: EuclideanEngine32
            The suited EuclideanEngine32 implementation.
        """
        X_is_sparse = issparse(X)
        Y_is_sparse = issparse(Y)
        dist_middle_terms_chunks_size = pda.Y_n_samples_chunk * pda.X_n_samples_chunk
        if not X_is_sparse and not Y_is_sparse:
            return DenseDenseEuclideanEngine32(
                X,
                Y,
                effective_n_threads=pda.effective_n_threads,
                chunks_n_threads=pda.chunks_n_threads,
                dist_middle_terms_chunks_size=dist_middle_terms_chunks_size,
                chunk_size=pda.chunk_size,
                metric_kwargs=pda.metric_kwargs,
            )
        if X_is_sparse and Y_is_sparse:
            return SparseSparseEuclideanEngine32(
                X,
                Y,
                effective_n_threads=pda.effective_n_threads,
                chunks_n_threads=pda.chunks_n_threads,
                dist_middle_terms_chunks_size=dist_middle_terms_chunks_size,
                chunk_size=pda.chunk_size,
                metric_kwargs=pda.metric_kwargs,
            )

        raise NotImplementedError(
            "X and Y must be both CSR sparse matrices or both numpy arrays."
        )


    @classmethod
    def unpack_csr_matrix(cls, X: csr_matrix):
        """Ensure that the CSR matrix is indexed with SPARSE_INDEX_TYPE."""
        X_data = np.asarray(X.data, dtype=DTYPE)
        X_indices = np.asarray(X.indices, dtype=SPARSE_INDEX_TYPE)
        X_indptr = np.asarray(X.indptr, dtype=SPARSE_INDEX_TYPE)
        return X_data, X_indices, X_indptr

    def __init__(
        self,
        X,
        Y,
        ITYPE_t effective_n_threads,
        ITYPE_t chunks_n_threads,
        ITYPE_t dist_middle_terms_chunks_size,
        ITYPE_t chunk_size,
        dict metric_kwargs=None,
    ):
        self.effective_n_threads = effective_n_threads
        self.chunks_n_threads = chunks_n_threads
        self.dist_middle_terms_chunks_size = dist_middle_terms_chunks_size
        self.n_features = X.shape[1]
        self.chunk_size = chunk_size

        self.dist_middle_terms_chunks = vector[vector[DTYPE_t]](self.effective_n_threads)

        if metric_kwargs is not None and "Y_norm_squared" in metric_kwargs:
            self.Y_norm_squared = check_array(
                metric_kwargs.pop("Y_norm_squared"),
                ensure_2d=False,
                input_name="Y_norm_squared",
                dtype=np.float64,
            )
        else:
            self.Y_norm_squared = _sqeuclidean_row_norms32(
                Y,
                self.effective_n_threads,
            )

        if metric_kwargs is not None and "X_norm_squared" in metric_kwargs:
            self.X_norm_squared = check_array(
                metric_kwargs.pop("X_norm_squared"),
                ensure_2d=False,
                input_name="X_norm_squared",
                dtype=np.float64,
            )
        else:
            # Do not recompute norms if datasets are identical.
            self.X_norm_squared = (
                self.Y_norm_squared if X is Y else
                _sqeuclidean_row_norms32(
                    X,
                    self.effective_n_threads,
                )
            )

    cdef void _parallel_on_X_parallel_init(self, ITYPE_t thread_num) nogil:
        self.dist_middle_terms_chunks[thread_num].resize(self.dist_middle_terms_chunks_size)

    cdef void _parallel_on_X_init_chunk(
        self,
        ITYPE_t thread_num,
        ITYPE_t X_start,
        ITYPE_t X_end,
    ) nogil:
        return

    cdef void _parallel_on_Y_init(self) nogil:
        for thread_num in range(self.chunks_n_threads):
            self.dist_middle_terms_chunks[thread_num].resize(
                self.dist_middle_terms_chunks_size
            )

    cdef DTYPE_t * _compute_dist_middle_terms(
        self,
        ITYPE_t X_start,
        ITYPE_t X_end,
        ITYPE_t Y_start,
        ITYPE_t Y_end,
        ITYPE_t thread_num,
    ) nogil:
        return NULL

    cdef DTYPE_t _compute_pair_distance(
        self,
        ITYPE_t i, # Index of X sample
        ITYPE_t j, # Index of Y sample
        ITYPE_t X_start, # Index offset
        ITYPE_t Y_start, # Index offset
        DTYPE_t * dist_middle_terms, # Array of pre-computeted middle terms
    ) nogil:

        cdef ITYPE_t n_Y = len(self.Y_norm_squared)
        # Index of middle term
        cdef ITYPE_t k = n_Y * i + j
        cdef DTYPE_t val = (
            self.X_norm_squared[i + X_start] +
            dist_middle_terms[i * n_Y + j] +
            self.Y_norm_squared[j + Y_start]
        )
        # Catastrophic cancellation might cause -0. to be present,
        # e.g. when computing d(x_i, y_i) when X is Y.
        return max(0., val)


cdef class DenseDenseEuclideanEngine32(EuclideanEngine32):
    """Computes the middle term of the Euclidean distance between two chunked dense matrices
    X_c and Y_c.

                        dist_middle_terms = - 2 X_c_i.Y_c_j^T

    This class use the BLAS gemm routine to perform the dot product of each chunks
    of the distance matrix with improved arithmetic intensity and vector instruction (SIMD).
    """

    def __init__(
        self,
        const cnp.float32_t[:, ::1] X,
        const cnp.float32_t[:, ::1] Y,
        ITYPE_t effective_n_threads,
        ITYPE_t chunks_n_threads,
        ITYPE_t dist_middle_terms_chunks_size,
        ITYPE_t n_features,
        ITYPE_t chunk_size,
        dict metric_kwargs=None,
    ):
        super().__init__(
            X, Y,
            effective_n_threads,
            chunks_n_threads,
            dist_middle_terms_chunks_size,
            n_features,
            chunk_size,
            metric_kwargs=None,
        )
        self.X = X
        self.Y = Y
        # We populate the buffer for upcasting chunks of X and Y from float32 to float64.
        self.X_c_upcast = vector[vector[DTYPE_t]](self.effective_n_threads)
        self.Y_c_upcast = vector[vector[DTYPE_t]](self.effective_n_threads)

        upcast_buffer_n_elements = self.chunk_size * n_features

        for thread_num in range(self.effective_n_threads):
            self.X_c_upcast[thread_num].resize(upcast_buffer_n_elements)
            self.Y_c_upcast[thread_num].resize(upcast_buffer_n_elements)

    cdef void _parallel_on_X_pre_compute_and_reduce_distances_on_chunks(
        self,
        ITYPE_t X_start,
        ITYPE_t X_end,
        ITYPE_t Y_start,
        ITYPE_t Y_end,
        ITYPE_t thread_num,
    ) nogil:
        cdef:
            ITYPE_t i, j
            ITYPE_t n_chunk_samples = Y_end - Y_start

        # Upcasting Y_c=Y[Y_start:Y_end, :] from float32 to float64
        for i in range(n_chunk_samples):
            for j in range(self.n_features):
                self.Y_c_upcast[thread_num][i * self.n_features + j] = <DTYPE_t> self.Y[Y_start + i, j]

    cdef void _parallel_on_X_init_chunk(
        self,
        ITYPE_t thread_num,
        ITYPE_t X_start,
        ITYPE_t X_end,
    ) nogil:
        cdef:
            ITYPE_t i, j
            ITYPE_t n_chunk_samples = X_end - X_start

        # Upcasting X_c=X[X_start:X_end, :] from float32 to float64
        for i in range(n_chunk_samples):
            for j in range(self.n_features):
                self.X_c_upcast[thread_num][i * self.n_features + j] = <DTYPE_t> self.X[X_start + i, j]

    cdef void _parallel_on_Y_parallel_init(
        self,
        ITYPE_t thread_num,
        ITYPE_t X_start,
        ITYPE_t X_end,
    ) nogil:
        cdef:
            ITYPE_t i, j
            ITYPE_t n_chunk_samples = X_end - X_start

        # Upcasting X_c=X[X_start:X_end, :] from float32 to float64
        for i in range(n_chunk_samples):
            for j in range(self.n_features):
                self.X_c_upcast[thread_num][i * self.n_features + j] = <DTYPE_t> self.X[X_start + i, j]

    cdef void _parallel_on_Y_pre_compute_and_reduce_distances_on_chunks(
        self,
        ITYPE_t X_start,
        ITYPE_t X_end,
        ITYPE_t Y_start,
        ITYPE_t Y_end,
        ITYPE_t thread_num
    ) nogil:
        cdef:
            ITYPE_t i, j
            ITYPE_t n_chunk_samples = Y_end - Y_start

        # Upcasting Y_c=Y[Y_start:Y_end, :] from float32 to float64
        for i in range(n_chunk_samples):
            for j in range(self.n_features):
                self.Y_c_upcast[thread_num][i * self.n_features + j] = <DTYPE_t> self.Y[Y_start + i, j]

    cdef DTYPE_t * _compute_dist_middle_terms(
        self,
        ITYPE_t X_start,
        ITYPE_t X_end,
        ITYPE_t Y_start,
        ITYPE_t Y_end,
        ITYPE_t thread_num,
    ) nogil:
        cdef:
            DTYPE_t *dist_middle_terms = self.dist_middle_terms_chunks[thread_num].data()

            # Careful: LDA, LDB and LDC are given for F-ordered arrays
            # in BLAS documentations, for instance:
            # https://www.netlib.org/lapack/explore-html/db/dc9/group__single__blas__level3_gafe51bacb54592ff5de056acabd83c260.html #noqa
            #
            # Here, we use their counterpart values to work with C-ordered arrays.
            BLAS_Order order = RowMajor
            BLAS_Trans ta = NoTrans
            BLAS_Trans tb = Trans
            ITYPE_t m = X_end - X_start
            ITYPE_t n = Y_end - Y_start
            ITYPE_t K = self.n_features
            DTYPE_t alpha = - 2.
            DTYPE_t * A = self.X_c_upcast[thread_num].data()
            DTYPE_t * B = self.Y_c_upcast[thread_num].data()
            ITYPE_t lda = self.n_features
            ITYPE_t ldb = self.n_features
            DTYPE_t beta = 0.
            ITYPE_t ldc = Y_end - Y_start

        # dist_middle_terms = `-2 * X[X_start:X_end] @ Y[Y_start:Y_end].T`
        _gemm(order, ta, tb, m, n, K, alpha, A, lda, B, ldb, beta, dist_middle_terms, ldc)

        return dist_middle_terms


cdef class SparseSparseEuclideanEngine32(EuclideanEngine32):
    """Middle term of the Euclidean distance between two chunked CSR matrices.

    The result is return as a contiguous array.

            dist_middle_terms = - 2 X_c_i.Y_c_j^T

    The logic of the computation is wrapped in the routine _middle_term_sparse_sparse_64.
    This routine iterates over the data, indices and indptr arrays of the sparse matrices without
    densifying them.
    """

    def __init__(
        self,
        X,
        Y,
        ITYPE_t effective_n_threads,
        ITYPE_t chunks_n_threads,
        ITYPE_t dist_middle_terms_chunks_size,
        ITYPE_t n_features,
        ITYPE_t chunk_size,
    ):
        super().__init__(
            X, Y,
            effective_n_threads,
            chunks_n_threads,
            dist_middle_terms_chunks_size,
            n_features,
            chunk_size,
            metric_kwargs=None,
        )
        self.X_data, self.X_indices, self.X_indptr = self.unpack_csr_matrix(X)
        self.Y_data, self.Y_indices, self.Y_indptr = self.unpack_csr_matrix(Y)

    cdef void _parallel_on_X_pre_compute_and_reduce_distances_on_chunks(
        self,
        ITYPE_t X_start,
        ITYPE_t X_end,
        ITYPE_t Y_start,
        ITYPE_t Y_end,
        ITYPE_t thread_num,
    ) nogil:
        # Flush the thread dist_middle_terms_chunks to 0.0
        fill(
            self.dist_middle_terms_chunks[thread_num].begin(),
            self.dist_middle_terms_chunks[thread_num].end(),
            0.0,
        )

    cdef void _parallel_on_Y_pre_compute_and_reduce_distances_on_chunks(
        self,
        ITYPE_t X_start,
        ITYPE_t X_end,
        ITYPE_t Y_start,
        ITYPE_t Y_end,
        ITYPE_t thread_num,
    ) nogil:
        # Flush the thread dist_middle_terms_chunks to 0.0
        fill(
            self.dist_middle_terms_chunks[thread_num].begin(),
            self.dist_middle_terms_chunks[thread_num].end(),
            0.0,
        )

    cdef DTYPE_t * _compute_dist_middle_terms(
        self,
        ITYPE_t X_start,
        ITYPE_t X_end,
        ITYPE_t Y_start,
        ITYPE_t Y_end,
        ITYPE_t thread_num,
    ) nogil:
        cdef:
            DTYPE_t *dist_middle_terms = (
                self.dist_middle_terms_chunks[thread_num].data()
            )

        _middle_term_sparse_sparse_64(
            self.X_data,
            self.X_indices,
            self.X_indptr,
            X_start,
            X_end,
            self.Y_data,
            self.Y_indices,
            self.Y_indptr,
            Y_start,
            Y_end,
            dist_middle_terms,
        )

        return dist_middle_terms
