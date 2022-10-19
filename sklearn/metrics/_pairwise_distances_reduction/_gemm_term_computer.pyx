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

import numpy as np
from scipy.sparse import issparse, csr_matrix
from ...utils._typedefs import DTYPE, SPARSE_INDEX_TYPE

# TODO: change for `libcpp.algorithm.fill` once Cython 3 is used
# Introduction in Cython:
# https://github.com/cython/cython/blob/05059e2a9b89bf6738a7750b905057e5b1e3fe2e/Cython/Includes/libcpp/algorithm.pxd#L50 #noqa
cdef extern from "<algorithm>" namespace "std" nogil:
    void fill[Iter, T](Iter first, Iter last, const T& value) except + #noqa



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


cdef class MiddleTermComputer64:
    """Component for `EuclideanDistance` specialisation wrapping the logic for the call to GEMM.

    `EuclideanDistance` subclasses internally compute the squared Euclidean distances
    between chunks of vectors X_c and Y_c using the following decomposition:


                ||X_c_i - Y_c_j||² = ||X_c_i||² - 2 X_c_i.Y_c_j^T + ||Y_c_j||²


    This helper class is in charge of wrapping the common logic to compute
    the middle term `- 2 X_c_i.Y_c_j^T` with a call to GEMM, which has a high
    arithmetic intensity.
    """

    @classmethod
    def get_for(
        cls,
        X,
        Y,
        effective_n_threads,
        chunks_n_threads,
        dist_middle_terms_chunks_size,
        n_features,
        chunk_size,
    ) -> MiddleTermComputer64:
        """Return the DatasetsPair implementation for the given arguments.

        Parameters
        ----------
        X : {ndarray, sparse matrix} of shape (n_samples_X, n_features)
            Input data.
            If provided as a ndarray, it must be C-contiguous.
            If provided as a sparse matrix, it must be in CSR format.

        Y : {ndarray, sparse matrix} of shape (n_samples_Y, n_features)
            Input data.
            If provided as a ndarray, it must be C-contiguous.
            If provided as a sparse matrix, it must be in CSR format.

        Returns
        -------
        middle_term_computer: MiddleTermComputer64
            The suited MiddleTermComputer64 implementation.
        """
        X_is_sparse = issparse(X)
        Y_is_sparse = issparse(Y)

        if not X_is_sparse and not Y_is_sparse:
            return DenseDenseMiddleTermComputer64(
                X,
                Y,
                effective_n_threads,
                chunks_n_threads,
                dist_middle_terms_chunks_size,
                n_features,
                chunk_size,
            )
        if X_is_sparse and Y_is_sparse:
            return SparseSparseMiddleTermComputer64(
                X,
                Y,
                effective_n_threads,
                chunks_n_threads,
                dist_middle_terms_chunks_size,
                n_features,
                chunk_size,
            )

        raise NotImplementedError("X and Y must be both sparse or dense")


    @classmethod
    def unpack_csr_matrix(cls, X: csr_matrix):
        """Ensure that the CSR matrix is indexed with SPARSE_INDEX_TYPE."""
        X_data = np.asarray(X.data, dtype=DTYPE)
        X_indices = np.asarray(X.indices, dtype=SPARSE_INDEX_TYPE)
        X_indptr = np.asarray(X.indptr, dtype=SPARSE_INDEX_TYPE)
        return X_data, X_indices, X_indptr

    def __init__(
        self,
        ITYPE_t effective_n_threads,
        ITYPE_t chunks_n_threads,
        ITYPE_t dist_middle_terms_chunks_size,
        ITYPE_t n_features,
        ITYPE_t chunk_size,
    ):
        self.effective_n_threads = effective_n_threads
        self.chunks_n_threads = chunks_n_threads
        self.dist_middle_terms_chunks_size = dist_middle_terms_chunks_size
        self.n_features = n_features
        self.chunk_size = chunk_size

        self.dist_middle_terms_chunks = vector[vector[DTYPE_t]](self.effective_n_threads)

    cdef void _parallel_on_X_pre_compute_and_reduce_distances_on_chunks(
        self,
        ITYPE_t X_start,
        ITYPE_t X_end,
        ITYPE_t Y_start,
        ITYPE_t Y_end,
        ITYPE_t thread_num,
    ) nogil:
        return

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
        return NULL


cdef class DenseDenseMiddleTermComputer64(MiddleTermComputer64):

    def __init__(
        self,
        const DTYPE_t[:, ::1] X,
        const DTYPE_t[:, ::1] Y,
        ITYPE_t effective_n_threads,
        ITYPE_t chunks_n_threads,
        ITYPE_t dist_middle_terms_chunks_size,
        ITYPE_t n_features,
        ITYPE_t chunk_size,
    ):
        super().__init__(
            effective_n_threads,
            chunks_n_threads,
            dist_middle_terms_chunks_size,
            n_features,
            chunk_size,
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
            const DTYPE_t[:, ::1] X_c = self.X[X_start:X_end, :]
            const DTYPE_t[:, ::1] Y_c = self.Y[Y_start:Y_end, :]
            DTYPE_t *dist_middle_terms = self.dist_middle_terms_chunks[thread_num].data()

            # Careful: LDA, LDB and LDC are given for F-ordered arrays
            # in BLAS documentations, for instance:
            # https://www.netlib.org/lapack/explore-html/db/dc9/group__single__blas__level3_gafe51bacb54592ff5de056acabd83c260.html #noqa
            #
            # Here, we use their counterpart values to work with C-ordered arrays.
            BLAS_Order order = RowMajor
            BLAS_Trans ta = NoTrans
            BLAS_Trans tb = Trans
            ITYPE_t m = X_c.shape[0]
            ITYPE_t n = Y_c.shape[0]
            ITYPE_t K = X_c.shape[1]
            DTYPE_t alpha = - 2.
            # Casting for A and B to remove the const is needed because APIs exposed via
            # scipy.linalg.cython_blas aren't reflecting the arguments' const qualifier.
            # See: https://github.com/scipy/scipy/issues/14262
            DTYPE_t * A = <DTYPE_t *> &X_c[0, 0]
            DTYPE_t * B = <DTYPE_t *> &Y_c[0, 0]
            ITYPE_t lda = X_c.shape[1]
            ITYPE_t ldb = X_c.shape[1]
            DTYPE_t beta = 0.
            ITYPE_t ldc = Y_c.shape[0]

        # dist_middle_terms = `-2 * X_c @ Y_c.T`
        _gemm(order, ta, tb, m, n, K, alpha, A, lda, B, ldb, beta, dist_middle_terms, ldc)

        return dist_middle_terms


cdef class SparseSparseMiddleTermComputer64(MiddleTermComputer64):

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
            effective_n_threads,
            chunks_n_threads,
            dist_middle_terms_chunks_size,
            n_features,
            chunk_size,
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
) nogil:
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


cdef class MiddleTermComputer32:
    """Component for `EuclideanDistance` specialisation wrapping the logic for the call to GEMM.

    `EuclideanDistance` subclasses internally compute the squared Euclidean distances
    between chunks of vectors X_c and Y_c using the following decomposition:


                ||X_c_i - Y_c_j||² = ||X_c_i||² - 2 X_c_i.Y_c_j^T + ||Y_c_j||²


    This helper class is in charge of wrapping the common logic to compute
    the middle term `- 2 X_c_i.Y_c_j^T` with a call to GEMM, which has a high
    arithmetic intensity.
    """

    @classmethod
    def get_for(
        cls,
        X,
        Y,
        effective_n_threads,
        chunks_n_threads,
        dist_middle_terms_chunks_size,
        n_features,
        chunk_size,
    ) -> MiddleTermComputer32:
        """Return the DatasetsPair implementation for the given arguments.

        Parameters
        ----------
        X : {ndarray, sparse matrix} of shape (n_samples_X, n_features)
            Input data.
            If provided as a ndarray, it must be C-contiguous.
            If provided as a sparse matrix, it must be in CSR format.

        Y : {ndarray, sparse matrix} of shape (n_samples_Y, n_features)
            Input data.
            If provided as a ndarray, it must be C-contiguous.
            If provided as a sparse matrix, it must be in CSR format.

        Returns
        -------
        middle_term_computer: MiddleTermComputer32
            The suited MiddleTermComputer32 implementation.
        """
        X_is_sparse = issparse(X)
        Y_is_sparse = issparse(Y)

        if not X_is_sparse and not Y_is_sparse:
            return DenseDenseMiddleTermComputer32(
                X,
                Y,
                effective_n_threads,
                chunks_n_threads,
                dist_middle_terms_chunks_size,
                n_features,
                chunk_size,
            )
        if X_is_sparse and Y_is_sparse:
            return SparseSparseMiddleTermComputer32(
                X,
                Y,
                effective_n_threads,
                chunks_n_threads,
                dist_middle_terms_chunks_size,
                n_features,
                chunk_size,
            )

        raise NotImplementedError("X and Y must be both sparse or dense")


    @classmethod
    def unpack_csr_matrix(cls, X: csr_matrix):
        """Ensure that the CSR matrix is indexed with SPARSE_INDEX_TYPE."""
        X_data = np.asarray(X.data, dtype=np.float32)
        X_indices = np.asarray(X.indices, dtype=SPARSE_INDEX_TYPE)
        X_indptr = np.asarray(X.indptr, dtype=SPARSE_INDEX_TYPE)
        return X_data, X_indices, X_indptr

    def __init__(
        self,
        ITYPE_t effective_n_threads,
        ITYPE_t chunks_n_threads,
        ITYPE_t dist_middle_terms_chunks_size,
        ITYPE_t n_features,
        ITYPE_t chunk_size,
    ):
        self.effective_n_threads = effective_n_threads
        self.chunks_n_threads = chunks_n_threads
        self.dist_middle_terms_chunks_size = dist_middle_terms_chunks_size
        self.n_features = n_features
        self.chunk_size = chunk_size

        self.dist_middle_terms_chunks = vector[vector[DTYPE_t]](self.effective_n_threads)

    cdef void _parallel_on_X_pre_compute_and_reduce_distances_on_chunks(
        self,
        ITYPE_t X_start,
        ITYPE_t X_end,
        ITYPE_t Y_start,
        ITYPE_t Y_end,
        ITYPE_t thread_num,
    ) nogil:
        return

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
        return NULL


cdef class DenseDenseMiddleTermComputer32(MiddleTermComputer32):

    def __init__(
        self,
        const cnp.float32_t[:, ::1] X,
        const cnp.float32_t[:, ::1] Y,
        ITYPE_t effective_n_threads,
        ITYPE_t chunks_n_threads,
        ITYPE_t dist_middle_terms_chunks_size,
        ITYPE_t n_features,
        ITYPE_t chunk_size,
    ):
        super().__init__(
            effective_n_threads,
            chunks_n_threads,
            dist_middle_terms_chunks_size,
            n_features,
            chunk_size,
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
            const cnp.float32_t[:, ::1] X_c = self.X[X_start:X_end, :]
            const cnp.float32_t[:, ::1] Y_c = self.Y[Y_start:Y_end, :]
            DTYPE_t *dist_middle_terms = self.dist_middle_terms_chunks[thread_num].data()

            # Careful: LDA, LDB and LDC are given for F-ordered arrays
            # in BLAS documentations, for instance:
            # https://www.netlib.org/lapack/explore-html/db/dc9/group__single__blas__level3_gafe51bacb54592ff5de056acabd83c260.html #noqa
            #
            # Here, we use their counterpart values to work with C-ordered arrays.
            BLAS_Order order = RowMajor
            BLAS_Trans ta = NoTrans
            BLAS_Trans tb = Trans
            ITYPE_t m = X_c.shape[0]
            ITYPE_t n = Y_c.shape[0]
            ITYPE_t K = X_c.shape[1]
            DTYPE_t alpha = - 2.
            DTYPE_t * A = self.X_c_upcast[thread_num].data()
            DTYPE_t * B = self.Y_c_upcast[thread_num].data()
            ITYPE_t lda = X_c.shape[1]
            ITYPE_t ldb = X_c.shape[1]
            DTYPE_t beta = 0.
            ITYPE_t ldc = Y_c.shape[0]

        # dist_middle_terms = `-2 * X_c @ Y_c.T`
        _gemm(order, ta, tb, m, n, K, alpha, A, lda, B, ldb, beta, dist_middle_terms, ldc)

        return dist_middle_terms


cdef class SparseSparseMiddleTermComputer32(MiddleTermComputer32):

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
            effective_n_threads,
            chunks_n_threads,
            dist_middle_terms_chunks_size,
            n_features,
            chunk_size,
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

        _middle_term_sparse_sparse_32(
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
