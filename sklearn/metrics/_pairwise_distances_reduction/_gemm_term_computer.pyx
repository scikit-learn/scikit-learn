from libcpp.vector cimport vector

from ...utils._typedefs cimport DTYPE_t, ITYPE_t

from ...utils._cython_blas cimport (
  BLAS_Order,
  BLAS_Trans,
  ColMajor,
  NoTrans,
  RowMajor,
  Trans,
  _gemm,
)

cdef class GEMMTermComputer64:
    """Component for `FastEuclidean*` variant wrapping the logic for the call to GEMM.

    `FastEuclidean*` classes internally compute the squared Euclidean distances between
    chunks of vectors X_c and Y_c using the following decomposition:


                ||X_c_i - Y_c_j||² = ||X_c_i||² - 2 X_c_i.Y_c_j^T + ||Y_c_j||²


    This helper class is in charge of wrapping the common logic to compute
    the middle term `- 2 X_c_i.Y_c_j^T` with a call to GEMM, which has a high
    arithmetic intensity.
    """

    def __init__(self,
        const DTYPE_t[:, ::1] X,
        const DTYPE_t[:, ::1] Y,
        ITYPE_t effective_n_threads,
        ITYPE_t chunks_n_threads,
        ITYPE_t dist_middle_terms_chunks_size,
        ITYPE_t n_features,
        ITYPE_t chunk_size,
    ):
        self.X = X
        self.Y = Y
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

    cdef DTYPE_t * _compute_distances_on_chunks(
        self,
        ITYPE_t X_start,
        ITYPE_t X_end,
        ITYPE_t Y_start,
        ITYPE_t Y_end,
        ITYPE_t thread_num,
    ) nogil:
        cdef:
            ITYPE_t i, j
            DTYPE_t squared_dist_i_j
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
