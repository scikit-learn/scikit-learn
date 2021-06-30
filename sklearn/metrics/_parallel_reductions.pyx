# cython: language_level=3
# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
# cython: profile=False
# cython: linetrace=False
# cython: initializedcheck=False
# cython: binding=False
# distutils: define_macros=CYTHON_TRACE_NOGIL=0

import numpy as np
cimport numpy as np
cimport openmp

from libc.math cimport sqrt
from libc.stdlib cimport free, malloc

from cython.parallel cimport parallel, prange

# from ..neighbors._dist_metrics cimport DistanceMetric

DEF CHUNK_SIZE = 256  # number of vectors

DEF MIN_CHUNK_SAMPLES = 20

DEF FLOAT_INF = 1e36

from ..utils._cython_blas cimport (
  BLAS_Order,
  BLAS_Trans,
  ColMajor,
  NoTrans,
  RowMajor,
  Trans,
  _gemm,
)

from ..utils._heap cimport _simultaneous_sort, _push
from ..utils._openmp_helpers import _openmp_effective_n_threads
from ..utils._typedefs cimport ITYPE_t, DTYPE_t
from ..utils._typedefs import ITYPE, DTYPE


cdef inline DTYPE_t _euclidean_dist(
    DTYPE_t[:, ::1] X,
    DTYPE_t[:, ::1] Y,
    ITYPE_t i,
    ITYPE_t j,
) nogil:
    cdef:
        DTYPE_t dist = 0
        ITYPE_t k
        ITYPE_t upper_unrolled_idx = (X.shape[1] // 4) * 4

    # Unrolling loop to help with vectorisation
    for k in range(0, upper_unrolled_idx, 4):
        dist += (X[i, k] - Y[j, k]) * (X[i, k] - Y[j, k])
        dist += (X[i, k + 1] - Y[j, k + 1]) * (X[i, k + 1] - Y[j, k + 1])
        dist += (X[i, k + 2] - Y[j, k + 2]) * (X[i, k + 2] - Y[j, k + 2])
        dist += (X[i, k + 3] - Y[j, k + 3]) * (X[i, k + 3] - Y[j, k + 3])

    for k in range(upper_unrolled_idx, X.shape[1]):
        dist += (X[i, k] - Y[j, k]) * (X[i, k] - Y[j, k])

    return sqrt(dist)

cdef int _exact_euclidean_dist(
    DTYPE_t[:, ::1] X,                  # IN
    DTYPE_t[:, ::1] Y,                  # IN
    ITYPE_t[:, ::1] Y_indices,          # IN
    ITYPE_t n_threads,                  # IN
    DTYPE_t[:, ::1] distances,          # OUT
) nogil:
    """
    Compute exact pairwise euclidean distances in parallel.

    The pairwise distances considered are X vectors
    and a subset of Y given for each row if X given in
    Y_indices.

    Notes: the body of this function could have been inlined,
    but we use a function to have a cdef nogil context.
    """
    cdef:
        ITYPE_t i, k

    for i in prange(X.shape[0], schedule='static',
                    nogil=True, num_threads=n_threads):
        for k in range(Y_indices.shape[1]):
            distances[i, k] = _euclidean_dist(X, Y, i,
                                              Y_indices[i, k])


cdef class ParallelReduction:
    """Abstract class to computes a reduction of a set of
    vectors (rows) of X on another set of vectors (rows) of Y

    The implementation of the reduction is done parallelised
    on chunks whose size can be set using ``chunk_size``.
    Parameters
    ----------
    X: ndarray of shape (n, d)
        Rows represent vectors
    Y: ndarray of shape (m, d)
        Rows represent vectors
    chunk_size: int
        The number of vectors per chunk.
    """

    cdef:
        ITYPE_t effective_omp_n_thread

        ITYPE_t k, d, sf, si
        ITYPE_t n_samples_chunk, chunk_size

        ITYPE_t n_Y, Y_n_samples_chunk, Y_n_samples_rem
        ITYPE_t n_X, X_n_samples_chunk, X_n_samples_rem

        # Counting remainder chunk in total number of chunks
        ITYPE_t Y_n_chunks, X_n_chunks, num_threads

        DTYPE_t[:, ::1] X
        DTYPE_t[:, ::1] Y

        # TODO: needs to move DistanceMetric
        # from neighbors to be able to use them
        # some adaptation
        # DistanceMetric distance_metric

    def __cinit__(self):
        # Initializing memory view to prevent memory errors and seg-faults
        # in rare cases where __init__ is not called
        self.X = np.empty((1, 1), dtype=DTYPE, order='c')
        self.Y = np.empty((1, 1), dtype=DTYPE, order='c')

    def __init__(self,
                  DTYPE_t[:, ::1] X,
                  DTYPE_t[:, ::1] Y,
                  ITYPE_t k,
                  ITYPE_t chunk_size = CHUNK_SIZE,
    ):
        cdef:
            ITYPE_t X_n_full_chunks, Y_n_full_chunks
        self.X = X
        self.Y = Y

        # TODO: use proper internals checks of scikit-learn
        assert X.shape[1] == Y.shape[1], (
            f"Vectors of X and Y must have the same "
            f"number of dimensions but are respectively "
            f"{X.shape[1]}-dimensional and {Y.shape[1]}-dimensional."
        )

        self.k = k
        self.d = X.shape[1]
        self.sf = sizeof(DTYPE_t)
        self.si = sizeof(ITYPE_t)
        self.chunk_size = chunk_size
        self.n_samples_chunk = max(MIN_CHUNK_SAMPLES, chunk_size)

        self.n_Y = Y.shape[0]
        self.Y_n_samples_chunk = min(self.n_Y, self.n_samples_chunk)
        Y_n_full_chunks = self.n_Y // self.Y_n_samples_chunk
        self.Y_n_samples_rem = self.n_Y % self.Y_n_samples_chunk

        self.n_X = X.shape[0]
        self.X_n_samples_chunk = min(self.n_X, self.n_samples_chunk)
        X_n_full_chunks = self.n_X // self.X_n_samples_chunk
        self.X_n_samples_rem = self.n_X % self.X_n_samples_chunk

        # Counting remainder chunk in total number of chunks
        self.Y_n_chunks = Y_n_full_chunks + (
            self.n_Y != (Y_n_full_chunks * self.Y_n_samples_chunk)
        )

        self.X_n_chunks = X_n_full_chunks + (
            self.n_X != (X_n_full_chunks * self.X_n_samples_chunk)
        )

        self.effective_omp_n_thread = _openmp_effective_n_threads()


    cdef int _reduce_on_chunks(self,
        DTYPE_t[:, ::1] X,                  # IN
        DTYPE_t[:, ::1] Y,                  # IN
        ITYPE_t X_start,
        ITYPE_t X_end,
        ITYPE_t Y_start,
        ITYPE_t Y_end,
        ITYPE_t thread_num,
    ) nogil except -1:
        """ Abstract method: Sub-classes implemented the reduction
        on a pair of chunks"""
        return -1

cdef class ArgKmin(ParallelReduction):

    cdef:
        DTYPE_t ** dist_middle_terms_chunks
        DTYPE_t ** heaps_red_distances_chunks
        ITYPE_t ** heaps_indices_chunks

    def __init__(self,
                  DTYPE_t[:, ::1] X,
                  DTYPE_t[:, ::1] Y,
                  ITYPE_t k,
                  ITYPE_t chunk_size = CHUNK_SIZE,
    ):
        ParallelReduction.__init__(self, X, Y, k)

        self.dist_middle_terms_chunks = <DTYPE_t **> malloc(sizeof(DTYPE_t *) * self.effective_omp_n_thread)
        self.heaps_red_distances_chunks = <DTYPE_t **> malloc(sizeof(DTYPE_t *) * self.effective_omp_n_thread)
        self.heaps_indices_chunks = <ITYPE_t **> malloc(sizeof(ITYPE_t *) * self.effective_omp_n_thread)

    def __dealloc__(self):
        if self.heaps_indices_chunks is not NULL:
            free(self.heaps_indices_chunks)
        else:
            raise RuntimeError("Trying to free heaps_indices_chunks which is NULL")

        if self.heaps_red_distances_chunks is not NULL:
            free(self.heaps_red_distances_chunks)
        else:
            raise RuntimeError("Trying to free heaps_red_distances_chunks which is NULL")

        if self.dist_middle_terms_chunks is not NULL:
            free(self.dist_middle_terms_chunks)
        else:
            raise RuntimeError("Trying to free dist_middle_terms_chunks which is NULL")

    cdef int _reduce_on_chunks(self,
        DTYPE_t[:, ::1] X,                  # IN
        DTYPE_t[:, ::1] Y,                  # IN
        ITYPE_t X_start,
        ITYPE_t X_end,
        ITYPE_t Y_start,
        ITYPE_t Y_end,
        ITYPE_t thread_num,
    ) nogil except -1:
        cdef:
            ITYPE_t i, j
            DTYPE_t[:, ::1] X_c = X[X_start:X_end, :]
            DTYPE_t[:, ::1] Y_c = Y[Y_start:Y_end, :]
            ITYPE_t k = self.k
            DTYPE_t *dist_middle_terms = self.dist_middle_terms_chunks[thread_num]
            DTYPE_t *heaps_red_distances = self.heaps_red_distances_chunks[thread_num]
            ITYPE_t *heaps_indices = self.heaps_indices_chunks[thread_num]

            ITYPE_t n_x = X_end - X_start
            ITYPE_t n_y = Y_end - Y_start

        for i in range(X_c.shape[0]):
            for j in range(Y_c.shape[0]):
                _push(heaps_red_distances + i * self.k,
                      heaps_indices + i * self.k,
                      k,
                      0,
                      # TODO: needs to move DistanceMetric
                      # from neighbors to be able to use them
                      # some adaptation
                      # self.distance_metric.rdist(&X_c[i, 0],
                      #                           &Y_c[j, 0],
                      #                           self.d),
                      Y_start + j)

        return 0

    cdef int _parallel_on_X(self,
        ITYPE_t[:, ::1] argkmin_indices,
        DTYPE_t[:, ::1] argkmin_red_distances,
    ) nogil:
        """Computes the argkmin of each vector (row) of X on Y
        by parallelising computation on chunks of X.
        """
        cdef:
            ITYPE_t Y_start, Y_end, X_start, X_end, X_chunk_idx, Y_chunk_idx, idx, jdx
            ITYPE_t num_threads = min(self.X_n_chunks, self.effective_omp_n_thread)
            ITYPE_t thread_num

            # in bytes
            ITYPE_t size_dist_middle_terms = self.Y_n_samples_chunk * self.X_n_samples_chunk * self.sf
            ITYPE_t heap_size = self.X_n_samples_chunk * self.k * self.sf

        with nogil, parallel(num_threads=num_threads):
            # Thread local buffers
            thread_num = openmp.omp_get_thread_num()
            # Temporary buffer for the -2 * X_c.dot(Y_c.T) term
            self.dist_middle_terms_chunks[thread_num] = <DTYPE_t*> malloc(size_dist_middle_terms)
            self.heaps_red_distances_chunks[thread_num] = <DTYPE_t*> malloc(heap_size)

            for X_chunk_idx in prange(self.X_n_chunks, schedule='static'):
                # We reset the heap between X chunks (memset can't be used here)
                for idx in range(self.X_n_samples_chunk * self.k):
                    self.heaps_red_distances_chunks[thread_num][idx] = FLOAT_INF

                X_start = X_chunk_idx * self.X_n_samples_chunk
                if X_chunk_idx == self.X_n_chunks - 1 and self.X_n_samples_rem > 0:
                    X_end = X_start + self.X_n_samples_rem
                else:
                    X_end = X_start + self.X_n_samples_chunk

                # Referencing the thread-local heaps via the thread-scope pointer
                # of pointers attached to the instance
                self.heaps_indices_chunks[thread_num] = &argkmin_indices[X_start, 0]

                for Y_chunk_idx in range(self.Y_n_chunks):
                    Y_start = Y_chunk_idx * self.Y_n_samples_chunk
                    if Y_chunk_idx == self.Y_n_chunks - 1 and self.Y_n_samples_rem > 0:
                        Y_end = Y_start + self.Y_n_samples_rem
                    else:
                        Y_end = Y_start + self.Y_n_samples_chunk

                    self._reduce_on_chunks(
                        self.X,
                        self.Y,
                        X_start, X_end,
                        Y_start, Y_end,
                        thread_num,
                    )

                # Sorting indices so that the closests' come first.
                for idx in range(X_end - X_start):
                    _simultaneous_sort(
                        self.heaps_red_distances_chunks[thread_num] + idx * self.k,
                        &argkmin_indices[X_start + idx, 0],
                        self.k
                    )

            # end: for X_chunk_idx
            free(self.dist_middle_terms_chunks[thread_num])
            free(self.heaps_red_distances_chunks[thread_num])

        # end: with nogil, parallel
        return self.X_n_chunks


    cdef int _parallel_on_Y(self,
        ITYPE_t[:, ::1] argkmin_indices,          # OUT
        DTYPE_t[:, ::1] argkmin_red_distances,   # OUT
    ) nogil:
        """Computes the argkmin of each vector (row) of X on Y
        by parallelising computation on chunks of Y.

        This parallelisation strategy is more costly (as we need
        extra heaps and synchronisation), yet it is useful in
        most contexts.
        """
        cdef:
            ITYPE_t Y_start, Y_end, X_start, X_end, X_chunk_idx, Y_chunk_idx, idx, jdx
            ITYPE_t num_threads = min(self.X_n_chunks, self.effective_omp_n_thread)
            ITYPE_t thread_num

            # in bytes
            ITYPE_t size_dist_middle_terms = self.Y_n_samples_chunk * self.X_n_samples_chunk * self.sf
            ITYPE_t int_heap_size = self.X_n_samples_chunk * self.k * self.si
            ITYPE_t float_heap_size = self.X_n_samples_chunk * self.k * self.sf

        for X_chunk_idx in range(self.X_n_chunks):
            X_start = X_chunk_idx * self.X_n_samples_chunk
            if X_chunk_idx == self.X_n_chunks - 1 and self.X_n_samples_rem > 0:
                X_end = X_start + self.X_n_samples_rem
            else:
                X_end = X_start + self.X_n_samples_chunk

            with nogil, parallel(num_threads=num_threads):
                # Thread local buffers
                thread_num = openmp.omp_get_thread_num()

                # Temporary buffer for the -2 * X_c.dot(Y_c.T) term
                self.dist_middle_terms_chunks[thread_num] = <DTYPE_t*> malloc(size_dist_middle_terms)
                self.heaps_red_distances_chunks[thread_num] = <DTYPE_t*> malloc(float_heap_size)

                # As chunks of X are shared across threads, so must their
                # heaps. To solve this, each thread has its own locals
                # heaps which are then synchronised back in the main ones.
                self.heaps_indices_chunks[thread_num] = <ITYPE_t*> malloc(int_heap_size)

                # Initialising heaps (memset can't be used here)
                for idx in range(self.X_n_samples_chunk * self.k):
                    self.heaps_red_distances_chunks[thread_num][idx] = FLOAT_INF
                    self.heaps_indices_chunks[thread_num][idx] = -1

                for Y_chunk_idx in prange(self.Y_n_chunks, schedule='static'):
                    Y_start = Y_chunk_idx * self.Y_n_samples_chunk
                    if Y_chunk_idx == self.Y_n_chunks - 1 \
                        and self.Y_n_samples_rem > 0:
                        Y_end = Y_start + self.Y_n_samples_rem
                    else:
                        Y_end = Y_start + self.Y_n_samples_chunk


                    self._reduce_on_chunks(
                        self.X,
                        self.Y,
                        X_start, X_end,
                        Y_start, Y_end,
                        thread_num,
                    )

                # end: for Y_chunk_idx
                with gil:
                    # Synchronising the thread local heaps
                    # with the main heaps
                    for idx in range(X_end - X_start):
                        for jdx in range(self.k):
                            _push(
                                &argkmin_red_distances[X_start + idx, 0],
                                &argkmin_indices[X_start + idx, 0],
                                self.k,
                                self.heaps_red_distances_chunks[thread_num][idx * self.k + jdx],
                                self.heaps_indices_chunks[thread_num][idx * self.k + jdx],
                            )

            # end: with nogil, parallel

            # Sorting indices of the argkmin for each query vector of X
            for idx in prange(self.n_X, schedule='static',
                              nogil=True, num_threads=num_threads):
                _simultaneous_sort(
                    &argkmin_red_distances[idx, 0],
                    &argkmin_indices[idx, 0],
                    self.k,
                )
            # end: prange

        # end: for X_chunk_idx
        return self.Y_n_chunks


    # Python interface
    def compute(self,
           str strategy = "auto",
           bint return_distance = False
    ):
        """Computes the reduction of vectors (rows) of X on Y.

        strategy: str, {'auto', 'parallel_on_X', 'parallel_on_Y'}
            The chunking strategy defining which dataset
            parallelisation are made on.

             - 'parallel_on_X' is embarassingly parallel but
            is less used in practice.
             - 'parallel_on_Y' comes with synchronisation but
            is more useful in practice.
             -'auto' relies on a simple heuristic to choose
            between 'parallel_on_X' and 'parallel_on_Y'.

        return_distance: boolean
            Return distances between each X vector and its
            argkmin if set to True.

        Returns
        -------
        distances: ndarray of shape (n, k)
            Distances between each X vector and its argkmin
            in Y. Only returned if ``return_distance=True``.

        indices: ndarray of shape (n, k)
            Indices of each X vector argkmin in Y.
        """
        cdef:
            ITYPE_t n_X = self.X.shape[0]
            ITYPE_t[:, ::1] argkmin_indices = np.full((n_X, self.k), 0,
                                                   dtype=ITYPE)
            DTYPE_t[:, ::1] argkmin_distances = np.full((n_X, self.k),
                                                      FLOAT_INF,
                                                      dtype=DTYPE)

        if strategy == 'auto':
            # This is a simple heuristic whose constant for the
            # comparison has been chosen based on experiments.
            if 4 * self.chunk_size * self.effective_omp_n_thread < n_X:
                strategy = 'parallel_on_X'
            else:
                strategy = 'parallel_on_Y'

        if strategy == 'parallel_on_Y':
            self._parallel_on_Y(
                argkmin_indices, argkmin_distances
            )
        elif strategy == 'parallel_on_X':
            self._parallel_on_X(
                argkmin_indices, argkmin_distances
            )
        else:
            raise RuntimeError(f"strategy '{strategy}' not supported.")

        if return_distance:
            # We need to recompute distances because we relied on
            # reduced distances using _gemm, which are missing a
            # term for squared norms and which are not the most
            # precise (catastrophic cancellation might have happened).
            _exact_euclidean_dist(self.X, self.Y, argkmin_indices,
                                  self.effective_omp_n_thread,
                                  argkmin_distances)
            return (np.asarray(argkmin_distances),
                    np.asarray(argkmin_indices))

        return np.asarray(argkmin_indices)

cdef class FastSquaredEuclideanArgKmin(ArgKmin):

    cdef:
        DTYPE_t[::1] Y_sq_norms

    def __init__(self,
                  DTYPE_t[:, ::1] X,
                  DTYPE_t[:, ::1] Y,
                  ITYPE_t k,
                  ITYPE_t chunk_size = CHUNK_SIZE,
    ):
        ArgKmin.__init__(self, X, Y, k)
        self.Y_sq_norms = np.einsum('ij,ij->i', self.Y, self.Y)


    cdef int _reduce_on_chunks(self,
        DTYPE_t[:, ::1] X,                  # IN
        DTYPE_t[:, ::1] Y,                  # IN
        ITYPE_t X_start,
        ITYPE_t X_end,
        ITYPE_t Y_start,
        ITYPE_t Y_end,
        ITYPE_t thread_num,
    ) nogil except -1:
        """
        Critical part of the computation of pairwise distances.

        "Fast Squared Euclidean" distances strategy relying
        on the gemm-trick.
        """
        cdef:
            ITYPE_t i, j
            DTYPE_t[:, ::1] X_c = X[X_start:X_end, :]
            DTYPE_t[:, ::1] Y_c = Y[Y_start:Y_end, :]
            ITYPE_t k = self.k
            DTYPE_t *dist_middle_terms = self.dist_middle_terms_chunks[thread_num]
            DTYPE_t *heaps_red_distances = self.heaps_red_distances_chunks[thread_num]
            ITYPE_t *heaps_indices = self.heaps_indices_chunks[thread_num]

        # Instead of computing the full pairwise squared distances matrix,
        # ||X_c - Y_c||² = ||X_c||² - 2 X_c.Y_c^T + ||Y_c||²,
        # we only need to store the - 2 X_c.Y_c^T + ||Y_c||²
        # term since the argmin for a given sample X_c^{i} does not depend on
        # ||X_c^{i}||²

        # Careful: LDA, LDB and LDC are given for F-ordered arrays.
        # Here, we use their counterpart values as indicated in the documentation.
        # See the documentation of parameters here:
        # https://www.netlib.org/lapack/explore-html/db/dc9/group__single__blas__level3_gafe51bacb54592ff5de056acabd83c260.html
        #
        # dist_middle_terms = -2 * X_c.dot(Y_c.T)
        _gemm(RowMajor, NoTrans, Trans,
              X_c.shape[0], Y_c.shape[0], X_c.shape[1],
              -2.0,
              &X_c[0, 0], X_c.shape[1],
              &Y_c[0, 0], X_c.shape[1], 0.0,
              dist_middle_terms, Y_c.shape[0])

        # Computing argmins here
        for i in range(X_c.shape[0]):
            for j in range(Y_c.shape[0]):
                _push(heaps_red_distances + i * k,
                      heaps_indices + i * k,
                      k,
                      # reduced distance: - 2 X_c_i.Y_c_j^T + ||Y_c_j||²
                      dist_middle_terms[i * Y_c.shape[0] + j] + self.Y_sq_norms[j + Y_start],
                      j + Y_start)
        return 0
