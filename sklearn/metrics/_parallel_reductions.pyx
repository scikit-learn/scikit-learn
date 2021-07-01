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

from ._dist_metrics cimport DistanceMetric
from ._dist_metrics import METRIC_MAPPING
from ..utils import check_array

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


cdef class ParallelReduction:
    """Abstract class to computes a reduction of a set of
    vectors (rows) of X on another set of vectors (rows) of Y.

    The implementation of the reduction is done parallelized
    on chunks whose size can be set using ``chunk_size``.
    Parameters
    ----------
    X: ndarray of shape (n, d)
        Rows represent vectors
    Y: ndarray of shape (m, d)
        Rows represent vectors
    distance_metric: DistanceMetric
        The distance to use
    chunk_size: int
        The number of vectors per chunk
    """

    cdef:
        const DTYPE_t[:, ::1] X  # shape: (n_X, d)
        const DTYPE_t[:, ::1] Y  # shape: (n_Y, d)

        DistanceMetric distance_metric

        ITYPE_t effective_omp_n_thread
        ITYPE_t n_samples_chunk, chunk_size

        ITYPE_t d

        # dtypes sizes
        ITYPE_t sf, si

        ITYPE_t n_X, X_n_samples_chunk, X_n_chunks, X_n_samples_rem
        ITYPE_t n_Y, Y_n_samples_chunk, Y_n_chunks, Y_n_samples_rem

    @classmethod
    def valid_metrics(cls):
        return {*METRIC_MAPPING.keys()}

    def __cinit__(self):
        # Initializing memory view to prevent memory errors and seg-faults
        # in rare cases where __init__ is not called
        self.X = np.empty((1, 1), dtype=DTYPE, order='c')
        self.Y = np.empty((1, 1), dtype=DTYPE, order='c')

    def __init__(self,
                 X,
                 Y,
                 DistanceMetric distance_metric,
                 ITYPE_t chunk_size = CHUNK_SIZE,
    ):
        cdef:
            ITYPE_t X_n_full_chunks, Y_n_full_chunks

        self.effective_omp_n_thread = _openmp_effective_n_threads()

        self.X = check_array(X, dtype=DTYPE)
        self.Y = check_array(Y, dtype=DTYPE)

        assert X.shape[1] == Y.shape[1], "Vectors of X and Y must have the " \
                                         "same dimension but currently are " \
                                         f"respectively {X.shape[1]}-dimensional " \
                                         f"and {Y.shape[1]}-dimensional."
        distance_metric._validate_data(X)
        distance_metric._validate_data(Y)

        self.d = X.shape[1]
        self.sf = sizeof(DTYPE_t)
        self.si = sizeof(ITYPE_t)
        self.chunk_size = chunk_size
        self.n_samples_chunk = max(MIN_CHUNK_SAMPLES, chunk_size)

        self.distance_metric = distance_metric

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

    def __dealloc__(self):
        pass

    cdef void _on_X_parallel_init(self,
            ITYPE_t thread_num,
    ) nogil:
        return

    cdef void _on_X_parallel_finalize(self,
            ITYPE_t thread_num
    ) nogil:
        return

    cdef void _on_X_prange_iter_init(self,
            ITYPE_t thread_num,
            ITYPE_t X_chunk_idx,
            ITYPE_t X_start,
            ITYPE_t X_end,
    ) nogil:
        return

    cdef void _on_X_prange_iter_finalize(self,
            ITYPE_t thread_num,
            ITYPE_t X_chunk_idx,
            ITYPE_t X_start,
            ITYPE_t X_end,
    ) nogil:
        return

    cdef void _parallel_on_X(self) nogil:
        """Computes the reduction of each vector (row) of X on Y
        by parallelizing computation on chunks of X.

        Private datastructures are modified internally by threads.

        Private template methods can be implemented on subclasses to
        interact with those datastructures at various stages.
        """
        cdef:
            ITYPE_t Y_start, Y_end, X_start, X_end, X_chunk_idx, Y_chunk_idx
            ITYPE_t num_threads = min(self.X_n_chunks, self.effective_omp_n_thread)
            ITYPE_t thread_num

        with nogil, parallel(num_threads=num_threads):
            thread_num = openmp.omp_get_thread_num()

            # Allocating thread local datastructures
            self._on_X_parallel_init(thread_num)

            for X_chunk_idx in prange(self.X_n_chunks, schedule='static'):
                X_start = X_chunk_idx * self.X_n_samples_chunk
                if X_chunk_idx == self.X_n_chunks - 1 and self.X_n_samples_rem > 0:
                    X_end = X_start + self.X_n_samples_rem
                else:
                    X_end = X_start + self.X_n_samples_chunk

                # Reinitializing thread local datastructures for the new X chunk
                self._on_X_prange_iter_init(thread_num, X_chunk_idx, X_start, X_end)

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

                # Adjusting thread local datastructures on the full pass on Y
                self._on_X_prange_iter_finalize(thread_num, X_chunk_idx, X_start, X_end)

            # end: for X_chunk_idx

            # Deallocating thread local datastructures
            self._on_X_parallel_finalize(thread_num)

        # end: with nogil, parallel
        return

    cdef void _on_Y_parallel_init(self,
        ITYPE_t thread_num,
    ) nogil:
        return

    cdef void _on_Y_parallel_finalize(self,
        ITYPE_t thread_num,
        ITYPE_t X_chunk_idx,
        ITYPE_t X_start,
        ITYPE_t X_end,
    ) nogil:
        return

    cdef void _on_Y_finalize(self,
        ITYPE_t thread_num,
    ) nogil:
        return

    cdef void _parallel_on_Y(self) nogil:
        """Computes the argkmin of each vector (row) of X on Y
        by parallelizing computation on chunks of Y.

        Private datastructures are modified internally by threads.

        Private template methods can be implemented on subclasses to
        interact with those datastructures at various stages.
        """
        cdef:
            ITYPE_t Y_start, Y_end, X_start, X_end, X_chunk_idx, Y_chunk_idx
            ITYPE_t num_threads = min(self.X_n_chunks, self.effective_omp_n_thread)
            ITYPE_t thread_num

        for X_chunk_idx in range(self.X_n_chunks):
            X_start = X_chunk_idx * self.X_n_samples_chunk
            if X_chunk_idx == self.X_n_chunks - 1 and self.X_n_samples_rem > 0:
                X_end = X_start + self.X_n_samples_rem
            else:
                X_end = X_start + self.X_n_samples_chunk

            with nogil, parallel(num_threads=num_threads):
                # Thread local buffers
                thread_num = openmp.omp_get_thread_num()

                # Allocating thread local datastructures
                self._on_Y_parallel_init(thread_num)

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
                # end: prange

                # Synchronizing thread local datastructures with the main ones
                # This can potentially block
                self._on_Y_parallel_finalize(thread_num, X_chunk_idx, X_start, X_end)
            # end: with nogil, parallel

        # end: for X_chunk_idx
        # Adjusting main datastructures before returning
        self._on_Y_finalize(num_threads)
        return

    cdef int _reduce_on_chunks(self,
        const DTYPE_t[:, ::1] X,
        const DTYPE_t[:, ::1] Y,
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
    """Computes the argkmin of vectors (rows) of a set of
    vectors (rows) of X on another set of vectors (rows) of Y.

    The implementation is parallelized on chunks whose size can
    be set using ``chunk_size``.

    Parameters
    ----------
    X: ndarray of shape (n, d)
        Rows represent vectors
    Y: ndarray of shape (m, d)
        Rows represent vectors
    distance_metric: DistanceMetric
        The distance to use
    k: int
        The k for the argkmin reduction
    chunk_size: int
        The number of vectors per chunk
    """

    cdef:
        ITYPE_t k

        DTYPE_t ** heaps_approx_distances_chunks
        ITYPE_t ** heaps_indices_chunks

        ITYPE_t[:, ::1] argkmin_indices
        DTYPE_t[:, ::1] argkmin_distances

    @classmethod
    def valid_metrics(cls):
        return {"fast_sqeuclidean", *METRIC_MAPPING.keys()}

    @classmethod
    def get_for(cls,
                X,
                Y,
                ITYPE_t k,
                str metric="fast_sqeuclidean",
                ITYPE_t chunk_size=CHUNK_SIZE,
                dict metric_kwargs=dict(),
        ):
        if metric == "fast_sqeuclidean":
            return FastSquaredEuclideanArgKmin(X=X, Y=Y, k=k, chunk_size=chunk_size)
        return ArgKmin(X=X, Y=Y,
                       distance_metric=DistanceMetric.get_metric(metric, **metric_kwargs),
                       k=k,
                       chunk_size=chunk_size)

    def __init__(self,
                 X,
                 Y,
                 DistanceMetric distance_metric,
                 ITYPE_t k,
                 ITYPE_t chunk_size = CHUNK_SIZE,
    ):
        ParallelReduction.__init__(self, X, Y, distance_metric, chunk_size)

        self.k = k

        # Results returned by ArgKmin.compute
        self.argkmin_indices = np.full((self.n_X, self.k), 0, dtype=ITYPE)
        self.argkmin_distances = np.full((self.n_X, self.k), FLOAT_INF, dtype=DTYPE)

        # Temporary datastructures used in threads
        self.heaps_approx_distances_chunks = <DTYPE_t **> malloc(sizeof(DTYPE_t *) * self.effective_omp_n_thread)
        self.heaps_indices_chunks = <ITYPE_t **> malloc(sizeof(ITYPE_t *) * self.effective_omp_n_thread)

    def __dealloc__(self):
        ParallelReduction.__dealloc__(self)
        if self.heaps_indices_chunks is not NULL:
            free(self.heaps_indices_chunks)
        else:
            raise RuntimeError("Trying to free heaps_indices_chunks which is NULL")

        if self.heaps_approx_distances_chunks is not NULL:
            free(self.heaps_approx_distances_chunks)
        else:
            raise RuntimeError("Trying to free heaps_approx_distances_chunks which is NULL")

    cdef int _reduce_on_chunks(self,
        const DTYPE_t[:, ::1] X,
        const DTYPE_t[:, ::1] Y,
        ITYPE_t X_start,
        ITYPE_t X_end,
        ITYPE_t Y_start,
        ITYPE_t Y_end,
        ITYPE_t thread_num,
    ) nogil except -1:
        cdef:
            ITYPE_t i, j
            const DTYPE_t[:, ::1] X_c = X[X_start:X_end, :]
            const DTYPE_t[:, ::1] Y_c = Y[Y_start:Y_end, :]
            ITYPE_t k = self.k
            DTYPE_t *heaps_approx_distances = self.heaps_approx_distances_chunks[thread_num]
            ITYPE_t *heaps_indices = self.heaps_indices_chunks[thread_num]

            ITYPE_t n_x = X_end - X_start
            ITYPE_t n_y = Y_end - Y_start

        for i in range(X_c.shape[0]):
            for j in range(Y_c.shape[0]):
                _push(heaps_approx_distances + i * self.k,
                      heaps_indices + i * self.k,
                      k,
                      self.distance_metric.rdist(&X_c[i, 0],
                                                 &Y_c[j, 0],
                                                 self.d),
                      Y_start + j)

        return 0

    cdef void _on_X_parallel_init(self,
            ITYPE_t thread_num,
    ) nogil:
        cdef:
            # in bytes
            ITYPE_t heap_size = self.X_n_samples_chunk * self.k * self.sf

        # Temporary buffer for the -2 * X_c.dot(Y_c.T) term
        self.heaps_approx_distances_chunks[thread_num] = <DTYPE_t *> malloc(heap_size)

    cdef void _on_X_prange_iter_init(self,
            ITYPE_t thread_num,
            ITYPE_t X_chunk_idx,
            ITYPE_t X_start,
            ITYPE_t X_end,
    ) nogil:

        # We reset the heap between X chunks (memset can't be used here)
        for idx in range(self.X_n_samples_chunk * self.k):
            self.heaps_approx_distances_chunks[thread_num][idx] = FLOAT_INF

        # Referencing the thread-local heaps via the thread-scope pointer
        # of pointers attached to the instance
        self.heaps_indices_chunks[thread_num] = &self.argkmin_indices[X_start, 0]

    cdef void _on_X_prange_iter_finalize(self,
            ITYPE_t thread_num,
            ITYPE_t X_chunk_idx,
            ITYPE_t X_start,
            ITYPE_t X_end,
    ) nogil:
        cdef:
            ITYPE_t idx, jdx

        # Sorting indices of the argkmin for each query vector of X
        for idx in range(X_end - X_start):
            _simultaneous_sort(
                self.heaps_approx_distances_chunks[thread_num] + idx * self.k,
                &self.argkmin_indices[X_start + idx, 0],
                self.k
            )

    cdef void _on_X_parallel_finalize(self,
            ITYPE_t thread_num
    ) nogil:
        free(self.heaps_approx_distances_chunks[thread_num])

    cdef void _on_Y_parallel_init(self,
            ITYPE_t thread_num,
    ) nogil:
        cdef:
            # in bytes
            ITYPE_t int_heap_size = self.X_n_samples_chunk * self.k * self.si
            ITYPE_t float_heap_size = self.X_n_samples_chunk * self.k * self.sf

        self.heaps_approx_distances_chunks[thread_num] = <DTYPE_t *> malloc(float_heap_size)

        # As chunks of X are shared across threads, so must their
        # heaps. To solve this, each thread has its own locals
        # heaps which are then synchronised back in the main ones.
        self.heaps_indices_chunks[thread_num] = <ITYPE_t *> malloc(int_heap_size)

        # Initialising heaps (memset can't be used here)
        for idx in range(self.X_n_samples_chunk * self.k):
            self.heaps_approx_distances_chunks[thread_num][idx] = FLOAT_INF
            self.heaps_indices_chunks[thread_num][idx] = -1

    cdef void _on_Y_parallel_finalize(self,
            ITYPE_t thread_num,
            ITYPE_t X_chunk_idx,
            ITYPE_t X_start,
            ITYPE_t X_end,
    ) nogil:
        cdef:
            ITYPE_t idx, jdx
        with gil:
            # Synchronising the thread local heaps
            # with the main heaps
            for idx in range(X_end - X_start):
                for jdx in range(self.k):
                    _push(
                        &self.argkmin_distances[X_start + idx, 0],
                        &self.argkmin_indices[X_start + idx, 0],
                        self.k,
                        self.heaps_approx_distances_chunks[thread_num][idx * self.k + jdx],
                        self.heaps_indices_chunks[thread_num][idx * self.k + jdx],
                    )

        free(self.heaps_approx_distances_chunks[thread_num])
        free(self.heaps_indices_chunks[thread_num])

    cdef void _on_Y_finalize(self,
            ITYPE_t thread_num,
    ) nogil:
        cdef:
            ITYPE_t num_threads = min(self.X_n_chunks, self.effective_omp_n_thread)
            ITYPE_t idx

        # Sorting indices of the argkmin for each query vector of X
        for idx in prange(self.n_X, schedule='static',
                          nogil=True, num_threads=num_threads):
            _simultaneous_sort(
                &self.argkmin_distances[idx, 0],
                &self.argkmin_indices[idx, 0],
                self.k,
            )
        return

    cdef void _exact_distances(self,
        ITYPE_t[:, ::1] Y_indices,  # IN
        DTYPE_t[:, ::1] distances,  # IN/OUT
    ) nogil:
        """Convert reduced distances to pairwise distances in parallel."""
        cdef:
            ITYPE_t i, j

        for i in prange(self.n_X, schedule='static', nogil=True,
                        num_threads=self.effective_omp_n_thread):
            for j in range(self.k):
                distances[i, j] = self.distance_metric.dist(&self.X[i, 0],
                                                 &self.Y[Y_indices[i, j], 0],
                                                 self.d)

    # Python interface
    def compute(self,
           str strategy = "auto",
           bint return_distance = False
    ):
        """Computes the reduction of vectors (rows) of X on Y.

        strategy: str, {'auto', 'parallel_on_X', 'parallel_on_Y'}
            The chunking strategy defining which dataset
            parallelization are made on.

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
        if strategy == 'auto':
            # This is a simple heuristic whose constant for the
            # comparison has been chosen based on experiments.
            if 4 * self.chunk_size * self.effective_omp_n_thread < self.n_X:
                strategy = 'parallel_on_X'
            else:
                strategy = 'parallel_on_Y'

        if strategy == 'parallel_on_Y':
            self._parallel_on_Y()
        elif strategy == 'parallel_on_X':
            self._parallel_on_X()
        else:
            raise RuntimeError(f"strategy '{strategy}' not supported.")

        if return_distance:
            # We need to recompute distances because we relied on
            # reduced distances.
            self._exact_distances(self.argkmin_indices, self.argkmin_distances)
            return np.asarray(self.argkmin_distances), np.asarray(self.argkmin_indices)

        return np.asarray(self.argkmin_indices)

cdef class FastSquaredEuclideanArgKmin(ArgKmin):
    """Fast specialized alternative for ArgKmin on
    EuclideanDistance.

    Computes the argkmin of vectors (rows) of a set of
    vectors (rows) of X on another set of vectors (rows) of Y
    using the GEMM-trick.

    This implementation has an superior arithmetic intensity
    and hence running time, but it can suffer from numerical
    instability. We recommend using ArgKmin with
    EuclideanDistance when exact precision is needed.
    """

    cdef:
        DTYPE_t[::1] Y_sq_norms

        # Buffers for GEMM
        DTYPE_t ** dist_middle_terms_chunks

    def __init__(self,
                  X,
                  Y,
                  ITYPE_t k,
                  ITYPE_t chunk_size = CHUNK_SIZE,
    ):
        ArgKmin.__init__(self, X, Y,
                         distance_metric=DistanceMetric.get_metric("euclidean"),
                         k=k,
                         chunk_size=chunk_size)
        self.Y_sq_norms = np.einsum('ij,ij->i', self.Y, self.Y)
        # Temporary datastructures used in threads
        self.dist_middle_terms_chunks = <DTYPE_t **> malloc(sizeof(DTYPE_t *) * self.effective_omp_n_thread)

    def __dealloc__(self):
        ArgKmin.__dealloc__(self)
        if self.dist_middle_terms_chunks is not NULL:
            free(self.dist_middle_terms_chunks)
        else:
            raise RuntimeError("Trying to free dist_middle_terms_chunks which is NULL")

    cdef void _on_X_parallel_init(self,
            ITYPE_t thread_num,
    ) nogil:
        ArgKmin._on_X_parallel_init(self, thread_num)
        # Temporary buffer for the -2 * X_c.dot(Y_c.T) term
        self.dist_middle_terms_chunks[thread_num] = <DTYPE_t *> malloc(
            self.Y_n_samples_chunk * self.X_n_samples_chunk * self.sf)

    cdef void _on_X_parallel_finalize(self,
            ITYPE_t thread_num
    ) nogil:
        ArgKmin._on_X_parallel_finalize(self, thread_num)
        free(self.dist_middle_terms_chunks[thread_num])

    cdef void _on_Y_parallel_init(self,
            ITYPE_t thread_num,
    ) nogil:
        ArgKmin._on_Y_parallel_init(self, thread_num)
        # Temporary buffer for the -2 * X_c.dot(Y_c.T) term
        self.dist_middle_terms_chunks[thread_num] = <DTYPE_t *> malloc(
            self.Y_n_samples_chunk * self.X_n_samples_chunk * self.sf)

    cdef void _on_Y_parallel_finalize(self,
            ITYPE_t thread_num,
            ITYPE_t X_chunk_idx,
            ITYPE_t X_start,
            ITYPE_t X_end,
    ) nogil:
        ArgKmin._on_Y_parallel_finalize(self, thread_num, X_chunk_idx, X_start, X_end)
        free(self.dist_middle_terms_chunks[thread_num])

    cdef int _reduce_on_chunks(self,
        const DTYPE_t[:, ::1] X,
        const DTYPE_t[:, ::1] Y,
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
            const DTYPE_t[:, ::1] X_c = X[X_start:X_end, :]
            const DTYPE_t[:, ::1] Y_c = Y[Y_start:Y_end, :]
            ITYPE_t k = self.k
            DTYPE_t *dist_middle_terms = self.dist_middle_terms_chunks[thread_num]
            DTYPE_t *heaps_approx_distances = self.heaps_approx_distances_chunks[thread_num]
            ITYPE_t *heaps_indices = self.heaps_indices_chunks[thread_num]

            # Instead of computing the full pairwise squared distances matrix,
            #
            #      ||X_c - Y_c||² = ||X_c||² - 2 X_c.Y_c^T + ||Y_c||²,
            #
            # we only need to store the
            #                                - 2 X_c.Y_c^T + ||Y_c||²
            #
            # term since the argkmin for a given sample X_c^{i} does not depend on
            # ||X_c^{i}||²
            #
            # This term gets computed efficiently bellow using GEMM from BLAS Level 3.
            #
            # Careful: LDA, LDB and LDC are given for F-ordered arrays in BLAS documentations,
            # for instance:
            # https://www.netlib.org/lapack/explore-html/db/dc9/group__single__blas__level3_gafe51bacb54592ff5de056acabd83c260.html
            #
            # Here, we use their counterpart values to work with C-ordered arrays.
            BLAS_Order order = RowMajor
            BLAS_Trans ta = NoTrans
            BLAS_Trans tb = Trans
            ITYPE_t m = X_c.shape[0]
            ITYPE_t n = Y_c.shape[0]
            ITYPE_t K = X_c.shape[1]
            DTYPE_t alpha = - 2.
            DTYPE_t * A = & X_c[0, 0]
            ITYPE_t lda = X_c.shape[1]
            DTYPE_t * B = & Y_c[0, 0]
            ITYPE_t ldb = X_c.shape[1]
            DTYPE_t beta = 0.
            DTYPE_t * C = dist_middle_terms
            ITYPE_t ldc = Y_c.shape[0]

        # dist_middle_terms = -2 * X_c.dot(Y_c.T)
        _gemm(order, ta, tb, m, n, K, alpha, A, lda, B, ldb, beta, C, ldc)

        # Pushing the distance and their associated indices on heaps
        # which keep tracks of the argkmin.
        for i in range(X_c.shape[0]):
            for j in range(Y_c.shape[0]):
                _push(heaps_approx_distances + i * k,
                      heaps_indices + i * k,
                      k,
                      # approximated distance: - 2 X_c_i.Y_c_j^T + ||Y_c_j||²
                      dist_middle_terms[i * Y_c.shape[0] + j] + self.Y_sq_norms[j + Y_start],
                      j + Y_start)
        return 0
