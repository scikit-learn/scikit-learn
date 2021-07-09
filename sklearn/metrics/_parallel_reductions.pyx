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

np.import_array()

from libc.stdlib cimport free, malloc
from libcpp.vector cimport vector
from cython.operator cimport dereference as deref

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
from ..utils._typedefs cimport ITYPE_t, DTYPE_t, DITYPE_t
from ..utils._typedefs cimport ITYPECODE, DTYPECODE
from ..utils._typedefs import ITYPE, DTYPE


cdef extern from "<algorithm>" namespace "std" nogil:
    OutputIt move[InputIt, OutputIt](InputIt first, InputIt last, OutputIt d_first)

######################
## std::vector to np.ndarray coercion
# TODO: for now using this simple solution: https://stackoverflow.com/a/23873586
# A better solution would make sure of using the same allocator implementations.
# Those implementations depend on the runtimes' allocator which can be different
# in some configuration and thus would make the program crash.

cdef extern from "numpy/arrayobject.h":
    void PyArray_ENABLEFLAGS(np.ndarray arr, int flags)

cdef np.ndarray[DITYPE_t, ndim=1] buffer_to_numpy_array(DITYPE_t * ptr, np.npy_intp size):
    """ Create a numpy ndarray given a buffer and its size. """
    typenum = DTYPECODE if DITYPE_t is DTYPE_t else ITYPECODE
    cdef np.ndarray[DITYPE_t, ndim=1] arr = np.PyArray_SimpleNewFromData(1, &size, typenum, ptr)

    # Makes the numpy array responsible to the life-cycle of its buffer.
    PyArray_ENABLEFLAGS(arr, np.NPY_OWNDATA)
    return arr

# TODO: this got duplicated because type covariance is not support; i.e. the following function
#
#       cdef np.ndarray[object, ndim=1] _coerce_vectors_to_np_nd_arrays(vector[vector[DITYPE_t]] * vecs)
#
# cannot be called dispatched for vector[vector[ITYPE_t]]* and vector[vector[DTYPE_t]]*
cdef np.ndarray[object, ndim=1] _coerce_vectors_to_np_nd_arrays_ITYPE(vector[vector[ITYPE_t]]* vecs):
    cdef ITYPE_t n = deref(vecs).size()
    np_arrays_of_np_arrays = np.empty(n, dtype=np.ndarray)

    for i in range(n):
        np_arrays_of_np_arrays[i] = buffer_to_numpy_array(deref(vecs)[i].data(),
                                                           deref(vecs)[i].size())

    return np_arrays_of_np_arrays

cdef np.ndarray[object, ndim=1] _coerce_vectors_to_np_nd_arrays_DTYPE(vector[vector[DTYPE_t]]* vecs):
    cdef ITYPE_t n = deref(vecs).size()
    np_arrays_of_np_arrays = np.empty(n, dtype=np.ndarray)

    for i in range(n):
        np_arrays_of_np_arrays[i] = buffer_to_numpy_array(deref(vecs)[i].data(),
                                                           deref(vecs)[i].size())

    return np_arrays_of_np_arrays

#####################


cdef class PairwiseDistancesReduction:
    """Abstract class to computes a reduction on pairwise
    distances between a set of vectors (rows) X and another
    set of vectors (rows) of Y.

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

        ITYPE_t n_X, X_n_samples_chunk, X_n_chunks, X_n_samples_rem
        ITYPE_t n_Y, Y_n_samples_chunk, Y_n_chunks, Y_n_samples_rem

    @classmethod
    def valid_metrics(cls):
        excluded = {"pyfunc", "sokalmichener", "matching", "jaccard"}
        return sorted({*METRIC_MAPPING.keys()}.difference(excluded))

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

            # Allocating thread datastructures
            self._on_X_parallel_init(thread_num)

            for X_chunk_idx in prange(self.X_n_chunks, schedule='static'):
                X_start = X_chunk_idx * self.X_n_samples_chunk
                if X_chunk_idx == self.X_n_chunks - 1 and self.X_n_samples_rem > 0:
                    X_end = X_start + self.X_n_samples_rem
                else:
                    X_end = X_start + self.X_n_samples_chunk

                # Reinitializing thread datastructures for the new X chunk
                self._on_X_prange_iter_init(thread_num, X_start, X_end)

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

                # Adjusting thread datastructures on the full pass on Y
                self._on_X_prange_iter_finalize(thread_num, X_start, X_end)

            # end: for X_chunk_idx

            # Deallocating thread datastructures
            self._on_X_parallel_finalize(thread_num)

        # end: with nogil, parallel
        return

    cdef void _parallel_on_Y(self) nogil:
        """Computes the reduction of each vector (row) of X on Y
        by parallelizing computation on chunks of Y.

        Private datastructures are modified internally by threads.

        Private template methods can be implemented on subclasses to
        interact with those datastructures at various stages.
        """
        cdef:
            ITYPE_t Y_start, Y_end, X_start, X_end, X_chunk_idx, Y_chunk_idx
            ITYPE_t num_threads = min(self.Y_n_chunks, self.effective_omp_n_thread)
            ITYPE_t thread_num

        # TODO: put the "with nogil, parallel"-context here
        # Allocating datastructures
        self._on_Y_init(num_threads)

        for X_chunk_idx in range(self.X_n_chunks):
            X_start = X_chunk_idx * self.X_n_samples_chunk
            if X_chunk_idx == self.X_n_chunks - 1 and self.X_n_samples_rem > 0:
                X_end = X_start + self.X_n_samples_rem
            else:
                X_end = X_start + self.X_n_samples_chunk

            with nogil, parallel(num_threads=num_threads):
                thread_num = openmp.omp_get_thread_num()

                # Initializing datastructures used in this thread
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

                # Synchronizing the thread datastructures with the main ones
                # This can potentially block
                self._on_Y_parallel_finalize(thread_num, X_start, X_end)
            # end: with nogil, parallel

        # end: for X_chunk_idx
        # Deallocating temporary datastructures
        # Adjusting main datastructures before returning
        self._on_Y_finalize(num_threads)
        return

    # Placeholder methods which have to be implemented
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

    # Placeholder methods which can be implemented

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
            ITYPE_t X_start,
            ITYPE_t X_end,
    ) nogil:
        return

    cdef void _on_X_prange_iter_finalize(self,
            ITYPE_t thread_num,
            ITYPE_t X_start,
            ITYPE_t X_end,
    ) nogil:
        return

    cdef void _on_Y_init(self,
        ITYPE_t num_threads,
    ) nogil:
        return

    cdef void _on_Y_parallel_init(self,
        ITYPE_t thread_num,
    ) nogil:
        return

    cdef void _on_Y_parallel_finalize(self,
        ITYPE_t thread_num,
        ITYPE_t X_start,
        ITYPE_t X_end,
    ) nogil:
        return

    cdef void _on_Y_finalize(self,
        ITYPE_t num_threads,
    ) nogil:
        return

cdef class ArgKmin(PairwiseDistancesReduction):
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
        return {"fast_sqeuclidean", *PairwiseDistancesReduction.valid_metrics()}

    @classmethod
    def get_for(cls,
                X,
                Y,
                ITYPE_t k,
                str metric="fast_sqeuclidean",
                ITYPE_t chunk_size=CHUNK_SIZE,
                dict metric_kwargs=dict(),
        ):
        # This factory comes to handle specialisation on fast_sqeuclidean.
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
        PairwiseDistancesReduction.__init__(self, X, Y, distance_metric, chunk_size)

        self.k = k

        # Results returned by ArgKmin.compute used as the main heaps
        self.argkmin_indices = np.full((self.n_X, self.k), 0, dtype=ITYPE)
        self.argkmin_distances = np.full((self.n_X, self.k), FLOAT_INF, dtype=DTYPE)

        # Pointers to thread heaps used in threads for `parallel_on_Y` solely
        self.heaps_approx_distances_chunks = <DTYPE_t **> malloc(sizeof(DTYPE_t *) * self.effective_omp_n_thread)
        self.heaps_indices_chunks = <ITYPE_t **> malloc(sizeof(ITYPE_t *) * self.effective_omp_n_thread)

    def __dealloc__(self):
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

        # Pushing the distance and their associated indices on heaps
        # which keep tracks of the argkmin.
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

    cdef void _on_X_prange_iter_init(self,
            ITYPE_t thread_num,
            ITYPE_t X_start,
            ITYPE_t X_end,
    ) nogil:

        # As this strategy is embarrassingly parallel, we can set the
        # thread heaps pointers to the proper position on the main heaps
        self.heaps_approx_distances_chunks[thread_num] = &self.argkmin_distances[X_start, 0]
        self.heaps_indices_chunks[thread_num] = &self.argkmin_indices[X_start, 0]

    cdef void _on_X_prange_iter_finalize(self,
            ITYPE_t thread_num,
            ITYPE_t X_start,
            ITYPE_t X_end,
    ) nogil:
        cdef:
            ITYPE_t idx, jdx

        # Sorting indices of the argkmin for each query vector of X
        for idx in range(X_end - X_start):
            _simultaneous_sort(
                self.heaps_approx_distances_chunks[thread_num] + idx * self.k,
                self.heaps_indices_chunks[thread_num] + idx * self.k,
                self.k
            )

    cdef void _on_Y_init(self,
            ITYPE_t num_threads,
    ) nogil:
        cdef:
            # number of scalar elements
            ITYPE_t heaps_size = self.X_n_samples_chunk * self.k
            ITYPE_t thread_num

        for thread_num in prange(num_threads, schedule='static', nogil=True,
                                 num_threads=num_threads):
            # As chunks of X are shared across threads, so must their
            # heaps. To solve this, each thread has its own heaps
            # which are then synchronised back in the main ones.
            self.heaps_approx_distances_chunks[thread_num] = <DTYPE_t *> malloc(heaps_size * sizeof(DTYPE_t))
            self.heaps_indices_chunks[thread_num] = <ITYPE_t *> malloc(heaps_size * sizeof(ITYPE_t))

    cdef void _on_Y_parallel_init(self,
            ITYPE_t thread_num,
    ) nogil:
        # Initialising heaps (memset can't be used here)
        for idx in range(self.X_n_samples_chunk * self.k):
            self.heaps_approx_distances_chunks[thread_num][idx] = FLOAT_INF
            self.heaps_indices_chunks[thread_num][idx] = -1

    cdef void _on_Y_parallel_finalize(self,
            ITYPE_t thread_num,
            ITYPE_t X_start,
            ITYPE_t X_end,
    ) nogil:
        cdef:
            ITYPE_t idx, jdx
        # TODO: see if synchronisation can be made in parallel samples-wise
        # thanks to the context which wraps the call to this method
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


    cdef void _on_Y_finalize(self,
        ITYPE_t num_threads,
    ) nogil:
        cdef:
            ITYPE_t idx, thread_num

        with nogil, parallel(num_threads=self.effective_omp_n_thread):
            # Deallocating temporary datastructures
            for thread_num in prange(num_threads, schedule='static'):
                free(self.heaps_approx_distances_chunks[thread_num])
                free(self.heaps_indices_chunks[thread_num])

            # Sort the main heaps into arrays in parallel
            # in ascending order w.r.t the distances
            for idx in prange(self.n_X, schedule='static'):
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
        """Convert approximate distances to pairwise distances in parallel."""
        cdef:
            ITYPE_t i, j

        for i in prange(self.n_X, schedule='static', nogil=True,
                        num_threads=self.effective_omp_n_thread):
            for j in range(self.k):
                distances[i, j] = self.distance_metric.dist(&self.X[i, 0],
                                                 &self.Y[Y_indices[i, j], 0],
                                                 self.d)

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

        return self._finalise_compute(return_distance)

    def _finalise_compute(self,
           bint return_distance
    ):
        if return_distance:
            # We need to recompute distances because we relied on
            # approximate distances.
            self._exact_distances(self.argkmin_indices, self.argkmin_distances)
            return np.asarray(self.argkmin_distances), np.asarray(self.argkmin_indices)

        return np.asarray(self.argkmin_indices)

cdef class FastSquaredEuclideanArgKmin(ArgKmin):
    """Fast specialized alternative for ArgKmin on
    EuclideanDistance.

    Computes the argkmin of vectors (rows) of a set of
    vectors (rows) of X on another set of vectors (rows) of Y
    using the GEMM-trick.

    Notes
    -----
    This implementation has an superior arithmetic intensity
    and hence running time, but it can suffer from numerical
    instability. ArgKmin with EuclideanDistance must be
    used when exact precision is needed.
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
                         # The distance metric here is used for exact distances computations
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
            self.Y_n_samples_chunk * self.X_n_samples_chunk * sizeof(DTYPE_t))

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
            self.Y_n_samples_chunk * self.X_n_samples_chunk * sizeof(DTYPE_t))

    cdef void _on_Y_parallel_finalize(self,
            ITYPE_t thread_num,
            ITYPE_t X_start,
            ITYPE_t X_end,
    ) nogil:
        ArgKmin._on_Y_parallel_finalize(self, thread_num, X_start, X_end)
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
            # TODO: necessarily casting because APIs exposed
            # via scipy.linalg.cython_blas aren't reflecting
            # the const-identifier for arguments
            DTYPE_t * A = <DTYPE_t*> & X_c[0, 0]
            ITYPE_t lda = X_c.shape[1]
            DTYPE_t * B = <DTYPE_t*> & Y_c[0, 0]
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
                      # approximate distance: - 2 X_c_i.Y_c_j^T + ||Y_c_j||²
                      dist_middle_terms[i * Y_c.shape[0] + j] + self.Y_sq_norms[j + Y_start],
                      j + Y_start)
        return 0


cdef class RadiusNeighborhood(PairwiseDistancesReduction):

    cdef:
        DTYPE_t radius

        # Distances metrics compute approximated distance
        # ("reduced distance" in the original wording),
        # which are proxies necessitating less computations.
        # We get the proxy for the radius to be able to compare

        # NOTE: not used for now.
        DTYPE_t radius_proxy

        # We want resizable buffers which we will to wrapped within numpy
        # arrays at the end.
        #
        # std::vector comes as a handy interface for efficient resizable
        # buffers.
        #
        # Though it is possible to access their buffer address with
        # std::vector::data, their buffer can't be stolen: their
        # life-time is tight to the buffer's.
        #
        # To solve this, we allocate dynamically allocate vectors which won't be
        # freed, but their buffer eventually will as the ownership will be
        # transferred to numpy arrays.
        #
        # TODO: Find a proper way to handle buffers' lifetime
        # It's "OK-ish" as numpy arrays are then responsible for their buffer
        # lifetime which consist of most of the vectors actual data (residual
        # metadata exist, don't account but won't be deleted).
        #
        # Still, meh.
        vector[vector[ITYPE_t]] * neigh_indices
        vector[vector[DTYPE_t]] * neigh_distances

        vector[vector[ITYPE_t]] ** neigh_indices_chunks
        vector[vector[DTYPE_t]] ** neigh_distances_chunks

        bint sort_results


    @classmethod
    def get_for(cls,
                X,
                Y,
                DTYPE_t radius,
                str metric="euclidean",
                ITYPE_t chunk_size=CHUNK_SIZE,
                dict metric_kwargs=dict(),
        ):
        return RadiusNeighborhood(X=X, Y=Y,
                       distance_metric=DistanceMetric.get_metric(metric, **metric_kwargs),
                       radius=radius,
                       chunk_size=chunk_size)

    def __init__(self,
                 X,
                 Y,
                 DistanceMetric distance_metric,
                 DTYPE_t radius,
                 ITYPE_t chunk_size = CHUNK_SIZE,
    ):
        PairwiseDistancesReduction.__init__(self, X, Y, distance_metric, chunk_size)

        self.radius = radius
        self.sort_results = False

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
            DTYPE_t dist_i_j

        for i in range(X_c.shape[0]):
            for j in range(Y_c.shape[0]):
                dist_i_j = self.distance_metric.dist(
                    &X_c[i, 0], &Y_c[j, 0], self.d)
                if dist_i_j <= self.radius:
                    deref(self.neigh_distances)[X_start + i].push_back(dist_i_j)
                    deref(self.neigh_indices)[X_start + i].push_back(Y_start + j)

        return 0

    cdef void _on_X_prange_iter_init(self,
            ITYPE_t thread_num,
            ITYPE_t X_start,
            ITYPE_t X_end,
    ) nogil:

        # As this strategy is embarrassingly parallel, we can set the
        # thread-local vectors' pointers to the main vectors'.
        self.neigh_distances_chunks[thread_num] = self.neigh_distances
        self.neigh_indices_chunks[thread_num] = self.neigh_indices

    cdef void _on_X_prange_iter_finalize(self,
            ITYPE_t thread_num,
            ITYPE_t X_start,
            ITYPE_t X_end,
    ) nogil:
        cdef:
            ITYPE_t idx, jdx

        # Sorting neighbors for each query vector of X
        if self.sort_results:
            for idx in range(X_start, X_end):
                _simultaneous_sort(
                    deref(self.neigh_distances)[idx].data(),
                    deref(self.neigh_indices)[idx].data(),
                    deref(self.neigh_indices)[idx].size()
                )
    cdef void _on_Y_parallel_init(self,
        ITYPE_t thread_num,
    ) nogil:
        self.neigh_distances_chunks[thread_num] = new vector[vector[DTYPE_t]](self.n_X)
        self.neigh_indices_chunks[thread_num] = new vector[vector[ITYPE_t]](self.n_X)

    cdef void _on_Y_finalize(self,
        ITYPE_t num_threads,
    ) nogil:
        cdef:
            ITYPE_t idx, jdx, thread_num, idx_n_element, idx_current

        # Merge associated vectors into one and sort in parallel
        # in ascending order w.r.t the distances if needed
        for idx in range(self.n_X): #, schedule='static'):
            idx_n_element = 0
            for jdx in range(num_threads):
                idx_n_element += self.neigh_indices_chunks[jdx][idx].size()

            deref(self.neigh_distances)[idx].reserve(idx_n_element)
            deref(self.neigh_indices)[idx].reserve(idx_n_element)

            for jdx in range(num_threads):
                move(deref(self.neigh_distances_chunks[jdx])[idx].begin(),
                     deref(self.neigh_distances_chunks[jdx])[idx].end(),
                     deref(self.neigh_distances)[idx].end())
                move(deref(self.neigh_indices_chunks[jdx])[idx].begin(),
                     deref(self.neigh_indices_chunks[jdx])[idx].end(),
                     deref(self.neigh_indices)[idx].end())

        for jdx in range(num_threads):
            free(self.neigh_distances_chunks[jdx])
            free(self.neigh_indices_chunks[jdx])

        if self.sort_results:
            for idx in prange(self.n_X, schedule='static',
                              num_threads=self.effective_omp_n_thread,
                              nogil=True):
                _simultaneous_sort(
                    deref(self.neigh_distances)[idx].data(),
                    deref(self.neigh_indices)[idx].data(),
                    deref(self.neigh_indices)[idx].size()
                )

        return

    def compute(self,
           str strategy = "auto",
           bint return_distance = False,
           bint sort_results = False
    ):
        # Won't be freed for reasons stated at their definition.
        self.neigh_distances = new vector[vector[DTYPE_t]](self.n_X)
        self.neigh_indices = new vector[vector[ITYPE_t]](self.n_X)

        # TODO: setup thread local datastructures for supporting _parallel_on_Y
        self.sort_results = sort_results
        self._parallel_on_X()
        return self._finalise_compute(return_distance)

    def _finalise_compute(self,
           bint return_distance
    ):
        if return_distance:
            return (_coerce_vectors_to_np_nd_arrays_DTYPE(self.neigh_distances),
                    _coerce_vectors_to_np_nd_arrays_ITYPE(self.neigh_indices))

        # We need to free the buffers here because they won't be managed
        # by a numpy array then.
        free(self.neigh_distances)
        return _coerce_vectors_to_np_nd_arrays_ITYPE(self.neigh_indices)
