# cython: language_level=3
# cython: annotate=False
# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
# cython: profile=False
# cython: linetrace=False
# cython: initializedcheck=False
# cython: binding=False
# distutils: language=c++
# distutils: define_macros=CYTHON_TRACE_NOGIL=0

# Pairwise Distances Reductions
# =============================
#
#    Author: Julien Jerphanion <git@jjerphan.xyz>
#
#
# The routines defined here are used in various algorithms realising
# the same structure of operations on distances between vectors
# of a datasets pair (X, Y).

import numpy as np
cimport numpy as np

np.import_array()

from libc.stdlib cimport free, malloc
from libc.float cimport DBL_MAX
from libcpp.vector cimport vector
from cython cimport final
from cpython.object cimport PyObject
from cython.operator cimport dereference as deref
from cython.parallel cimport parallel, prange
from cpython.ref cimport Py_INCREF

from ._dist_metrics cimport DatasetsPair, DenseDenseDatasetsPair
from ..utils._cython_blas cimport (
  BLAS_Order,
  BLAS_Trans,
  ColMajor,
  NoTrans,
  RowMajor,
  Trans,
  _gemm,
)
from ..utils._heap cimport simultaneous_sort, heap_push
from ..utils._openmp_helpers cimport _openmp_thread_num
from ..utils._typedefs cimport ITYPE_t, DTYPE_t, DITYPE_t
from ..utils._typedefs cimport ITYPECODE, DTYPECODE

from numbers import Integral, Real
from typing import List
from scipy.sparse import issparse
from threadpoolctl import threadpool_limits
from ._dist_metrics import BOOL_METRICS, METRIC_MAPPING
from ..utils import check_array, check_scalar, _in_unstable_openblas_configuration
from ..utils._openmp_helpers import _openmp_effective_n_threads
from ..utils._typedefs import ITYPE, DTYPE

# Those constants have been chosen for modern laptops' caches and architecture.
DEF CHUNK_SIZE = 256  # number of vectors
DEF MIN_CHUNK_SAMPLES = 20


# TODO: change for `libcpp.algorithm.move` once Cython 3 is used
# Introduction in Cython:
# https://github.com/cython/cython/blob/05059e2a9b89bf6738a7750b905057e5b1e3fe2e/Cython/Includes/libcpp/algorithm.pxd#L47 #noqa
cdef extern from "<algorithm>" namespace "std" nogil:
    OutputIt move[InputIt, OutputIt](InputIt first, InputIt last, OutputIt d_first) except + #noqa

######################
## std::vector to np.ndarray coercion
# As type covariance is not supported for C++ containers via Cython,
# we need to redefine fused types.
ctypedef fused vector_DITYPE_t:
    vector[ITYPE_t]
    vector[DTYPE_t]


ctypedef fused vector_vector_DITYPE_t:
    vector[vector[ITYPE_t]]
    vector[vector[DTYPE_t]]


cdef extern from "numpy/arrayobject.h":
    int PyArray_SetBaseObject(np.ndarray arr, PyObject *obj) nogil except -1


cdef class StdVectorSentinel:
    """Wraps a reference to a vector which will be deallocated with this object.

    When created, the StdVectorSentinel swaps the reference of its internal
    vectors with the provided one (vec_ptr), thus making the StdVectorSentinel
    manage the provided one's lifetime.
    """
    pass

# We necessarily need to define two extension types extending StdVectorSentinel
# because we need to provide the dtype of the vector but can't use numeric fused types.
cdef class StdVectorSentinelDTYPE(StdVectorSentinel):
    cdef vector[DTYPE_t] vec

    @staticmethod
    cdef StdVectorSentinel create_for(vector[DTYPE_t] * vec_ptr):
        # This initializes the object directly without calling __init__
        cdef StdVectorSentinelDTYPE sentinel = StdVectorSentinelDTYPE.__new__(StdVectorSentinelDTYPE)
        sentinel.vec.swap(deref(vec_ptr))
        return sentinel


cdef class StdVectorSentinelITYPE(StdVectorSentinel):
    cdef vector[ITYPE_t] vec

    @staticmethod
    cdef StdVectorSentinel create_for(vector[ITYPE_t] * vec_ptr):
        # This initializes the object directly without calling __init__
        cdef StdVectorSentinelITYPE sentinel = StdVectorSentinelITYPE.__new__(StdVectorSentinelITYPE)
        sentinel.vec.swap(deref(vec_ptr))
        return sentinel


cdef np.ndarray vector_to_nd_array(vector_DITYPE_t * vect_ptr):
    """Create a numpy ndarray given a C++ vector.

    The numpy array buffer is the one of the C++ vector.
    A StdVectorSentinel is registers as the base object for the numpy array,
    freeing the C++ vector it encapsulates when the numpy array is freed.
    """
    typenum = DTYPECODE if vector_DITYPE_t is vector[DTYPE_t] else ITYPECODE
    cdef:
        np.npy_intp size = deref(vect_ptr).size()
        # TODO: use PyArray_SimpleNewFromData from the Numpy C API directly
        # I've tried, but Cython fails when parsing the C API
        np.ndarray arr = np.PyArray_SimpleNewFromData(1, &size, typenum,
                                                      deref(vect_ptr).data())
        StdVectorSentinel sentinel

    if vector_DITYPE_t is vector[DTYPE_t]:
        sentinel = StdVectorSentinelDTYPE.create_for(vect_ptr)
    else:
        sentinel = StdVectorSentinelITYPE.create_for(vect_ptr)

    # Makes the numpy array responsible to the life-cycle of its buffer.
    # A reference to the StdVectorSentinel will be stolen by the call bellow,
    # so we increase its reference counter.
    # See: https://docs.python.org/3/c-api/intro.html#reference-count-details
    Py_INCREF(sentinel)
    PyArray_SetBaseObject(arr, <PyObject*>sentinel)
    return arr


cdef np.ndarray[object, ndim=1] coerce_vectors_to_nd_arrays(
    vector_vector_DITYPE_t* vecs
    ):
    """Coerce a std::vector of std::vector to a ndarray of ndarray."""
    cdef:
        ITYPE_t n = deref(vecs).size()
        np.ndarray[object, ndim=1] nd_arrays_of_nd_arrays = np.empty(n,
                                                                     dtype=np.ndarray)

    for i in range(n):
        nd_arrays_of_nd_arrays[i] = vector_to_nd_array(&(deref(vecs)[i]))

    return nd_arrays_of_nd_arrays

#####################

cdef class PairwiseDistancesReduction:
    f"""Abstract class computing a reduction on pairwise
    distances between a set of vectors (rows) X and another
    set of vectors (rows) of Y.

    The implementation of the reduction is done parallelized
    on chunks whose size can be set using ``chunk_size``.

    Parameters
    ----------
    datasets_pair: DatasetsPair
        The pair of dataset to use.

    chunk_size: int, default={CHUNK_SIZE}
        The number of vectors per chunk.

    n_threads: int, default=-1
        The number of OpenMP threads to use for the reduction.
        Parallelism is done on chunks and the sharding of chunks
        depends on the `strategy` set on :method:`~PairwiseDistancesReduction.compute`.

        -1 means using all processors.
    """

    cdef:
        DatasetsPair _datasets_pair

        ITYPE_t n_threads
        ITYPE_t effective_omp_n_thread
        ITYPE_t n_samples_chunk, chunk_size

        ITYPE_t n_X, X_n_samples_chunk, X_n_chunks, X_n_samples_remainder
        ITYPE_t n_Y, Y_n_samples_chunk, Y_n_chunks, Y_n_samples_remainder

    @classmethod
    def valid_metrics(cls) -> List[str]:
        excluded = {
            "pyfunc",  # is relatively slow because we need to coerce data as np arrays
            "mahalanobis", # is numerically unstable
            # TODO: In order to support discrete distance metrics, we need to have a
            # simultaneous sort which breaks ties on indices when distances are identical.
            # The best might be using a std::sort and a Comparator whic might need
            # AoS instead of SoA (currently used).
            "hamming",
            *BOOL_METRICS,
        }
        return sorted({"fast_euclidean", "fast_sqeuclidean", *METRIC_MAPPING.keys()}.difference(excluded))

    @classmethod
    def is_usable_for(cls, X, Y, metric) -> bool:
        # TODO: what's the best coercion for lists?
        X = np.asarray(X) if isinstance(X, (tuple, list)) else X
        Y = np.asarray(Y) if isinstance(Y, (tuple, list)) else Y
        # TODO: support sparse arrays and 32 bits
        return (not issparse(X) and X.dtype.itemsize == 8 and X.ndim == 2 and
                not issparse(Y) and Y.dtype.itemsize == 8 and Y.ndim == 2 and
                metric in cls.valid_metrics())

    @property
    def datasets_pair(self) ->DatasetsPair:
        return self._datasets_pair

    def __init__(self,
        DatasetsPair datasets_pair,
        ITYPE_t chunk_size=CHUNK_SIZE,
        n_threads=-1,
     ):
        cdef:
            ITYPE_t X_n_full_chunks, Y_n_full_chunks

        check_scalar(chunk_size, "chunk_size", Integral, min_val=1)
        self.chunk_size = chunk_size

        if n_threads is None:
            # By convention.
            n_threads = -1

        self.n_threads = n_threads

        if self.n_threads == -1:
            # Using all possible cores
            self.effective_omp_n_thread = _openmp_effective_n_threads()
        else:
            check_scalar(self.n_threads, "n_threads", Integral, min_val=1)
            self.effective_omp_n_thread = self.n_threads

        self.n_samples_chunk = max(MIN_CHUNK_SAMPLES, chunk_size)

        self._datasets_pair = datasets_pair

        self.n_Y = datasets_pair.n_Y()
        self.Y_n_samples_chunk = min(self.n_Y, self.n_samples_chunk)
        Y_n_full_chunks = self.n_Y // self.Y_n_samples_chunk
        self.Y_n_samples_remainder = self.n_Y % self.Y_n_samples_chunk

        self.n_X = datasets_pair.n_X()
        self.X_n_samples_chunk = min(self.n_X, self.n_samples_chunk)
        X_n_full_chunks = self.n_X // self.X_n_samples_chunk
        self.X_n_samples_remainder = self.n_X % self.X_n_samples_chunk

        # Counting remainder chunk in total number of chunks
        self.Y_n_chunks = Y_n_full_chunks + (
            self.n_Y != (Y_n_full_chunks * self.Y_n_samples_chunk)
        )

        self.X_n_chunks = X_n_full_chunks + (
            self.n_X != (X_n_full_chunks * self.X_n_samples_chunk)
        )

    @final
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
            thread_num = _openmp_thread_num()

            # Allocating thread datastructures
            self._on_X_parallel_init(thread_num)

            for X_chunk_idx in prange(self.X_n_chunks, schedule='static'):
                X_start = X_chunk_idx * self.X_n_samples_chunk
                if (X_chunk_idx == self.X_n_chunks - 1
                    and self.X_n_samples_remainder > 0):
                    X_end = X_start + self.X_n_samples_remainder
                else:
                    X_end = X_start + self.X_n_samples_chunk

                # Reinitializing thread datastructures for the new X chunk
                self._on_X_prange_iter_init(thread_num, X_start, X_end)

                for Y_chunk_idx in range(self.Y_n_chunks):
                    Y_start = Y_chunk_idx * self.Y_n_samples_chunk
                    if (Y_chunk_idx == self.Y_n_chunks - 1
                        and self.Y_n_samples_remainder > 0):
                        Y_end = Y_start + self.Y_n_samples_remainder
                    else:
                        Y_end = Y_start + self.Y_n_samples_chunk

                    self._reduce_on_chunks(
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

    @final
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

        # Allocating datastructures
        self._on_Y_init(num_threads)

        for X_chunk_idx in range(self.X_n_chunks):
            X_start = X_chunk_idx * self.X_n_samples_chunk
            if X_chunk_idx == self.X_n_chunks - 1 and self.X_n_samples_remainder > 0:
                X_end = X_start + self.X_n_samples_remainder
            else:
                X_end = X_start + self.X_n_samples_chunk

            with nogil, parallel(num_threads=num_threads):
                thread_num = _openmp_thread_num()

                # Initializing datastructures used in this thread
                self._on_Y_parallel_init(thread_num)

                for Y_chunk_idx in prange(self.Y_n_chunks, schedule='static'):
                    Y_start = Y_chunk_idx * self.Y_n_samples_chunk
                    if Y_chunk_idx == self.Y_n_chunks - 1 \
                            and self.Y_n_samples_remainder > 0:
                        Y_end = Y_start + self.Y_n_samples_remainder
                    else:
                        Y_end = Y_start + self.Y_n_samples_chunk

                    self._reduce_on_chunks(
                        X_start, X_end,
                        Y_start, Y_end,
                        thread_num,
                    )
                # end: prange
            # end: with nogil, parallel

            # Synchronizing the thread datastructures with the main ones
            self._on_Y_after_parallel(num_threads, X_start, X_end)

        # end: for X_chunk_idx
        # Deallocating temporary datastructures and adjusting main datastructures
        self._on_Y_finalize(num_threads)
        return

    # Placeholder methods which have to be implemented

    cdef void _reduce_on_chunks(self,
        ITYPE_t X_start,
        ITYPE_t X_end,
        ITYPE_t Y_start,
        ITYPE_t Y_end,
        ITYPE_t thread_num,
    ) nogil:
        """Implemented the reduction on a pair of chunks."""
        return

    # Placeholder methods which can be implemented

    cdef void compute_exact_distances(self) nogil:
        """Convert proxy distances to exact distances or recompute them."""
        return

    cdef void _on_X_parallel_init(self,
        ITYPE_t thread_num,
    ) nogil:
        """Allocate datastructures used in a thread given its number."""
        return

    cdef void _on_X_prange_iter_init(self,
        ITYPE_t thread_num,
        ITYPE_t X_start,
        ITYPE_t X_end,
    ) nogil:
        """Initialise datastructures used in a thread given its number."""
        return

    cdef void _on_X_prange_iter_finalize(self,
        ITYPE_t thread_num,
        ITYPE_t X_start,
        ITYPE_t X_end,
    ) nogil:
        """Interact with datastructures after a reduction on chunks."""
        return

    cdef void _on_X_parallel_finalize(self,
        ITYPE_t thread_num
    ) nogil:
        """Interact with datastructures after executing all the reductions."""
        return

    cdef void _on_Y_init(self,
        ITYPE_t num_threads,
    ) nogil:
        """Allocate datastructures used in threads."""
        return

    cdef void _on_Y_parallel_init(self,
        ITYPE_t thread_num,
    ) nogil:
        """Initialise datastructures used in a thread given its number."""
        return

    cdef void _on_Y_after_parallel(self,
        ITYPE_t num_threads,
        ITYPE_t X_start,
        ITYPE_t X_end,
    ) nogil:
        """Interact with datastructures after a threads parallel region."""
        return

    cdef void _on_Y_finalize(self,
        ITYPE_t num_threads,
    ) nogil:
        """Interact with datastructures after executing all the reductions."""
        return

cdef class ArgKmin(PairwiseDistancesReduction):
    """Computes the argkmin of vectors (rows) of a set of
    vectors (rows) of X on another set of vectors (rows) of Y.

    The implementation is parallelized on chunks whose size can
    be set using ``chunk_size``.

    Parameters
    ----------
    datasets_pair: DatasetsPair
        The dataset pairs (X, Y) for the reduction.

    k: int
        The k for the argkmin reduction.

    chunk_size: int
        The number of vectors per chunk.

    n_threads: int, default=-1
        The number of OpenMP threads to use for the reduction.
        Parallelism is done on chunks and the sharding of chunks
        depends on the `strategy` set on :method:`~ArgKmin.compute`.

        -1 means using all processors.
    """

    cdef:
        ITYPE_t k

        ITYPE_t[:, ::1] argkmin_indices
        DTYPE_t[:, ::1] argkmin_distances

        # Used as array of pointers to private datastructures used in threads.
        DTYPE_t ** heaps_proxy_distances_chunks
        ITYPE_t ** heaps_indices_chunks

    @classmethod
    def get_for(cls,
        X,
        Y,
        ITYPE_t k,
        str metric="fast_euclidean",
        ITYPE_t chunk_size=CHUNK_SIZE,
        dict metric_kwargs=dict(),
        n_threads=-1,
    ) -> ArgKmin:
        f"""Return the ArgKmin implementation for the given arguments.

        Parameters
        ----------
        X : array-like of shape (n_X, d)
            Input data.

        Y : array-like of shape (n_Y, d)
            Input data.

        k : int
            The k for the argkmin reduction.

        metric : str, default='fast_euclidean'
            The distance metric to use for argkmin. The default metric is
            a fast implementation of the standard Euclidean metric.
            For a list of available metrics, see the documentation of
            :class:`~sklearn.metrics.DistanceMetric`.

        chunk_size: int, default={CHUNK_SIZE},
            The number of vectors per chunk.

        metric_kwargs : dict, default=None
            Keyword arguments to pass to specified metric function.

        n_threads: int, default=-1
            The number of OpenMP threads to use for the reduction.
            Parallelism is done on chunks and the sharding of chunks
            depends on the `strategy` set on :method:`~ArgKmin.compute`.

            -1 means using all processors.

        Returns
        -------
        argkmin: ArgKmin
            The suited ArgKmin implementation.
        """
        # This factory comes to handle specialisations.
        if metric in {"fast_euclidean", "fast_sqeuclidean"} and not issparse(X) and not issparse(Y):
            use_squared_distances = metric == "fast_sqeuclidean"
            return FastEuclideanArgKmin(X=X, Y=Y, k=k,
                                        use_squared_distances=use_squared_distances,
                                        chunk_size=chunk_size)

        return ArgKmin(
            datasets_pair=DatasetsPair.get_for(X, Y, metric, metric_kwargs),
            k=k,
            chunk_size=chunk_size,
        )

    def __init__(self,
        DatasetsPair datasets_pair,
        ITYPE_t k,
        ITYPE_t chunk_size=CHUNK_SIZE,
        n_threads=-1,
    ):
        PairwiseDistancesReduction.__init__(self, datasets_pair, chunk_size, n_threads)

        check_scalar(k, "k", Integral, min_val=1)
        self.k = k

        # Allocating pointers to datastructures but not the datastructures themselves.
        # There as many pointers as available threads.
        # When reducing on small datasets, there can be more pointers than actual
        # threads used for the reduction but there won't be allocated but unused
        # datastructures.
        self.heaps_proxy_distances_chunks = <DTYPE_t **> malloc(
            sizeof(DTYPE_t *) * self.effective_omp_n_thread)
        self.heaps_indices_chunks = <ITYPE_t **> malloc(
            sizeof(ITYPE_t *) * self.effective_omp_n_thread)

    def __dealloc__(self):
        if self.heaps_indices_chunks is not NULL:
            free(self.heaps_indices_chunks)

        if self.heaps_proxy_distances_chunks is not NULL:
            free(self.heaps_proxy_distances_chunks)

    cdef void _reduce_on_chunks(self,
        ITYPE_t X_start,
        ITYPE_t X_end,
        ITYPE_t Y_start,
        ITYPE_t Y_end,
        ITYPE_t thread_num,
    ) nogil:
        cdef:
            ITYPE_t i, j
            ITYPE_t n_X = X_end - X_start
            ITYPE_t n_Y = Y_end - Y_start
            ITYPE_t k = self.k
            DTYPE_t *heaps_proxy_distances = self.heaps_proxy_distances_chunks[thread_num]
            ITYPE_t *heaps_indices = self.heaps_indices_chunks[thread_num]

        # Pushing the distance and their associated indices on heaps
        # which keep tracks of the argkmin.
        for i in range(n_X):
            for j in range(n_Y):
                heap_push(
                    heaps_proxy_distances + i * self.k,
                    heaps_indices + i * self.k,
                    k,
                    self._datasets_pair.proxy_dist(X_start + i, Y_start + j),
                    Y_start + j,
                )

    @final
    cdef void _on_X_prange_iter_init(self,
        ITYPE_t thread_num,
        ITYPE_t X_start,
        ITYPE_t X_end,
    ) nogil:
        # As this strategy is embarrassingly parallel, we can set the
        # thread heaps pointers to the proper position on the main heaps
        self.heaps_proxy_distances_chunks[thread_num] = &self.argkmin_distances[X_start, 0]
        self.heaps_indices_chunks[thread_num] = &self.argkmin_indices[X_start, 0]

    @final
    cdef void _on_X_prange_iter_finalize(self,
        ITYPE_t thread_num,
        ITYPE_t X_start,
        ITYPE_t X_end,
    ) nogil:
        cdef:
            ITYPE_t idx, jdx

        # Sorting indices of the argkmin for each query vector of X
        for idx in range(X_end - X_start):
            simultaneous_sort(
                self.heaps_proxy_distances_chunks[thread_num] + idx * self.k,
                self.heaps_indices_chunks[thread_num] + idx * self.k,
                self.k
            )

    cdef void _on_Y_init(self,
        ITYPE_t num_threads,
    ) nogil:
        cdef:
            # Maximum number of scalar elements (the last chunks can be smaller)
            ITYPE_t heaps_size = self.X_n_samples_chunk * self.k
            ITYPE_t thread_num

        for thread_num in prange(num_threads, schedule='static', nogil=True,
                                 num_threads=num_threads):
            # As chunks of X are shared across threads, so must their
            # heaps. To solve this, each thread has its own heaps
            # which are then synchronised back in the main ones.
            self.heaps_proxy_distances_chunks[thread_num] = <DTYPE_t *> malloc(
                heaps_size * sizeof(DTYPE_t))
            self.heaps_indices_chunks[thread_num] = <ITYPE_t *> malloc(
                heaps_size * sizeof(ITYPE_t))

    @final
    cdef void _on_Y_parallel_init(self,
        ITYPE_t thread_num,
    ) nogil:
        # Initialising heaps (memset can't be used here)
        for idx in range(self.X_n_samples_chunk * self.k):
            self.heaps_proxy_distances_chunks[thread_num][idx] = DBL_MAX
            self.heaps_indices_chunks[thread_num][idx] = -1

    @final
    cdef void _on_Y_after_parallel(self,
        ITYPE_t num_threads,
        ITYPE_t X_start,
        ITYPE_t X_end,
    ) nogil:
        cdef:
            ITYPE_t idx, jdx, thread_num
        with nogil, parallel(num_threads=self.effective_omp_n_thread):
            # Synchronising the thread heaps with the main heaps
            # This is done in parallel samples-wise (no need for locks)
            for idx in prange(X_end - X_start, schedule="static"):
                for thread_num in range(num_threads):
                    for jdx in range(self.k):
                        heap_push(
                            &self.argkmin_distances[X_start + idx, 0],
                            &self.argkmin_indices[X_start + idx, 0],
                            self.k,
                            self.heaps_proxy_distances_chunks[thread_num][idx * self.k + jdx],
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
                free(self.heaps_proxy_distances_chunks[thread_num])
                free(self.heaps_indices_chunks[thread_num])

            # Sort the main heaps into arrays in parallel
            # in ascending order w.r.t the distances
            for idx in prange(self.n_X, schedule='static'):
                simultaneous_sort(
                    &self.argkmin_distances[idx, 0],
                    &self.argkmin_indices[idx, 0],
                    self.k,
                )
        return

    cdef void compute_exact_distances(self) nogil:
        cdef:
            ITYPE_t i, j
            ITYPE_t[:, ::1] Y_indices = self.argkmin_indices
            DTYPE_t[:, ::1] distances = self.argkmin_distances
        for i in prange(self.n_X, schedule='static', nogil=True,
                        num_threads=self.effective_omp_n_thread):
            for j in range(self.k):
                distances[i, j] = self._datasets_pair.distance_metric._rdist_to_dist(
                    # Guard against eventual -0., causing nan production.
                    distances[i, j] if distances[i, j] > 0. else 0.
                )

    @final
    def compute(self,
       str strategy = "auto",
       bint return_distance = False
    ):
        """Computes the reduction of vectors (rows) of X on Y.

        Parameters
        ----------
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

        # Results returned by ArgKmin.compute used as the main heaps
        self.argkmin_indices = np.full((self.n_X, self.k), 0, dtype=ITYPE)
        self.argkmin_distances = np.full((self.n_X, self.k), DBL_MAX, dtype=DTYPE)

        if strategy == 'auto':
            # This is a simple heuristic whose constant for the
            # comparison has been chosen based on experiments.
            if 4 * self.chunk_size * self.effective_omp_n_thread < self.n_X:
                strategy = 'parallel_on_X'
            else:
                strategy = 'parallel_on_Y'

        # Limit the number of threads in second level of nested parallelism for BLAS
        # to avoid threads over-subscription (in GEMM for instance).
        with threadpool_limits(limits=1, user_api="blas"):
            if strategy == 'parallel_on_Y':
                self._parallel_on_Y()
            elif strategy == 'parallel_on_X':
                self._parallel_on_X()
            else:
                raise RuntimeError(f"strategy '{strategy}' not supported.")

        if return_distance:
            # We eventually need to recompute distances because we relied on proxies.
            self.compute_exact_distances()
            return np.asarray(self.argkmin_distances), np.asarray(self.argkmin_indices)

        return np.asarray(self.argkmin_indices)


cdef class FastEuclideanArgKmin(ArgKmin):
    """Fast specialized alternative for ArgKmin on EuclideanDistance.

    Notes
    -----
    This implementation has a superior arithmetic intensity and hence
    better running time when the alternative is IO bound, but it can suffer
    from numerical instability.

    ArgKmin with EuclideanDistance must be used when higher numerical precision
    is needed.
    """

    cdef:
        const DTYPE_t[:, ::1] X
        const DTYPE_t[:, ::1] Y
        const DTYPE_t[::1] X_sq_norms
        const DTYPE_t[::1] Y_sq_norms

        # Buffers for GEMM
        DTYPE_t ** dist_middle_terms_chunks
        bint use_squared_distances

    @classmethod
    def is_usable_for(cls, X, Y, metric) -> bool:
        return (ArgKmin.is_usable_for(X, Y, metric) and
                not _in_unstable_openblas_configuration())

    def __init__(self,
        X,
        Y,
        ITYPE_t k,
        bint use_squared_distances=False,
        ITYPE_t chunk_size=CHUNK_SIZE,
    ):
        ArgKmin.__init__(self,
            # The datasets pair here is used for exact distances computations
            datasets_pair=DatasetsPair.get_for(X, Y, metric="euclidean"),
            k=k,
            chunk_size=chunk_size)
        # X and Y are checked by the DatasetsPair implemented as a DenseDenseDatasetsPair
        cdef:
            DenseDenseDatasetsPair datasets_pair = <DenseDenseDatasetsPair> self.datasets_pair
        self.X, self.Y = datasets_pair.X, datasets_pair.Y
        self.X_sq_norms = np.einsum('ij,ij->i', self.X, self.X)
        self.Y_sq_norms = np.einsum('ij,ij->i', self.Y, self.Y)
        self.use_squared_distances = use_squared_distances

        # Temporary datastructures used in threads
        self.dist_middle_terms_chunks = <DTYPE_t **> malloc(
            sizeof(DTYPE_t *) * self.effective_omp_n_thread)

    def __dealloc__(self):
        if self.dist_middle_terms_chunks is not NULL:
            free(self.dist_middle_terms_chunks)

    @final
    cdef void compute_exact_distances(self) nogil:
        if not self.use_squared_distances:
            ArgKmin.compute_exact_distances(self)

    @final
    cdef void _on_X_parallel_init(self,
            ITYPE_t thread_num,
    ) nogil:
        ArgKmin._on_X_parallel_init(self, thread_num)

        # Temporary buffer for the -2 * X_c.dot(Y_c.T) term
        self.dist_middle_terms_chunks[thread_num] = <DTYPE_t *> malloc(
            self.Y_n_samples_chunk * self.X_n_samples_chunk * sizeof(DTYPE_t))

    @final
    cdef void _on_X_parallel_finalize(self,
            ITYPE_t thread_num
    ) nogil:
        ArgKmin._on_X_parallel_finalize(self, thread_num)
        free(self.dist_middle_terms_chunks[thread_num])

    @final
    cdef void _on_Y_init(self,
            ITYPE_t num_threads,
    ) nogil:
        cdef ITYPE_t thread_num
        ArgKmin._on_Y_init(self, num_threads)

        for thread_num in range(num_threads):
            # Temporary buffer for the -2 * X_c.dot(Y_c.T) term
            self.dist_middle_terms_chunks[thread_num] = <DTYPE_t *> malloc(
                self.Y_n_samples_chunk * self.X_n_samples_chunk * sizeof(DTYPE_t))

    @final
    cdef void _on_Y_finalize(self,
            ITYPE_t num_threads,
    ) nogil:
        cdef ITYPE_t thread_num
        ArgKmin._on_Y_finalize(self, num_threads)

        for thread_num in range(num_threads):
            free(self.dist_middle_terms_chunks[thread_num])

    @final
    cdef void _reduce_on_chunks(self,
        ITYPE_t X_start,
        ITYPE_t X_end,
        ITYPE_t Y_start,
        ITYPE_t Y_end,
        ITYPE_t thread_num,
    ) nogil:
        cdef:
            ITYPE_t i, j
            ITYPE_t k = self.k

            const DTYPE_t[:, ::1] X_c = self.X[X_start:X_end, :]
            const DTYPE_t[:, ::1] Y_c = self.Y[Y_start:Y_end, :]
            DTYPE_t *dist_middle_terms = self.dist_middle_terms_chunks[thread_num]
            DTYPE_t *heaps_proxy_distances = self.heaps_proxy_distances_chunks[thread_num]
            ITYPE_t *heaps_indices = self.heaps_indices_chunks[thread_num]

            # We compute the full pairwise squared distances matrix as follows
            #
            #      ||X_c - Y_c||² = ||X_c||² - 2 X_c.Y_c^T + ||Y_c||²,
            #
            # The middle term gets computed efficiently bellow using BLAS Level 3 GEMM.
            #
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
                heap_push(
                    heaps_proxy_distances + i * k,
                    heaps_indices + i * k,
                    k,
                    # proxy distance: ||X_c_i||² - 2 X_c_i.Y_c_j^T + ||Y_c_j||²
                    (
                        self.X_sq_norms[i + X_start] +
                        dist_middle_terms[i * Y_c.shape[0] + j] +
                        self.Y_sq_norms[j + Y_start]
                    ),
                    j + Y_start,
                )


cdef class RadiusNeighborhood(PairwiseDistancesReduction):
    """Returns radius-based neighbors vectors' indices in a dataset Y of
    of vectors in a dataset X.

    Parameters
    ----------
    datasets_pair: DatasetsPair
        The dataset pairs (X, Y) for the reduction.

    radius: float
        The radius defining the neighborhood.

    chunk_size: int
        The number of vectors per chunk.

    n_threads: int, default=-1
        The number of OpenMP threads to use for the reduction.
        Parallelism is done on chunks and the sharding of chunks
        depends on the `strategy` set on :method:`~RadiusNeighborhood.compute`.

        -1 means using all processors.
    """

    cdef:
        DTYPE_t radius

        # DistanceMetric compute rank preserving distance via rdist
        # ("reduced distance" in the original wording),
        # which are proxies necessitating less computations.
        # We get the proxy for the radius to be able to compare
        DTYPE_t proxy_radius

        # Neighbors informations are returned as np.ndarray or np.ndarray.
        #
        # We want resizable buffers which we will to wrapped within numpy
        # arrays at the end. std::vector comes as a handy interface for
        # interacting efficiently with resizable buffers.
        #
        # Though it is possible to access their buffer address with
        # std::vector::data, they can't be stolen: buffers lifetime
        # is tight to their std::vector and are deallocated when
        # std::vectors are.
        #
        # To solve this, we dynamically allocate std::vectors and then
        # encapsulate them in a StdVectorSentinel responsible for
        # freeing them when the associated np.ndarray is freed.
        vector[vector[ITYPE_t]] * neigh_indices
        vector[vector[DTYPE_t]] * neigh_distances

        # Used as array of pointers to private datastructures used in threads.
        vector[vector[ITYPE_t]] ** neigh_indices_chunks
        vector[vector[DTYPE_t]] ** neigh_distances_chunks

        bint sort_results

    @classmethod
    def get_for(cls,
        X,
        Y,
        DTYPE_t radius,
        str metric="fast_euclidean",
        ITYPE_t chunk_size=CHUNK_SIZE,
        dict metric_kwargs=dict(),
        n_threads=-1,
    ) -> RadiusNeighborhood:
        f"""Return the RadiusNeighborhood implementation for the given arguments.

        Parameters
        ----------
        X : array-like of shape (n_X, d)
            Input data.

        Y : array-like of shape (n_Y, d)
            Input data.

        radius : float
            The radius defining the neighborhood.

        metric : str, default='fast_euclidean'
            The distance metric to use for argkmin. The default metric is
            a fast implementation of the standard Euclidean metric.
            For a list of available metrics, see the documentation of
            :class:`~sklearn.metrics.DistanceMetric`.

        chunk_size: int, default={CHUNK_SIZE},
            The number of vectors per chunk.

        metric_kwargs : dict, default=None
            Keyword arguments to pass to specified metric function.

        n_threads: int, default=-1
            The number of OpenMP threads to use for the reduction.
            Parallelism is done on chunks and the sharding of chunks
            depends on the `strategy` set on :method:`~RadiusNeighborhood.compute`.

            -1 means using all processors.

        Returns
        -------
        radius_neighborhood: RadiusNeighborhood
            The suited RadiusNeighborhood implementation.
        """
        # This factory comes to handle specialisations.
        if metric in {"fast_euclidean", "fast_sqeuclidean"} and not issparse(X) and not issparse(Y):
            use_squared_distances = metric == "fast_sqeuclidean"
            return FastEuclideanRadiusNeighborhood(X=X, Y=Y,
                                                   radius=radius,
                                                   use_squared_distances=use_squared_distances,
                                                   chunk_size=chunk_size)

        return RadiusNeighborhood(
            datasets_pair=DatasetsPair.get_for(X, Y, metric, metric_kwargs),
            radius=radius,
            chunk_size=chunk_size,
        )

    def __init__(self,
        DatasetsPair datasets_pair,
        DTYPE_t radius,
        ITYPE_t chunk_size=CHUNK_SIZE,
        n_threads=-1,
    ):
        PairwiseDistancesReduction.__init__(self, datasets_pair, chunk_size, n_threads)

        check_scalar(radius, "radius", Real, min_val=0)
        self.radius = radius
        self.proxy_radius = self._datasets_pair.distance_metric._dist_to_rdist(radius)
        self.sort_results = False

        # Allocating pointers to datastructures but not the datastructures themselves.
        # There as many pointers as available threads.
        # When reducing on small datasets, there can be more pointers than actual
        # threads used for the reduction but there won't be allocated but unused
        # datastructures.
        self.neigh_distances_chunks = <vector[vector[DTYPE_t]] **> malloc(
            sizeof(self.neigh_distances) * self.effective_omp_n_thread)
        self.neigh_indices_chunks = <vector[vector[ITYPE_t]] **> malloc(
            sizeof(self.neigh_indices) * self.effective_omp_n_thread)

    def __dealloc__(self):
        if self.neigh_distances_chunks is not NULL:
            free(self.neigh_distances_chunks)

        if self.neigh_indices_chunks is not NULL:
            free(self.neigh_indices_chunks)

    cdef void _reduce_on_chunks(self,
        ITYPE_t X_start,
        ITYPE_t X_end,
        ITYPE_t Y_start,
        ITYPE_t Y_end,
        ITYPE_t thread_num,
    ) nogil:
        cdef:
            ITYPE_t i, j
            DTYPE_t proxy_dist_i_j

        for i in range(X_start, X_end):
            for j in range(Y_start, Y_end):
                proxy_dist_i_j = self._datasets_pair.proxy_dist(i, j)
                if proxy_dist_i_j <= self.proxy_radius:
                    deref(self.neigh_distances_chunks[thread_num])[i].push_back(proxy_dist_i_j)
                    deref(self.neigh_indices_chunks[thread_num])[i].push_back(j)

    @final
    cdef void _on_X_prange_iter_init(self,
        ITYPE_t thread_num,
        ITYPE_t X_start,
        ITYPE_t X_end,
    ) nogil:

        # As this strategy is embarrassingly parallel, we can set the
        # thread vectors' pointers to the main vectors'.
        self.neigh_distances_chunks[thread_num] = self.neigh_distances
        self.neigh_indices_chunks[thread_num] = self.neigh_indices

    @final
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
                simultaneous_sort(
                    deref(self.neigh_distances)[idx].data(),
                    deref(self.neigh_indices)[idx].data(),
                    deref(self.neigh_indices)[idx].size()
                )

    cdef void _on_Y_init(self,
        ITYPE_t num_threads,
    ) nogil:
        cdef:
            ITYPE_t thread_num
        # As chunks of X are shared across threads, so must datastructures
        # to avoid race conditions.
        # Each thread has its own vectors of n_X vectors which are then merged
        # back in the main n_X vectors.
        for thread_num in range(num_threads):
            self.neigh_distances_chunks[thread_num] = new vector[vector[DTYPE_t]](self.n_X)
            self.neigh_indices_chunks[thread_num] = new vector[vector[ITYPE_t]](self.n_X)

    @final
    cdef void _merge_vectors(self,
        ITYPE_t idx,
        ITYPE_t num_threads,
    ) nogil:
        cdef:
            ITYPE_t thread_num
            ITYPE_t idx_n_elements = 0
            ITYPE_t last_element_idx = deref(self.neigh_indices)[idx].size()

        # Resizing buffers only once for the given
        for thread_num in range(num_threads):
            idx_n_elements += deref(self.neigh_distances_chunks[thread_num])[idx].size()

        deref(self.neigh_distances)[idx].resize(last_element_idx + idx_n_elements)
        deref(self.neigh_indices)[idx].resize(last_element_idx + idx_n_elements)

        # Moving the elements by range using the range first element
        # as the reference for the insertion
        for thread_num in range(num_threads):
            move(
                deref(self.neigh_distances_chunks[thread_num])[idx].begin(),
                deref(self.neigh_distances_chunks[thread_num])[idx].end(),
                deref(self.neigh_distances)[idx].begin() + last_element_idx
            )
            move(
                deref(self.neigh_indices_chunks[thread_num])[idx].begin(),
                deref(self.neigh_indices_chunks[thread_num])[idx].end(),
                deref(self.neigh_indices)[idx].begin() + last_element_idx
            )
            last_element_idx += deref(self.neigh_distances_chunks[thread_num])[idx].size()


    cdef void _on_Y_finalize(self,
        ITYPE_t num_threads,
    ) nogil:
        cdef:
            ITYPE_t idx, jdx, thread_num, idx_n_element, idx_current

        with nogil, parallel(num_threads=self.effective_omp_n_thread):
            # Merge vectors used in threads into the main ones.
            # This is done in parallel sample-wise (no need for locks)
            # using dynamic scheduling because we generally do not have
            # the same number of neighbors for each query vectors.
            # TODO: compare 'dynamic' vs 'static' vs 'guided'
            for idx in prange(self.n_X, schedule='dynamic'):
                self._merge_vectors(idx, num_threads)

            # The content of the vector have been std::moved,
            # Hence they can't be used anymore and can only be deleted.
            for thread_num in prange(num_threads, schedule='static'):
                del self.neigh_distances_chunks[thread_num]
                del self.neigh_indices_chunks[thread_num]

            # Sort in parallel in ascending order w.r.t the distances if needed
            if self.sort_results:
                for idx in prange(self.n_X, schedule='static'):
                    simultaneous_sort(
                        deref(self.neigh_distances)[idx].data(),
                        deref(self.neigh_indices)[idx].data(),
                        deref(self.neigh_indices)[idx].size()
                    )

        return

    cdef void compute_exact_distances(self) nogil:
        """Convert proxy distances to pairwise distances in parallel."""
        cdef:
            ITYPE_t i, j

        for i in prange(self.n_X, nogil=True, schedule='static',
                        num_threads=self.effective_omp_n_thread):
            for j in range(deref(self.neigh_indices)[i].size()):
                deref(self.neigh_distances)[i][j] = (
                        self._datasets_pair.distance_metric._rdist_to_dist(
                            # Guard against eventual -0., causing nan production.
                            deref(self.neigh_distances)[i][j]
                            if deref(self.neigh_distances)[i][j] > 0.
                            else 0
                        )
                )

    @final
    def compute(self,
        str strategy = "auto",
        bint return_distance = False,
        bint sort_results = False
    ):
        if sort_results and not return_distance:
            raise ValueError("return_distance must be True "
                             "if sort_results is True.")

        # Temporary datastructures which will be coerced to
        # numpy arrays on return and then freed.
        self.neigh_indices = new vector[vector[ITYPE_t]](self.n_X)
        self.neigh_distances = new vector[vector[DTYPE_t]](self.n_X)

        self.sort_results = sort_results

        if strategy == 'auto':
            # This is a simple heuristic whose constant for the
            # comparison has been chosen based on experiments.
            if 4 * self.chunk_size * self.effective_omp_n_thread < self.n_X:
                strategy = 'parallel_on_X'
            else:
                strategy = 'parallel_on_Y'

        # Limit the number of threads in second level of nested parallelism for BLAS
        # to avoid threads over-subscription (in GEMM for instance).
        with threadpool_limits(limits=1, user_api="blas"):
            if strategy == 'parallel_on_Y':
                self._parallel_on_Y()
            elif strategy == 'parallel_on_X':
                self._parallel_on_X()
            else:
                raise RuntimeError(f"strategy '{strategy}' not supported.")

        if return_distance:
            self.compute_exact_distances()
            res = (
                coerce_vectors_to_nd_arrays(self.neigh_distances),
                coerce_vectors_to_nd_arrays(self.neigh_indices),
            )
        else:
            res = coerce_vectors_to_nd_arrays(self.neigh_indices)

        del self.neigh_distances
        del self.neigh_indices

        return res


cdef class FastEuclideanRadiusNeighborhood(RadiusNeighborhood):
    """Fast specialized alternative for RadiusNeighborhood on EuclideanDistance.

    Notes
    -----
    This implementation has a superior arithmetic intensity and hence
    better running time when the alternative is IO bound, but it can suffer
    from numerical instability.

    RadiusNeighborhood with EuclideanDistance must be used when higher
    numerical precision is needed.
    """

    cdef:
        const DTYPE_t[:, ::1] X
        const DTYPE_t[:, ::1] Y
        const DTYPE_t[::1] X_sq_norms
        const DTYPE_t[::1] Y_sq_norms

        # Buffers for GEMM
        DTYPE_t ** dist_middle_terms_chunks
        bint use_squared_distances

    @classmethod
    def is_usable_for(cls, X, Y, metric) -> bool:
        return (RadiusNeighborhood.is_usable_for(X, Y, metric)
                and not _in_unstable_openblas_configuration())

    def __init__(self,
        X,
        Y,
        DTYPE_t radius,
        bint use_squared_distances=False,
        ITYPE_t chunk_size=CHUNK_SIZE,
    ):
        RadiusNeighborhood.__init__(self,
            # The datasets pair here is used for exact distances computations
            datasets_pair=DatasetsPair.get_for(X, Y, metric="euclidean"),
            radius=radius,
            chunk_size=chunk_size)
        # X and Y are checked by the DatasetsPair implemented as a DenseDenseDatasetsPair
        cdef:
            DenseDenseDatasetsPair datasets_pair = <DenseDenseDatasetsPair> self.datasets_pair
        self.X, self.Y = datasets_pair.X, datasets_pair.Y
        self.X_sq_norms = np.einsum('ij,ij->i', self.X, self.X)
        self.Y_sq_norms = np.einsum('ij,ij->i', self.Y, self.Y)
        self.use_squared_distances = use_squared_distances

        if use_squared_distances:
            # In this specialisation and this setup, the value passed to the radius is
            # already considered to be the the proxy radius, so we overwrite it.
            self.proxy_radius = radius

        # Temporary datastructures used in threads
        self.dist_middle_terms_chunks = <DTYPE_t **> malloc(
            sizeof(DTYPE_t *) * self.effective_omp_n_thread)

    def __dealloc__(self):
        if self.dist_middle_terms_chunks is not NULL:
            free(self.dist_middle_terms_chunks)

    @final
    cdef void compute_exact_distances(self) nogil:
        if not self.use_squared_distances:
            RadiusNeighborhood.compute_exact_distances(self)

    @final
    cdef void _on_X_parallel_init(self,
        ITYPE_t thread_num,
    ) nogil:
        RadiusNeighborhood._on_X_parallel_init(self, thread_num)

        # Temporary buffer for the -2 * X_c.dot(Y_c.T) term
        self.dist_middle_terms_chunks[thread_num] = <DTYPE_t *> malloc(
            self.Y_n_samples_chunk * self.X_n_samples_chunk * sizeof(DTYPE_t))

    @final
    cdef void _on_X_parallel_finalize(self,
        ITYPE_t thread_num
    ) nogil:
        RadiusNeighborhood._on_X_parallel_finalize(self, thread_num)
        free(self.dist_middle_terms_chunks[thread_num])

    @final
    cdef void _on_Y_init(self,
        ITYPE_t num_threads,
    ) nogil:
        cdef ITYPE_t thread_num
        RadiusNeighborhood._on_Y_init(self, num_threads)

        for thread_num in range(num_threads):
            # Temporary buffer for the -2 * X_c.dot(Y_c.T) term
            self.dist_middle_terms_chunks[thread_num] = <DTYPE_t *> malloc(
                self.Y_n_samples_chunk * self.X_n_samples_chunk * sizeof(DTYPE_t))

    @final
    cdef void _on_Y_finalize(self,
        ITYPE_t num_threads,
    ) nogil:
        cdef ITYPE_t thread_num
        RadiusNeighborhood._on_Y_finalize(self, num_threads)

        for thread_num in range(num_threads):
            free(self.dist_middle_terms_chunks[thread_num])

    @final
    cdef void _reduce_on_chunks(self,
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
            DTYPE_t *dist_middle_terms = self.dist_middle_terms_chunks[thread_num]

            # We compute the full pairwise squared distances matrix as follows
            #
            #      ||X_c - Y_c||² = ||X_c||² - 2 X_c.Y_c^T + ||Y_c||²,
            #
            # The middle term gets computed efficiently bellow using BLAS Level 3 GEMM.
            #
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
                # ||X_c_i||² - 2 X_c_i.Y_c_j^T + ||Y_c_j||²
                squared_dist_i_j = (self.X_sq_norms[i + X_start]
                            + dist_middle_terms[i * Y_c.shape[0] + j]
                            + self.Y_sq_norms[j + Y_start])
                if squared_dist_i_j <= self.proxy_radius:
                    deref(self.neigh_distances_chunks[thread_num])[i + X_start].push_back(squared_dist_i_j)
                    deref(self.neigh_indices_chunks[thread_num])[i + X_start].push_back(j + Y_start)
