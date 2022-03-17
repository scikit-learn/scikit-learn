# Pairwise Distances Reductions
# =============================
#
#    Author: Julien Jerphanion <git@jjerphan.xyz>
#
#
# The abstractions defined here are used in various algorithms performing
# the same structure of operations on distances between row vectors
# of a datasets pair (X, Y).
#
# Importantly, the core of the computation is chunked to make sure that the pairwise
# distance chunk matrices stay in CPU cache before applying the final reduction step.
# Furthermore, the chunking strategy is also used to leverage OpenMP-based parallelism
# (using Cython prange loops) which gives another multiplicative speed-up in
# favorable cases on many-core machines.
cimport numpy as np
import numpy as np
import warnings

from .. import get_config
from libc.stdlib cimport free, malloc
from libc.float cimport DBL_MAX
from libcpp.memory cimport shared_ptr, make_shared
from libcpp.vector cimport vector
from cython cimport final
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
  _dot,
  _gemm,
)
from ..utils._heap cimport simultaneous_sort, heap_push
from ..utils._openmp_helpers cimport _openmp_thread_num
from ..utils._typedefs cimport ITYPE_t, DTYPE_t
from ..utils._typedefs cimport ITYPECODE, DTYPECODE

from numbers import Integral, Real
from typing import List
from scipy.sparse import issparse
from ._dist_metrics import BOOL_METRICS, METRIC_MAPPING
from ..utils import check_scalar, _in_unstable_openblas_configuration
from ..utils.fixes import threadpool_limits
from ..utils._openmp_helpers import _openmp_effective_n_threads
from ..utils._typedefs import ITYPE, DTYPE


np.import_array()

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
        # See: https://cython.readthedocs.io/en/latest/src/userguide/extension_types.html#instantiation-from-existing-c-c-pointers # noqa
        cdef StdVectorSentinelDTYPE sentinel = StdVectorSentinelDTYPE.__new__(StdVectorSentinelDTYPE)
        sentinel.vec.swap(deref(vec_ptr))
        return sentinel


cdef class StdVectorSentinelITYPE(StdVectorSentinel):
    cdef vector[ITYPE_t] vec

    @staticmethod
    cdef StdVectorSentinel create_for(vector[ITYPE_t] * vec_ptr):
        # This initializes the object directly without calling __init__
        # See: https://cython.readthedocs.io/en/latest/src/userguide/extension_types.html#instantiation-from-existing-c-c-pointers # noqa
        cdef StdVectorSentinelITYPE sentinel = StdVectorSentinelITYPE.__new__(StdVectorSentinelITYPE)
        sentinel.vec.swap(deref(vec_ptr))
        return sentinel


cdef np.ndarray vector_to_nd_array(vector_DITYPE_t * vect_ptr):
    """Create a numpy ndarray given a C++ vector.

    The numpy array buffer is the one of the C++ vector.
    A StdVectorSentinel is registered as the base object for the numpy array,
    freeing the C++ vector it encapsulates when the numpy array is freed.
    """
    typenum = DTYPECODE if vector_DITYPE_t is vector[DTYPE_t] else ITYPECODE
    cdef:
        np.npy_intp size = deref(vect_ptr).size()
        np.ndarray arr = np.PyArray_SimpleNewFromData(1, &size, typenum,
                                                      deref(vect_ptr).data())
        StdVectorSentinel sentinel

    if vector_DITYPE_t is vector[DTYPE_t]:
        sentinel = StdVectorSentinelDTYPE.create_for(vect_ptr)
    else:
        sentinel = StdVectorSentinelITYPE.create_for(vect_ptr)

    # Makes the numpy array responsible of the life-cycle of its buffer.
    # A reference to the StdVectorSentinel will be stolen by the call to
    # `PyArray_SetBaseObject` below, so we increase its reference counter.
    # See: https://docs.python.org/3/c-api/intro.html#reference-count-details
    Py_INCREF(sentinel)
    np.PyArray_SetBaseObject(arr, sentinel)
    return arr


cdef np.ndarray[object, ndim=1] coerce_vectors_to_nd_arrays(
    shared_ptr[vector_vector_DITYPE_t] vecs
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

cpdef DTYPE_t[::1] _sqeuclidean_row_norms(
    const DTYPE_t[:, ::1] X,
    ITYPE_t num_threads,
):
    """Compute the squared euclidean norm of the rows of X in parallel.

    This is faster than using np.einsum("ij, ij->i") even when using a single thread.
    """
    cdef:
        # Casting for X to remove the const qualifier is needed because APIs
        # exposed via scipy.linalg.cython_blas aren't reflecting the arguments'
        # const qualifier.
        # See: https://github.com/scipy/scipy/issues/14262
        DTYPE_t * X_ptr = <DTYPE_t *> &X[0, 0]
        ITYPE_t idx = 0
        ITYPE_t n = X.shape[0]
        ITYPE_t d = X.shape[1]
        DTYPE_t[::1] squared_row_norms = np.empty(n, dtype=DTYPE)

    for idx in prange(n, schedule='static', nogil=True, num_threads=num_threads):
        squared_row_norms[idx] = _dot(d, X_ptr + idx * d, 1, X_ptr + idx * d, 1)

    return squared_row_norms

#####################

cdef class PairwiseDistancesReduction:
    """Abstract base class for pairwise distance computation & reduction.

    Subclasses of this class compute pairwise distances between a set of
    row vectors of X and another set of row vectors of Y and apply a reduction on top.
    The reduction takes a matrix of pairwise distances between rows of X and Y
    as input and outputs an aggregate data-structure for each row of X.
    The aggregate values are typically smaller than the number of rows in Y,
    hence the term reduction.

    For computational reasons, it is interesting to perform the reduction on
    the fly on chunks of rows of X and Y so as to keep intermediate
    data-structures in CPU cache and avoid unnecessary round trips of large
    distance arrays with the RAM that would otherwise severely degrade the
    speed by making the overall processing memory-bound.

    The base class provides the generic chunked parallelization template using
    OpenMP loops (Cython prange), either on rows of X or rows of Y depending on
    their respective sizes.

    The subclasses are specialized for reduction.

    The actual distance computation for a given pair of rows of X and Y are
    delegated to format-specific subclasses of the DatasetsPair companion base
    class.

    Parameters
    ----------
    datasets_pair: DatasetsPair
        The pair of dataset to use.

    chunk_size: int, default=None
        The number of vectors per chunk. If None (default) looks-up in
        scikit-learn configuration for `pairwise_dist_chunk_size`,
        and use 256 if it is not set.

    strategy : str, {'auto', 'parallel_on_X', 'parallel_on_Y'}, default=None
        The chunking strategy defining which dataset parallelization are made on.

        For both strategies the computations happens with two nested loops,
        respectively on chunks of X and chunks of Y.
        Strategies differs on which loop (outer or inner) is made to run
        in parallel with the Cython `prange` construct:

          - 'parallel_on_X' dispatches chunks of X uniformly on threads.
          Each thread then iterates on all the chunks of Y. This strategy is
          embarrassingly parallel and comes with no datastructures synchronisation.

          - 'parallel_on_Y' dispatches chunks of Y uniformly on threads.
          Each thread processes all the chunks of X in turn. This strategy is
          a sequence of embarrassingly parallel subtasks (the inner loop on Y
          chunks) with intermediate datastructures synchronisation at each
          iteration of the sequential outer loop on X chunks.

          - 'auto' relies on a simple heuristic to choose between
          'parallel_on_X' and 'parallel_on_Y': when `X.shape[0]` is large enough,
          'parallel_on_X' is usually the most efficient strategy. When `X.shape[0]`
          is small but `Y.shape[0]` is large, 'parallel_on_Y' brings more opportunity
          for parallelism and is therefore more efficient despite the synchronization
          step at each iteration of the outer loop on chunks of `X`.

          - None (default) looks-up in scikit-learn configuration for
          `pairwise_dist_parallel_strategy`, and use 'auto' if it is not set.
    """

    cdef:
        readonly DatasetsPair datasets_pair

        # The number of threads that can be used is stored in effective_n_threads.
        #
        # The number of threads to use in the parallelisation strategy
        # (i.e. parallel_on_X or parallel_on_Y) can be smaller than effective_n_threads:
        # for small datasets, less threads might be needed to loop over pair of chunks.
        #
        # Hence the number of threads that _will_ be used for looping over chunks
        # is stored in chunks_n_threads, allowing solely using what we need.
        #
        # Thus, an invariant is:
        #
        #                 chunks_n_threads <= effective_n_threads
        #
        ITYPE_t effective_n_threads
        ITYPE_t chunks_n_threads

        ITYPE_t n_samples_chunk, chunk_size

        ITYPE_t n_samples_X, X_n_samples_chunk, X_n_chunks, X_n_samples_last_chunk
        ITYPE_t n_samples_Y, Y_n_samples_chunk, Y_n_chunks, Y_n_samples_last_chunk

        bint execute_in_parallel_on_Y

    @classmethod
    def valid_metrics(cls) -> List[str]:
        excluded = {
            "pyfunc",  # is relatively slow because we need to coerce data as np arrays
            "mahalanobis", # is numerically unstable
            # TODO: In order to support discrete distance metrics, we need to have a
            # stable simultaneous sort which preserves the order of the input.
            # The best might be using std::stable_sort and a Comparator taking an
            # Arrays of Structures instead of Structure of Arrays (currently used).
            "hamming",
            *BOOL_METRICS,
        }
        return sorted(set(METRIC_MAPPING.keys()) - excluded)

    @classmethod
    def is_usable_for(cls, X, Y, metric) -> bool:
        """Return True if the PairwiseDistancesReduction can be used for the given parameters.

        Parameters
        ----------
        X : {ndarray, sparse matrix} of shape (n_samples_X, n_features)
            Input data.

        Y : {ndarray, sparse matrix} of shape (n_samples_Y, n_features)
            Input data.

        metric : str, default='euclidean'
            The distance metric to use.
            For a list of available metrics, see the documentation of
            :class:`~sklearn.metrics.DistanceMetric`.

        Returns
        -------
        True if the PairwiseDistancesReduction can be used, else False.
        """
        # TODO: support sparse arrays and 32 bits
        return (get_config().get("enable_cython_pairwise_dist", True) and
                not issparse(X) and X.dtype == np.float64 and
                not issparse(Y) and Y.dtype == np.float64 and
                metric in cls.valid_metrics())

    def __init__(
        self,
        DatasetsPair datasets_pair,
        chunk_size=None,
        strategy=None,
     ):
        cdef:
            ITYPE_t n_samples_chunk, X_n_full_chunks, Y_n_full_chunks

        if chunk_size is None:
            chunk_size = get_config().get("pairwise_dist_chunk_size", 256)

        self.chunk_size = check_scalar(chunk_size, "chunk_size", Integral, min_val=20)

        self.effective_n_threads = _openmp_effective_n_threads()

        self.datasets_pair = datasets_pair

        self.n_samples_X = datasets_pair.n_samples_X()
        self.X_n_samples_chunk = min(self.n_samples_X, self.chunk_size)
        X_n_full_chunks = self.n_samples_X // self.X_n_samples_chunk
        X_n_samples_remainder = self.n_samples_X % self.X_n_samples_chunk
        self.X_n_chunks = X_n_full_chunks + (X_n_samples_remainder != 0)

        if X_n_samples_remainder != 0:
            self.X_n_samples_last_chunk = X_n_samples_remainder
        else:
            self.X_n_samples_last_chunk = self.X_n_samples_chunk

        self.n_samples_Y = datasets_pair.n_samples_Y()
        self.Y_n_samples_chunk = min(self.n_samples_Y, self.chunk_size)
        Y_n_full_chunks = self.n_samples_Y // self.Y_n_samples_chunk
        Y_n_samples_remainder = self.n_samples_Y % self.Y_n_samples_chunk
        self.Y_n_chunks = Y_n_full_chunks + (Y_n_samples_remainder != 0)

        if Y_n_samples_remainder != 0:
            self.Y_n_samples_last_chunk = Y_n_samples_remainder
        else:
            self.Y_n_samples_last_chunk = self.Y_n_samples_chunk

        if strategy is None:
            strategy = get_config().get("pairwise_dist_parallel_strategy", 'auto')

        if strategy not in ('parallel_on_X', 'parallel_on_Y', 'auto'):
            raise RuntimeError(f"strategy must be 'parallel_on_X, 'parallel_on_Y', "
                               f"or 'auto', but currently strategy='{self.strategy}'.")

        if strategy == 'auto':
            # This is a simple heuristic whose constant for the
            # comparison has been chosen based on experiments.
            if 4 * self.chunk_size * self.effective_n_threads < self.n_samples_X:
                strategy = 'parallel_on_X'
            else:
                strategy = 'parallel_on_Y'

        self.execute_in_parallel_on_Y = strategy == "parallel_on_Y"

        # Not using less, not using more.
        self.chunks_n_threads = min(
            self.Y_n_chunks if self.execute_in_parallel_on_Y else self.X_n_chunks,
            self.effective_n_threads,
        )

    @final
    cdef void _parallel_on_X(self) nogil:
        """Compute the pairwise distances of each row vector of X on Y
        by parallelizing computation on the outer loop on chunks of X
        and reduce them.

        This strategy dispatches chunks of Y uniformly on threads.
        Each thread processes all the chunks of X in turn. This strategy is
        a sequence of embarrassingly parallel subtasks (the inner loop on Y
        chunks) with intermediate datastructures synchronisation at each
        iteration of the sequential outer loop on X chunks.

        Private datastructures are modified internally by threads.

        Private template methods can be implemented on subclasses to
        interact with those datastructures at various stages.
        """
        cdef:
            ITYPE_t Y_start, Y_end, X_start, X_end, X_chunk_idx, Y_chunk_idx
            ITYPE_t thread_num

        with nogil, parallel(num_threads=self.chunks_n_threads):
            thread_num = _openmp_thread_num()

            # Allocating thread datastructures
            self._parallel_on_X_parallel_init(thread_num)

            for X_chunk_idx in prange(self.X_n_chunks, schedule='static'):
                X_start = X_chunk_idx * self.X_n_samples_chunk
                if X_chunk_idx == self.X_n_chunks - 1:
                    X_end = X_start + self.X_n_samples_last_chunk
                else:
                    X_end = X_start + self.X_n_samples_chunk

                # Reinitializing thread datastructures for the new X chunk
                self._parallel_on_X_init_chunk(thread_num, X_start)

                for Y_chunk_idx in range(self.Y_n_chunks):
                    Y_start = Y_chunk_idx * self.Y_n_samples_chunk
                    if Y_chunk_idx == self.Y_n_chunks - 1:
                        Y_end = Y_start + self.Y_n_samples_last_chunk
                    else:
                        Y_end = Y_start + self.Y_n_samples_chunk

                    self._compute_and_reduce_distances_on_chunks(
                        X_start, X_end,
                        Y_start, Y_end,
                        thread_num,
                    )

                # Adjusting thread datastructures on the full pass on Y
                self._parallel_on_X_prange_iter_finalize(thread_num, X_start, X_end)

            # end: for X_chunk_idx

            # Deallocating thread datastructures
            self._parallel_on_X_parallel_finalize(thread_num)

        # end: with nogil, parallel
        return

    @final
    cdef void _parallel_on_Y(self) nogil:
        """Compute the pairwise distances of each row vector of X on Y
        by parallelizing computation on the inner loop on chunks of Y
        and reduce them.

        This strategy dispatches chunks of Y uniformly on threads.
        Each thread processes all the chunks of X in turn. This strategy is
        a sequence of embarrassingly parallel subtasks (the inner loop on Y
        chunks) with intermediate datastructures synchronisation at each
        iteration of the sequential outer loop on X chunks.

        Private datastructures are modified internally by threads.

        Private template methods can be implemented on subclasses to
        interact with those datastructures at various stages.
        """
        cdef:
            ITYPE_t Y_start, Y_end, X_start, X_end, X_chunk_idx, Y_chunk_idx
            ITYPE_t thread_num

        # Allocating datastructures shared by all threads
        self._parallel_on_Y_init()

        for X_chunk_idx in range(self.X_n_chunks):
            X_start = X_chunk_idx * self.X_n_samples_chunk
            if X_chunk_idx == self.X_n_chunks - 1:
                X_end = X_start + self.X_n_samples_last_chunk
            else:
                X_end = X_start + self.X_n_samples_chunk

            with nogil, parallel(num_threads=self.chunks_n_threads):
                thread_num = _openmp_thread_num()

                # Initializing datastructures used in this thread
                self._parallel_on_Y_parallel_init(thread_num)

                for Y_chunk_idx in prange(self.Y_n_chunks, schedule='static'):
                    Y_start = Y_chunk_idx * self.Y_n_samples_chunk
                    if Y_chunk_idx == self.Y_n_chunks - 1:
                        Y_end = Y_start + self.Y_n_samples_last_chunk
                    else:
                        Y_end = Y_start + self.Y_n_samples_chunk

                    self._compute_and_reduce_distances_on_chunks(
                        X_start, X_end,
                        Y_start, Y_end,
                        thread_num,
                    )
                # end: prange

                # Note: we don't need a _parallel_on_Y_finalize similarly.
                # This can be introduced if needed.

            # end: with nogil, parallel

            # Synchronizing the thread datastructures with the main ones
            self._parallel_on_Y_synchronize(X_start, X_end)

        # end: for X_chunk_idx
        # Deallocating temporary datastructures and adjusting main datastructures
        self._parallel_on_Y_finalize()
        return

    # Placeholder methods which have to be implemented

    cdef void _compute_and_reduce_distances_on_chunks(
        self,
        ITYPE_t X_start,
        ITYPE_t X_end,
        ITYPE_t Y_start,
        ITYPE_t Y_end,
        ITYPE_t thread_num,
    ) nogil:
        """Compute the pairwise distances on two chunks of X and Y and reduce them.

        This is THE core computational method of PairwiseDistanceReductions.
        This must be implemented in subclasses.
        """
        return

    def _finalize_results(self, bint return_distance):
        """Callback adapting datastructures before returning results.

        This must be implemented in subclasses.
        """
        return None

    # Placeholder methods which can be implemented

    cdef void compute_exact_distances(self) nogil:
        """Convert rank-preserving distances to exact distances or recompute them."""
        return

    cdef void _parallel_on_X_parallel_init(
        self,
        ITYPE_t thread_num,
    ) nogil:
        """Allocate datastructures used in a thread given its number."""
        return

    cdef void _parallel_on_X_init_chunk(
        self,
        ITYPE_t thread_num,
        ITYPE_t X_start,
    ) nogil:
        """Initialise datastructures used in a thread given its number."""
        return

    cdef void _parallel_on_X_prange_iter_finalize(
        self,
        ITYPE_t thread_num,
        ITYPE_t X_start,
        ITYPE_t X_end,
    ) nogil:
        """Interact with datastructures after a reduction on chunks."""
        return

    cdef void _parallel_on_X_parallel_finalize(
        self,
        ITYPE_t thread_num
    ) nogil:
        """Interact with datastructures after executing all the reductions."""
        return

    cdef void _parallel_on_Y_init(
        self,
    ) nogil:
        """Allocate datastructures used in all threads."""
        return

    cdef void _parallel_on_Y_parallel_init(
        self,
        ITYPE_t thread_num,
    ) nogil:
        """Initialise datastructures used in a thread given its number."""
        return

    cdef void _parallel_on_Y_synchronize(
        self,
        ITYPE_t X_start,
        ITYPE_t X_end,
    ) nogil:
        """Update thread datastructures before leaving a parallel region."""
        return

    cdef void _parallel_on_Y_finalize(
        self,
    ) nogil:
        """Update datastructures after executing all the reductions."""
        return

cdef class PairwiseDistancesArgKmin(PairwiseDistancesReduction):
    """Compute the argkmin of row vectors of X on the ones of Y.

    For each row vector of X, computes the indices of k first the rows
    vectors of Y with the smallest distances.

    PairwiseDistancesArgKmin is typically used to perform
    bruteforce k-nearest neighbors queries.

    Parameters
    ----------
    datasets_pair: DatasetsPair
        The dataset pairs (X, Y) for the reduction.

    chunk_size: int, default=None,
        The number of vectors per chunk. If None (default) looks-up in
        scikit-learn configuration for `pairwise_dist_chunk_size`,
        and use 256 if it is not set.

    k: int, default=1
        The k for the argkmin reduction.
    """

    cdef:
        ITYPE_t k

        ITYPE_t[:, ::1] argkmin_indices
        DTYPE_t[:, ::1] argkmin_distances

        # Used as array of pointers to private datastructures used in threads.
        DTYPE_t ** heaps_r_distances_chunks
        ITYPE_t ** heaps_indices_chunks

    @classmethod
    def compute(
        cls,
        X,
        Y,
        ITYPE_t k,
        str metric="euclidean",
        chunk_size=None,
        dict metric_kwargs=None,
        str strategy=None,
        bint return_distance=False,
    ):
        """Return the results of the reduction for the given arguments.

        Parameters
        ----------
        X : ndarray or CSR matrix of shape (n_samples_X, n_features)
            Input data.

        Y : ndarray or CSR matrix of shape (n_samples_Y, n_features)
            Input data.

        k : int
            The k for the argkmin reduction.

        metric : str, default='euclidean'
            The distance metric to use for argkmin.
            For a list of available metrics, see the documentation of
            :class:`~sklearn.metrics.DistanceMetric`.

        chunk_size : int, default=None,
            The number of vectors per chunk. If None (default) looks-up in
            scikit-learn configuration for `pairwise_dist_chunk_size`,
            and use 256 if it is not set.

        metric_kwargs : dict, default=None
            Keyword arguments to pass to specified metric function.

        strategy : str, {'auto', 'parallel_on_X', 'parallel_on_Y'}, default=None
            The chunking strategy defining which dataset parallelization are made on.

            For both strategies the computations happens with two nested loops,
            respectively on chunks of X and chunks of Y.
            Strategies differs on which loop (outer or inner) is made to run
            in parallel with the Cython `prange` construct:

              - 'parallel_on_X' dispatches chunks of X uniformly on threads.
              Each thread then iterates on all the chunks of Y. This strategy is
              embarrassingly parallel and comes with no datastructures synchronisation.

              - 'parallel_on_Y' dispatches chunks of Y uniformly on threads.
              Each thread processes all the chunks of X in turn. This strategy is
              a sequence of embarrassingly parallel subtasks (the inner loop on Y
              chunks) with intermediate datastructures synchronisation at each
              iteration of the sequential outer loop on X chunks.

              - 'auto' relies on a simple heuristic to choose between
              'parallel_on_X' and 'parallel_on_Y': when `X.shape[0]` is large enough,
              'parallel_on_X' is usually the most efficient strategy. When `X.shape[0]`
              is small but `Y.shape[0]` is large, 'parallel_on_Y' brings more opportunity
              for parallelism and is therefore more efficient despite the synchronization
              step at each iteration of the outer loop on chunks of `X`.

              - None (default) looks-up in scikit-learn configuration for
              `pairwise_dist_parallel_strategy`, and use 'auto' if it is not set.

        return_distance : boolean, default=False
            Return distances between each X vector and its
            argkmin if set to True.

        Returns
        -------
            If return_distance=False:
              - argkmin_indices : ndarray of shape (n_samples_X, k)
                Indices of the argkmin for each vector in X.

            If return_distance=True:
              - argkmin_distances : ndarray of shape (n_samples_X, k)
                Distances to the argkmin for each vector in X.
              - argkmin_indices : ndarray of shape (n_samples_X, k)
                Indices of the argkmin for each vector in X.

        Notes
        -----
            This public classmethod is responsible for introspecting the arguments
            values to dispatch to the private :meth:`PairwiseDistancesArgKmin._compute`
            instance method of the most appropriate :class:`PairwiseDistancesArgKmin`
            concrete implementation.

            All temporarily allocated datastructures necessary for the concrete
            implementation are therefore freed when this classmethod returns.

            This allows entirely decoupling the interface entirely from the
            implementation details whilst maintaining RAII.
        """
        # Note (jjerphan): Some design thoughts for future extensions.
        # This factory comes to handle specialisations for the given arguments.
        # For future work, this might can be an entrypoint to specialise operations
        # for various backend and/or hardware and/or datatypes, and/or fused
        # {sparse, dense}-datasetspair etc.
        if (
            metric in ("euclidean", "sqeuclidean")
            and not issparse(X)
            and not issparse(Y)
        ):
            # Specialized implementation with improved arithmetic intensity
            # and vector instructions (SIMD) by processing several vectors
            # at time to leverage a call to the BLAS GEMM routine as explained
            # in more details in the docstring.
            use_squared_distances = metric == "sqeuclidean"
            pda = FastEuclideanPairwiseDistancesArgKmin(
                X=X, Y=Y, k=k,
                use_squared_distances=use_squared_distances,
                chunk_size=chunk_size,
                strategy=strategy,
                metric_kwargs=metric_kwargs,
            )
        else:
             # Fall back on a generic implementation that handles most scipy
             # metrics by computing the distances between 2 vectors at a time.
            pda = PairwiseDistancesArgKmin(
                datasets_pair=DatasetsPair.get_for(X, Y, metric, metric_kwargs),
                k=k,
                chunk_size=chunk_size,
                strategy=strategy,
            )

        # Limit the number of threads in second level of nested parallelism for BLAS
        # to avoid threads over-subscription (in GEMM for instance).
        with threadpool_limits(limits=1, user_api="blas"):
            if pda.execute_in_parallel_on_Y:
                pda._parallel_on_Y()
            else:
                pda._parallel_on_X()

        return pda._finalize_results(return_distance)

    def __init__(
        self,
        DatasetsPair datasets_pair,
        chunk_size=None,
        strategy=None,
        ITYPE_t k=1,
    ):
        super().__init__(
            datasets_pair=datasets_pair,
            chunk_size=chunk_size,
            strategy=strategy,
        )
        self.k = check_scalar(k, "k", Integral, min_val=1)

        # Allocating pointers to datastructures but not the datastructures themselves.
        # There are as many pointers as effective threads.
        #
        # For the sake of explicitness:
        #   - when parallelizing on X, the pointers of those heaps are referencing
        #   (with proper offsets) addresses of the two main heaps (see bellow)
        #   - when parallelizing on Y, the pointers of those heaps are referencing
        #   small heaps which are thread-wise-allocated and whose content will be
        #   merged with the main heaps'.
        self.heaps_r_distances_chunks = <DTYPE_t **> malloc(
            sizeof(DTYPE_t *) * self.chunks_n_threads
        )
        self.heaps_indices_chunks = <ITYPE_t **> malloc(
            sizeof(ITYPE_t *) * self.chunks_n_threads
        )

        # Main heaps which will be returned as results by `PairwiseDistancesArgKmin.compute`.
        self.argkmin_indices = np.full((self.n_samples_X, self.k), 0, dtype=ITYPE)
        self.argkmin_distances = np.full((self.n_samples_X, self.k), DBL_MAX, dtype=DTYPE)

    def __dealloc__(self):
        if self.heaps_indices_chunks is not NULL:
            free(self.heaps_indices_chunks)

        if self.heaps_r_distances_chunks is not NULL:
            free(self.heaps_r_distances_chunks)

    cdef void _compute_and_reduce_distances_on_chunks(
        self,
        ITYPE_t X_start,
        ITYPE_t X_end,
        ITYPE_t Y_start,
        ITYPE_t Y_end,
        ITYPE_t thread_num,
    ) nogil:
        cdef:
            ITYPE_t i, j
            ITYPE_t n_samples_X = X_end - X_start
            ITYPE_t n_samples_Y = Y_end - Y_start
            DTYPE_t *heaps_r_distances = self.heaps_r_distances_chunks[thread_num]
            ITYPE_t *heaps_indices = self.heaps_indices_chunks[thread_num]

        # Pushing the distances and their associated indices on a heap
        # which by construction will keep track of the argkmin.
        for i in range(n_samples_X):
            for j in range(n_samples_Y):
                heap_push(
                    heaps_r_distances + i * self.k,
                    heaps_indices + i * self.k,
                    self.k,
                    self.datasets_pair.surrogate_dist(X_start + i, Y_start + j),
                    Y_start + j,
                )

    @final
    cdef void _parallel_on_X_init_chunk(
        self,
        ITYPE_t thread_num,
        ITYPE_t X_start,
    ) nogil:
        # As this strategy is embarrassingly parallel, we can set each
        # thread's heaps pointer to the proper position on the main heaps.
        self.heaps_r_distances_chunks[thread_num] = &self.argkmin_distances[X_start, 0]
        self.heaps_indices_chunks[thread_num] = &self.argkmin_indices[X_start, 0]

    @final
    cdef void _parallel_on_X_prange_iter_finalize(
        self,
        ITYPE_t thread_num,
        ITYPE_t X_start,
        ITYPE_t X_end,
    ) nogil:
        cdef:
            ITYPE_t idx, jdx

        # Sorting the main heaps portion associated to `X[X_start:X_end]`
        # in ascending order w.r.t the distances.
        for idx in range(X_end - X_start):
            simultaneous_sort(
                self.heaps_r_distances_chunks[thread_num] + idx * self.k,
                self.heaps_indices_chunks[thread_num] + idx * self.k,
                self.k
            )

    cdef void _parallel_on_Y_init(
        self,
    ) nogil:
        cdef:
            # Maximum number of scalar elements (the last chunks can be smaller)
            ITYPE_t heaps_size = self.X_n_samples_chunk * self.k
            ITYPE_t thread_num

        # The allocation is done in parallel for data locality purposes: this way
        # the heaps used in each threads are allocated in pages which are closer
        # to the CPU core used by the thread.
        # See comments about First Touch Placement Policy:
        # https://www.openmp.org/wp-content/uploads/openmp-webinar-vanderPas-20210318.pdf #noqa
        for thread_num in prange(self.chunks_n_threads, schedule='static', nogil=True,
                                 num_threads=self.chunks_n_threads):
            # As chunks of X are shared across threads, so must their
            # heaps. To solve this, each thread has its own heaps
            # which are then synchronised back in the main ones.
            self.heaps_r_distances_chunks[thread_num] = <DTYPE_t *> malloc(
                heaps_size * sizeof(DTYPE_t)
            )
            self.heaps_indices_chunks[thread_num] = <ITYPE_t *> malloc(
                heaps_size * sizeof(ITYPE_t)
            )

    @final
    cdef void _parallel_on_Y_parallel_init(
        self,
        ITYPE_t thread_num,
    ) nogil:
        # Initialising heaps (memset can't be used here)
        for idx in range(self.X_n_samples_chunk * self.k):
            self.heaps_r_distances_chunks[thread_num][idx] = DBL_MAX
            self.heaps_indices_chunks[thread_num][idx] = -1

    @final
    cdef void _parallel_on_Y_synchronize(
        self,
        ITYPE_t X_start,
        ITYPE_t X_end,
    ) nogil:
        cdef:
            ITYPE_t idx, jdx, thread_num
        with nogil, parallel(num_threads=self.effective_n_threads):
            # Synchronising the thread heaps with the main heaps.
            # This is done in parallel sample-wise (no need for locks).
            #
            # This might break each thread's data locality as each heap which
            # was allocated in a thread is being now being used in several threads.
            #
            # Still, this parallel pattern has shown to be efficient in practice.
            for idx in prange(X_end - X_start, schedule="static"):
                for thread_num in range(self.chunks_n_threads):
                    for jdx in range(self.k):
                        heap_push(
                            &self.argkmin_distances[X_start + idx, 0],
                            &self.argkmin_indices[X_start + idx, 0],
                            self.k,
                            self.heaps_r_distances_chunks[thread_num][idx * self.k + jdx],
                            self.heaps_indices_chunks[thread_num][idx * self.k + jdx],
                        )

    cdef void _parallel_on_Y_finalize(
        self,
    ) nogil:
        cdef:
            ITYPE_t idx, thread_num

        with nogil, parallel(num_threads=self.chunks_n_threads):
            # Deallocating temporary datastructures
            for thread_num in prange(self.chunks_n_threads, schedule='static'):
                free(self.heaps_r_distances_chunks[thread_num])
                free(self.heaps_indices_chunks[thread_num])

            # Sorting the main in ascending order w.r.t the distances.
            # This is done in parallel sample-wise (no need for locks).
            for idx in prange(self.n_samples_X, schedule='static'):
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
        for i in prange(self.n_samples_X, schedule='static', nogil=True,
                        num_threads=self.effective_n_threads):
            for j in range(self.k):
                distances[i, j] = self.datasets_pair.distance_metric._rdist_to_dist(
                    # Guard against eventual -0., causing nan production.
                    max(distances[i, j], 0.)
                )

    def _finalize_results(self, bint return_distance=False):
        if return_distance:
            # We need to recompute distances because we relied on
            # surrogate distances for the reduction.
            self.compute_exact_distances()

            # Values are returned identically to the way `KNeighborsMixin.kneighbors`
            # returns values. This is counter-intuitive but this allows not using
            # complex adaptations where `PairwiseDistancesArgKmin.compute` is called.
            return np.asarray(self.argkmin_distances), np.asarray(self.argkmin_indices)

        return np.asarray(self.argkmin_indices)


cdef class GEMMTermComputer:
    """Component for `FastEuclidean*` variant wrapping the logic for the call to GEMM.

    `FastEuclidean*` classes internally compute the squared Euclidean distances between
    chunks of vectors X_c and Y_c using using the decomposition:


                ||X_c_i - Y_c_j||² = ||X_c_i||² - 2 X_c_i.Y_c_j^T + ||Y_c_j||²


    This helper class is in charge of wrapping the common logic to compute
    the middle term `- 2 X_c_i.Y_c_j^T` with a call to GEMM, which has a high
    arithmetic intensity.
    """

    cdef:
        const DTYPE_t[:, ::1] X
        const DTYPE_t[:, ::1] Y

        ITYPE_t effective_n_threads
        ITYPE_t chunks_n_threads
        ITYPE_t dist_middle_terms_chunks_size

        # Buffers for the `-2 * X_c @ Y_c.T` term computed via GEMM
        vector[vector[DTYPE_t]] dist_middle_terms_chunks

    def __init__(self,
        DTYPE_t[:, ::1] X,
        DTYPE_t[:, ::1] Y,
        ITYPE_t effective_n_threads,
        ITYPE_t chunks_n_threads,
        ITYPE_t dist_middle_terms_chunks_size,
    ):
        self.X = X
        self.Y = Y
        self.effective_n_threads = effective_n_threads
        self.chunks_n_threads = chunks_n_threads
        self.dist_middle_terms_chunks_size = dist_middle_terms_chunks_size

        self.dist_middle_terms_chunks = vector[vector[DTYPE_t]](self.effective_n_threads)

    cdef void _parallel_on_X_parallel_init(self, ITYPE_t thread_num) nogil:
        self.dist_middle_terms_chunks[thread_num].resize(self.dist_middle_terms_chunks_size)

    cdef void _parallel_on_Y_init(self) nogil:
        for thread_num in range(self.chunks_n_threads):
            self.dist_middle_terms_chunks[thread_num].resize(
                self.dist_middle_terms_chunks_size
            )

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
            ITYPE_t lda = X_c.shape[1]
            DTYPE_t * B = <DTYPE_t *> &Y_c[0, 0]
            ITYPE_t ldb = X_c.shape[1]
            DTYPE_t beta = 0.
            ITYPE_t ldc = Y_c.shape[0]

        # dist_middle_terms = `-2 * X_c @ Y_c.T`
        _gemm(order, ta, tb, m, n, K, alpha, A, lda, B, ldb, beta, dist_middle_terms, ldc)

        return dist_middle_terms


cdef class FastEuclideanPairwiseDistancesArgKmin(PairwiseDistancesArgKmin):
    """Fast specialized variant for PairwiseDistancesArgKmin on EuclideanDistance.

    The full pairwise squared distances matrix is computed as follows:

                  ||X - Y||² = ||X||² - 2 X.Y^T + ||Y||²

    The middle term gets computed efficiently bellow using BLAS Level 3 GEMM.

    Notes
    -----
    This implementation has a superior arithmetic intensity and hence
    better running time when the variant is IO bound, but it can suffer
    from numerical instability caused by catastrophic cancellation potentially
    introduced by the subtraction in the arithmetic expression above.
    """

    cdef:
        GEMMTermComputer gemm_term_computer
        const DTYPE_t[::1] X_norm_squared
        const DTYPE_t[::1] Y_norm_squared

        bint use_squared_distances

    @classmethod
    def is_usable_for(cls, X, Y, metric) -> bool:
        return (PairwiseDistancesArgKmin.is_usable_for(X, Y, metric) and
                not _in_unstable_openblas_configuration())

    def __init__(
        self,
        X,
        Y,
        ITYPE_t k,
        bint use_squared_distances=False,
        chunk_size=None,
        strategy=None,
        metric_kwargs=None,
    ):
        if (
            metric_kwargs is not None and
            len(metric_kwargs) > 0 and
            "Y_norm_squared" not in metric_kwargs
        ):
            warnings.warn(
                f"Some metric_kwargs have been passed ({metric_kwargs}) but aren't "
                f"usable for this case ({self.__class__.__name__}) and will be ignored.",
                UserWarning,
                stacklevel=3,
            )

        super().__init__(
            # The datasets pair here is used for exact distances computations
            datasets_pair=DatasetsPair.get_for(X, Y, metric="euclidean"),
            chunk_size=chunk_size,
            strategy=strategy,
            k=k,
        )
        # X and Y are checked by the DatasetsPair implemented as a DenseDenseDatasetsPair
        cdef:
            DenseDenseDatasetsPair datasets_pair = <DenseDenseDatasetsPair> self.datasets_pair
            ITYPE_t dist_middle_terms_chunks_size = self.Y_n_samples_chunk * self.X_n_samples_chunk

        self.gemm_term_computer = GEMMTermComputer(
            datasets_pair.X,
            datasets_pair.Y,
            self.effective_n_threads,
            self.chunks_n_threads,
            dist_middle_terms_chunks_size,
        )

        if metric_kwargs is not None and "Y_norm_squared" in metric_kwargs:
            self.Y_norm_squared = metric_kwargs.pop("Y_norm_squared")
        else:
            self.Y_norm_squared = _sqeuclidean_row_norms(datasets_pair.Y, self.effective_n_threads)

        # Do not recompute norms if datasets are identical.
        self.X_norm_squared = (
            self.Y_norm_squared if X is Y else
            _sqeuclidean_row_norms(datasets_pair.X, self.effective_n_threads)
        )
        self.use_squared_distances = use_squared_distances

    @final
    cdef void compute_exact_distances(self) nogil:
        if not self.use_squared_distances:
            PairwiseDistancesArgKmin.compute_exact_distances(self)

    @final
    cdef void _parallel_on_X_parallel_init(
        self,
        ITYPE_t thread_num,
    ) nogil:
        PairwiseDistancesArgKmin._parallel_on_X_parallel_init(self, thread_num)
        self.gemm_term_computer._parallel_on_X_parallel_init(thread_num)

    @final
    cdef void _parallel_on_Y_init(
        self,
    ) nogil:
        cdef ITYPE_t thread_num
        PairwiseDistancesArgKmin._parallel_on_Y_init(self)
        self.gemm_term_computer._parallel_on_Y_init()

    @final
    cdef void _compute_and_reduce_distances_on_chunks(
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
            ITYPE_t n_X = X_end - X_start
            ITYPE_t n_Y = Y_end - Y_start
            DTYPE_t * dist_middle_terms = self.gemm_term_computer._compute_distances_on_chunks(
                X_start, X_end, Y_start, Y_end, thread_num
            )
            DTYPE_t * heaps_r_distances = self.heaps_r_distances_chunks[thread_num]
            ITYPE_t * heaps_indices = self.heaps_indices_chunks[thread_num]


        # Pushing the distance and their associated indices on heaps
        # which keep tracks of the argkmin.
        for i in range(n_X):
            for j in range(n_Y):
                heap_push(
                    heaps_r_distances + i * self.k,
                    heaps_indices + i * self.k,
                    self.k,
                    # Using the squared euclidean distance as the rank-preserving distance:
                    #
                    #             ||X_c_i||² - 2 X_c_i.Y_c_j^T + ||Y_c_j||²
                    #
                    (
                        self.X_norm_squared[i + X_start] +
                        dist_middle_terms[i * n_Y + j] +
                        self.Y_norm_squared[j + Y_start]
                    ),
                    j + Y_start,
                )


cdef class PairwiseDistancesRadiusNeighborhood(PairwiseDistancesReduction):
    """Compute radius-based neighbors for two sets of vectors.

    For each row-vector X[i] of the queries X, find all the indices j of
    row-vectors in Y such that:

                        dist(X[i], Y[j]) <= radius

    The distance function `dist` depends on the values of the `metric`
    and `metric_kwargs` parameters.

    Parameters
    ----------
    datasets_pair: DatasetsPair
        The dataset pair (X, Y) for the reduction.

    chunk_size: int, default=None,
        The number of vectors per chunk. If None (default) looks-up in
        scikit-learn configuration for `pairwise_dist_chunk_size`,
        and use 256 if it is not set.

    radius: float
        The radius defining the neighborhood.
    """

    cdef:
        DTYPE_t radius

        # DistanceMetric compute rank-preserving surrogate distance via rdist
        # which are proxies necessitating less computations.
        # We get the equivalent for the radius to be able to compare it against
        # vectors' rank-preserving surrogate distances.
        DTYPE_t r_radius

        # Neighbors indices and distances are returned as np.ndarrays of np.ndarrays.
        #
        # For this implementation, we want resizable buffers which we will wrap
        # into numpy arrays at the end. std::vector comes as a handy interface
        # for interacting efficiently with resizable buffers.
        #
        # Though it is possible to access their buffer address with
        # std::vector::data, they can't be stolen: buffers lifetime
        # is tied to their std::vector and are deallocated when
        # std::vectors are.
        #
        # To solve this, we dynamically allocate std::vectors and then
        # encapsulate them in a StdVectorSentinel responsible for
        # freeing them when the associated np.ndarray is freed.
        #
        # Shared pointers (defined via shared_ptr) are use for safer memory management.
        # Unique pointers (defined via unique_ptr) can't be used as datastructures
        # are shared across threads for parallel_on_X; see _parallel_on_X_init_chunk.
        shared_ptr[vector[vector[ITYPE_t]]] neigh_indices
        shared_ptr[vector[vector[DTYPE_t]]] neigh_distances

        # Used as array of pointers to private datastructures used in threads.
        vector[shared_ptr[vector[vector[ITYPE_t]]]] neigh_indices_chunks
        vector[shared_ptr[vector[vector[DTYPE_t]]]] neigh_distances_chunks

        bint sort_results

    @classmethod
    def compute(
        cls,
        X,
        Y,
        DTYPE_t radius,
        str metric="euclidean",
        chunk_size=None,
        dict metric_kwargs=None,
        str strategy=None,
        bint return_distance=False,
        bint sort_results=False,
    ):
        """Return the results of the reduction for the given arguments.

        Parameters
        ----------
        X : ndarray or CSR matrix of shape (n_samples_X, n_features)
            Input data.

        Y : ndarray or CSR matrix of shape (n_samples_Y, n_features)
            Input data.

        radius : float
            The radius defining the neighborhood.

        metric : str, default='euclidean'
            The distance metric to use.
            For a list of available metrics, see the documentation of
            :class:`~sklearn.metrics.DistanceMetric`.

        chunk_size : int, default=None,
            The number of vectors per chunk. If None (default) looks-up in
            scikit-learn configuration for `pairwise_dist_chunk_size`,
            and use 256 if it is not set.

        metric_kwargs : dict, default=None
            Keyword arguments to pass to specified metric function.

        strategy : str, {'auto', 'parallel_on_X', 'parallel_on_Y'}, default=None
            The chunking strategy defining which dataset parallelization are made on.

            For both strategies the computations happens with two nested loops,
            respectively on chunks of X and chunks of Y.
            Strategies differs on which loop (outer or inner) is made to run
            in parallel with the Cython `prange` construct:

              - 'parallel_on_X' dispatches chunks of X uniformly on threads.
              Each thread then iterates on all the chunks of Y. This strategy is
              embarrassingly parallel and comes with no datastructures synchronisation.

              - 'parallel_on_Y' dispatches chunks of Y uniformly on threads.
              Each thread processes all the chunks of X in turn. This strategy is
              a sequence of embarrassingly parallel subtasks (the inner loop on Y
              chunks) with intermediate datastructures synchronisation at each
              iteration of the sequential outer loop on X chunks.

              - 'auto' relies on a simple heuristic to choose between
              'parallel_on_X' and 'parallel_on_Y': when `X.shape[0]` is large enough,
              'parallel_on_X' is usually the most efficient strategy. When `X.shape[0]`
              is small but `Y.shape[0]` is large, 'parallel_on_Y' brings more opportunity
              for parallelism and is therefore more efficient despite the synchronization
              step at each iteration of the outer loop on chunks of `X`.

              - None (default) looks-up in scikit-learn configuration for
              `pairwise_dist_parallel_strategy`, and use 'auto' if it is not set.

        return_distance : boolean, default=False
            Return distances between each X vector and its neighbors if set to True.

        sort_results : boolean, default=False
            Sort results with respect to distances between each X vector and its
            neighbors if set to True.

        Returns
        -------
        If return_distance=False:
          - neighbors_indices : ndarray of n_samples_X ndarray
            Indices of the neighbors for each vector in X.

        If return_distance=True:
          - neighbors_indices : ndarray of n_samples_X ndarray
            Indices of the neighbors for each vector in X.
          - neighbors_distances : ndarray of n_samples_X ndarray
            Distances to the neighbors for each vector in X.

        Notes
        -----
        This public classmethod is responsible for introspecting the arguments
        values to dispatch to the private
        :meth:`PairwiseDistancesRadiusNeighborhood._compute` instance method of
        the most appropriate :class:`PairwiseDistancesRadiusNeighborhood`
        concrete implementation.

        All temporarily allocated datastructures necessary for the concrete
        implementation are therefore freed when this classmethod returns.

        This allows entirely decoupling the interface entirely from the
        implementation details whilst maintaining RAII.
        """
        # Note (jjerphan): Some design thoughts for future extensions.
        # This factory comes to handle specialisations for the given arguments.
        # For future work, this might can be an entrypoint to specialise operations
        # for various backend and/or hardware and/or datatypes, and/or fused
        # {sparse, dense}-datasetspair etc.
        if (
            metric in ("euclidean", "sqeuclidean")
            and not issparse(X)
            and not issparse(Y)
        ):
            # Specialized implementation with improved arithmetic intensity
            # and vector instructions (SIMD) by processing several vectors
            # at time to leverage a call to the BLAS GEMM routine as explained
            # in more details in the docstring.
            use_squared_distances = metric == "sqeuclidean"
            pda = FastEuclideanPairwiseDistancesRadiusNeighborhood(
                X=X, Y=Y, radius=radius,
                use_squared_distances=use_squared_distances,
                chunk_size=chunk_size,
                metric_kwargs=metric_kwargs,
                strategy=strategy,
                sort_results=sort_results,
            )
        else:
             # Fall back on a generic implementation that handles most scipy
             # metrics by computing the distances between 2 vectors at a time.
            pda = PairwiseDistancesRadiusNeighborhood(
                datasets_pair=DatasetsPair.get_for(X, Y, metric, metric_kwargs),
                radius=radius,
                chunk_size=chunk_size,
                metric_kwargs=metric_kwargs,
                strategy=strategy,
                sort_results=sort_results,
            )

        # Limit the number of threads in second level of nested parallelism for BLAS
        # to avoid threads over-subscription (in GEMM for instance).
        with threadpool_limits(limits=1, user_api="blas"):
            if pda.execute_in_parallel_on_Y:
                pda._parallel_on_Y()
            else:
                pda._parallel_on_X()

        return pda._finalize_results(return_distance)


    def __init__(
        self,
        DatasetsPair datasets_pair,
        DTYPE_t radius,
        chunk_size=None,
        strategy=None,
        sort_results=False,
        metric_kwargs=None,
    ):
        super().__init__(
            datasets_pair=datasets_pair,
            chunk_size=chunk_size,
            strategy=strategy,
        )

        self.radius = check_scalar(radius, "radius", Real, min_val=0)
        self.r_radius = self.datasets_pair.distance_metric._dist_to_rdist(radius)
        self.sort_results = sort_results

        # Allocating pointers to datastructures but not the datastructures themselves.
        # There are as many pointers as effective threads.
        #
        # For the sake of explicitness:
        #   - when parallelizing on X, the pointers of those heaps are referencing
        #   self.neigh_distances and self.neigh_indices
        #   - when parallelizing on Y, the pointers of those heaps are referencing
        #   std::vectors of std::vectors which are thread-wise-allocated and whose
        #   content will be merged into self.neigh_distances and self.neigh_indices.
        self.neigh_distances_chunks = vector[shared_ptr[vector[vector[DTYPE_t]]]](
            self.chunks_n_threads
        )
        self.neigh_indices_chunks = vector[shared_ptr[vector[vector[ITYPE_t]]]](
            self.chunks_n_threads
        )

        # Temporary datastructures which will be coerced to numpy arrays on before
        # PairwiseDistancesRadiusNeighborhood.compute "return" and will be then freed.
        self.neigh_distances = make_shared[vector[vector[DTYPE_t]]](self.n_samples_X)
        self.neigh_indices = make_shared[vector[vector[ITYPE_t]]](self.n_samples_X)

    cdef void _compute_and_reduce_distances_on_chunks(
        self,
        ITYPE_t X_start,
        ITYPE_t X_end,
        ITYPE_t Y_start,
        ITYPE_t Y_end,
        ITYPE_t thread_num,
    ) nogil:
        cdef:
            ITYPE_t i, j
            DTYPE_t r_dist_i_j

        for i in range(X_start, X_end):
            for j in range(Y_start, Y_end):
                r_dist_i_j = self.datasets_pair.surrogate_dist(i, j)
                if r_dist_i_j <= self.r_radius:
                    deref(self.neigh_distances_chunks[thread_num])[i].push_back(r_dist_i_j)
                    deref(self.neigh_indices_chunks[thread_num])[i].push_back(j)

    def _finalize_results(self, bint return_distance=False):
        if return_distance:
            # We need to recompute distances because we relied on
            # surrogate distances for the reduction.
            self.compute_exact_distances()
            return (
                coerce_vectors_to_nd_arrays(self.neigh_distances),
                coerce_vectors_to_nd_arrays(self.neigh_indices),
            )

        return coerce_vectors_to_nd_arrays(self.neigh_indices)

    @final
    cdef void _parallel_on_X_init_chunk(
        self,
        ITYPE_t thread_num,
        ITYPE_t X_start,
    ) nogil:

        # As this strategy is embarrassingly parallel, we can set the
        # thread vectors' pointers to the main vectors'.
        self.neigh_distances_chunks[thread_num] = self.neigh_distances
        self.neigh_indices_chunks[thread_num] = self.neigh_indices

    @final
    cdef void _parallel_on_X_prange_iter_finalize(
        self,
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

    cdef void _parallel_on_Y_init(
        self,
    ) nogil:
        cdef:
            ITYPE_t thread_num
        # As chunks of X are shared across threads, so must datastructures to avoid race
        # conditions: each thread has its own vectors of n_samples_X vectors which are
        # then merged back in the main n_samples_X vectors.
        for thread_num in range(self.chunks_n_threads):
            self.neigh_distances_chunks[thread_num] = make_shared[vector[vector[DTYPE_t]]](self.n_samples_X)
            self.neigh_indices_chunks[thread_num] = make_shared[vector[vector[ITYPE_t]]](self.n_samples_X)

    @final
    cdef void _merge_vectors(
        self,
        ITYPE_t idx,
        ITYPE_t num_threads,
    ) nogil:
        cdef:
            ITYPE_t thread_num
            ITYPE_t idx_n_elements = 0
            ITYPE_t last_element_idx = deref(self.neigh_indices)[idx].size()

        # Resizing buffers only once for the given number of elements.
        for thread_num in range(num_threads):
            idx_n_elements += deref(self.neigh_distances_chunks[thread_num])[idx].size()

        deref(self.neigh_distances)[idx].resize(last_element_idx + idx_n_elements)
        deref(self.neigh_indices)[idx].resize(last_element_idx + idx_n_elements)

        # Moving the elements by range using the range first element
        # as the reference for the insertion.
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


    cdef void _parallel_on_Y_finalize(
        self,
    ) nogil:
        cdef:
            ITYPE_t idx, jdx, thread_num, idx_n_element, idx_current

        with nogil, parallel(num_threads=self.effective_n_threads):
            # Merge vectors used in threads into the main ones.
            # This is done in parallel sample-wise (no need for locks)
            # using dynamic scheduling because we might not have
            # the same number of neighbors for each query vector.
            # TODO: compare 'dynamic' vs 'static' vs 'guided'
            for idx in prange(self.n_samples_X, schedule='dynamic'):
                self._merge_vectors(idx, self.chunks_n_threads)

            # The content of the vector have been std::moved.
            # Hence they can't be used anymore and can be deleted.
            # Their deletion is carried out automatically as the
            # implementation relies on shared pointers.

            # Sort in parallel in ascending order w.r.t the distances if requested.
            if self.sort_results:
                for idx in prange(self.n_samples_X, schedule='static'):
                    simultaneous_sort(
                        deref(self.neigh_distances)[idx].data(),
                        deref(self.neigh_indices)[idx].data(),
                        deref(self.neigh_indices)[idx].size()
                    )

        return

    cdef void compute_exact_distances(self) nogil:
        """Convert rank-preserving distances to pairwise distances in parallel."""
        cdef:
            ITYPE_t i, j

        for i in prange(self.n_samples_X, nogil=True, schedule='dynamic',
                        num_threads=self.effective_n_threads):
            for j in range(deref(self.neigh_indices)[i].size()):
                deref(self.neigh_distances)[i][j] = (
                        self.datasets_pair.distance_metric._rdist_to_dist(
                            # Guard against eventual -0., causing nan production.
                            max(deref(self.neigh_distances)[i][j], 0.)
                        )
                )


cdef class FastEuclideanPairwiseDistancesRadiusNeighborhood(PairwiseDistancesRadiusNeighborhood):
    """Fast specialized variant for PairwiseDistancesRadiusNeighborhood on EuclideanDistance.

    The full pairwise squared distances matrix is computed as follows:

                  ||X - Y||² = ||X||² - 2 X.Y^T + ||Y||²

    The middle term gets computed efficiently bellow using BLAS Level 3 GEMM.

    Notes
    -----
    This implementation has a superior arithmetic intensity and hence
    better running time when the variant is IO bound, but it can suffer
    from numerical instability caused by catastrophic cancellation potentially
    introduced by the subtraction in the arithmetic expression above.
    numerical precision is needed.
    """

    cdef:
        GEMMTermComputer gemm_term_computer
        const DTYPE_t[::1] X_norm_squared
        const DTYPE_t[::1] Y_norm_squared

        bint use_squared_distances

    @classmethod
    def is_usable_for(cls, X, Y, metric) -> bool:
        return (PairwiseDistancesRadiusNeighborhood.is_usable_for(X, Y, metric)
                and not _in_unstable_openblas_configuration())

    def __init__(
        self,
        X,
        Y,
        DTYPE_t radius,
        bint use_squared_distances=False,
        chunk_size=None,
        strategy=None,
        sort_results=False,
        metric_kwargs=None,
    ):
        if (
            metric_kwargs is not None and
            len(metric_kwargs) > 0 and
            "Y_norm_squared" not in metric_kwargs
        ):
            warnings.warn(
                f"Some metric_kwargs have been passed ({metric_kwargs}) but aren't "
                f"usable for this case ({self.__class__.__name__}) and will be ignored.",
                UserWarning,
                stacklevel=3,
            )

        super().__init__(
            # The datasets pair here is used for exact distances computations
            datasets_pair=DatasetsPair.get_for(X, Y, metric="euclidean"),
            radius=radius,
            chunk_size=chunk_size,
            strategy=strategy,
            sort_results=sort_results,
            metric_kwargs=metric_kwargs,
        )
        # X and Y are checked by the DatasetsPair implemented as a DenseDenseDatasetsPair
        cdef:
            DenseDenseDatasetsPair datasets_pair = <DenseDenseDatasetsPair> self.datasets_pair
            ITYPE_t dist_middle_terms_chunks_size = self.Y_n_samples_chunk * self.X_n_samples_chunk

        self.gemm_term_computer = GEMMTermComputer(
            datasets_pair.X,
            datasets_pair.Y,
            self.effective_n_threads,
            self.chunks_n_threads,
            dist_middle_terms_chunks_size,
        )

        if metric_kwargs is not None and "Y_norm_squared" in metric_kwargs:
            self.Y_norm_squared = metric_kwargs.pop("Y_norm_squared")
        else:
            self.Y_norm_squared = _sqeuclidean_row_norms(datasets_pair.Y, self.effective_n_threads)

        # Do not recompute norms if datasets are identical.
        self.X_norm_squared = (
            self.Y_norm_squared if X is Y else
            _sqeuclidean_row_norms(datasets_pair.X, self.effective_n_threads)
        )
        self.use_squared_distances = use_squared_distances

        if use_squared_distances:
            # In this specialisation and this setup, the value passed to the radius is
            # already considered to be the adapted radius, so we overwrite it.
            self.r_radius = radius

    @final
    cdef void compute_exact_distances(self) nogil:
        if not self.use_squared_distances:
            PairwiseDistancesRadiusNeighborhood.compute_exact_distances(self)

    @final
    cdef void _parallel_on_X_parallel_init(
        self,
        ITYPE_t thread_num,
    ) nogil:
        PairwiseDistancesRadiusNeighborhood._parallel_on_X_parallel_init(self, thread_num)
        self.gemm_term_computer._parallel_on_X_parallel_init(thread_num)

    @final
    cdef void _parallel_on_Y_init(
        self,
    ) nogil:
        cdef ITYPE_t thread_num
        PairwiseDistancesRadiusNeighborhood._parallel_on_Y_init(self)
        self.gemm_term_computer._parallel_on_Y_init()

    @final
    cdef void _compute_and_reduce_distances_on_chunks(
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
            ITYPE_t n_X = X_end - X_start
            ITYPE_t n_Y = Y_end - Y_start
            DTYPE_t *dist_middle_terms = self.gemm_term_computer._compute_distances_on_chunks(
                X_start, X_end, Y_start, Y_end, thread_num
            )

        # Pushing the distance and their associated indices in vectors.
        for i in range(n_X):
            for j in range(n_Y):
                # Using the squared euclidean distance as the rank-preserving distance:
                #
                #             ||X_c_i||² - 2 X_c_i.Y_c_j^T + ||Y_c_j||²
                #
                squared_dist_i_j = (
                    self.X_norm_squared[i + X_start]
                    + dist_middle_terms[i * n_Y + j]
                    + self.Y_norm_squared[j + Y_start]
                )
                if squared_dist_i_j <= self.r_radius:
                    deref(self.neigh_distances_chunks[thread_num])[i + X_start].push_back(squared_dist_i_j)
                    deref(self.neigh_indices_chunks[thread_num])[i + X_start].push_back(j + Y_start)
