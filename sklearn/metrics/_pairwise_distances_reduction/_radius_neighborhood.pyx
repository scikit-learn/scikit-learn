cimport numpy as cnp
import numpy as np
import warnings

from libcpp.memory cimport shared_ptr, make_shared
from libcpp.vector cimport vector
from cython cimport final
from cython.operator cimport dereference as deref
from cython.parallel cimport parallel, prange

from ._base cimport (
    PairwiseDistancesReduction64,
    _sqeuclidean_row_norms64
)

from ._datasets_pair cimport (
    DatasetsPair,
    DenseDenseDatasetsPair,
)

from ._gemm_term_computer cimport GEMMTermComputer64

from ...utils._sorting cimport simultaneous_sort
from ...utils._typedefs cimport ITYPE_t, DTYPE_t
from ...utils._vector_sentinel cimport vector_to_nd_array

from numbers import Real
from scipy.sparse import issparse
from sklearn.utils import check_scalar, _in_unstable_openblas_configuration
from sklearn.utils.fixes import threadpool_limits

from ._base import PairwiseDistancesReduction

cnp.import_array()

# TODO: change for `libcpp.algorithm.move` once Cython 3 is used
# Introduction in Cython:
# https://github.com/cython/cython/blob/05059e2a9b89bf6738a7750b905057e5b1e3fe2e/Cython/Includes/libcpp/algorithm.pxd#L47 #noqa
cdef extern from "<algorithm>" namespace "std" nogil:
    OutputIt move[InputIt, OutputIt](InputIt first, InputIt last, OutputIt d_first) except + #noqa

######################

cdef cnp.ndarray[object, ndim=1] coerce_vectors_to_nd_arrays(
    shared_ptr[vector_vector_DITYPE_t] vecs
):
    """Coerce a std::vector of std::vector to a ndarray of ndarray."""
    cdef:
        ITYPE_t n = deref(vecs).size()
        cnp.ndarray[object, ndim=1] nd_arrays_of_nd_arrays = np.empty(n, dtype=np.ndarray)

    for i in range(n):
        nd_arrays_of_nd_arrays[i] = vector_to_nd_array(&(deref(vecs)[i]))

    return nd_arrays_of_nd_arrays

#####################

class PairwiseDistancesRadiusNeighborhood(PairwiseDistancesReduction):
    """Compute radius-based neighbors for two sets of vectors.

    For each row-vector X[i] of the queries X, find all the indices j of
    row-vectors in Y such that:

                        dist(X[i], Y[j]) <= radius

    The distance function `dist` depends on the values of the `metric`
    and `metric_kwargs` parameters.

    This class is not meant to be instanciated, one should only use
    its :meth:`compute` classmethod which handles allocation and
    deallocation consistently.
    """

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
        values to dispatch to the private dtype-specialized implementation of
        :class:`PairwiseDistancesRadiusNeighborhood`.

        All temporarily allocated datastructures necessary for the concrete
        implementation are therefore freed when this classmethod returns.

        This allows entirely decoupling the API entirely from the
        implementation details whilst maintaining RAII.
        """
        if X.dtype == Y.dtype == np.float64:
            return PairwiseDistancesRadiusNeighborhood64.compute(
                X=X,
                Y=Y,
                radius=radius,
                metric=metric,
                chunk_size=chunk_size,
                metric_kwargs=metric_kwargs,
                strategy=strategy,
                sort_results=sort_results,
                return_distance=return_distance,
            )
        raise ValueError(
            f"Only 64bit float datasets are supported at this time, "
            f"got: X.dtype={X.dtype} and Y.dtype={Y.dtype}."
        )


cdef class PairwiseDistancesRadiusNeighborhood64(PairwiseDistancesReduction64):
    """64bit implementation of PairwiseDistancesArgKmin."""

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
        """Compute the radius-neighbors reduction.

        This classmethod is responsible for introspecting the arguments
        values to dispatch to the most appropriate implementation of
        :class:`PairwiseDistancesRadiusNeighborhood64`.

        This allows decoupling the API entirely from the implementation details
        whilst maintaining RAII: all temporarily allocated datastructures necessary
        for the concrete implementation are therefore freed when this classmethod
        returns.

        No instance should directly be created outside of this class method.
        """
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
            pda = FastEuclideanPairwiseDistancesRadiusNeighborhood64(
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
            pda = PairwiseDistancesRadiusNeighborhood64(
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

    cdef void _parallel_on_X_init_chunk(
        self,
        ITYPE_t thread_num,
        ITYPE_t X_start,
        ITYPE_t X_end,
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
            for idx in prange(self.n_samples_X, schedule='static'):
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

        for i in prange(self.n_samples_X, nogil=True, schedule='static',
                        num_threads=self.effective_n_threads):
            for j in range(deref(self.neigh_indices)[i].size()):
                deref(self.neigh_distances)[i][j] = (
                        self.datasets_pair.distance_metric._rdist_to_dist(
                            # Guard against eventual -0., causing nan production.
                            max(deref(self.neigh_distances)[i][j], 0.)
                        )
                )


cdef class FastEuclideanPairwiseDistancesRadiusNeighborhood64(PairwiseDistancesRadiusNeighborhood64):
    """EuclideanDistance-specialized 64bit implementation for PairwiseDistancesRadiusNeighborhood."""

    @classmethod
    def is_usable_for(cls, X, Y, metric) -> bool:
        return (PairwiseDistancesRadiusNeighborhood64.is_usable_for(X, Y, metric)
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
                f"usable for this case (FastEuclideanPairwiseDistancesRadiusNeighborhood) and will be ignored.",
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

        self.gemm_term_computer = GEMMTermComputer64(
            datasets_pair.X,
            datasets_pair.Y,
            self.effective_n_threads,
            self.chunks_n_threads,
            dist_middle_terms_chunks_size,
            n_features=datasets_pair.X.shape[1],
            chunk_size=self.chunk_size,
        )

        if metric_kwargs is not None and "Y_norm_squared" in metric_kwargs:
            self.Y_norm_squared = metric_kwargs.pop("Y_norm_squared")
        else:
            self.Y_norm_squared = _sqeuclidean_row_norms64(datasets_pair.Y, self.effective_n_threads)

        # Do not recompute norms if datasets are identical.
        self.X_norm_squared = (
            self.Y_norm_squared if X is Y else
            _sqeuclidean_row_norms64(datasets_pair.X, self.effective_n_threads)
        )
        self.use_squared_distances = use_squared_distances

        if use_squared_distances:
            # In this specialisation and this setup, the value passed to the radius is
            # already considered to be the adapted radius, so we overwrite it.
            self.r_radius = radius

    @final
    cdef void _parallel_on_X_parallel_init(
        self,
        ITYPE_t thread_num,
    ) nogil:
        PairwiseDistancesRadiusNeighborhood64._parallel_on_X_parallel_init(self, thread_num)
        self.gemm_term_computer._parallel_on_X_parallel_init(thread_num)

    @final
    cdef void _parallel_on_X_init_chunk(
        self,
        ITYPE_t thread_num,
        ITYPE_t X_start,
        ITYPE_t X_end,
    ) nogil:
        PairwiseDistancesRadiusNeighborhood64._parallel_on_X_init_chunk(self, thread_num, X_start, X_end)
        self.gemm_term_computer._parallel_on_X_init_chunk(thread_num, X_start, X_end)

    @final
    cdef void _parallel_on_X_pre_compute_and_reduce_distances_on_chunks(
        self,
        ITYPE_t X_start,
        ITYPE_t X_end,
        ITYPE_t Y_start,
        ITYPE_t Y_end,
        ITYPE_t thread_num,
    ) nogil:
        PairwiseDistancesRadiusNeighborhood64._parallel_on_X_pre_compute_and_reduce_distances_on_chunks(
            self,
            X_start, X_end,
            Y_start, Y_end,
            thread_num,
        )
        self.gemm_term_computer._parallel_on_X_pre_compute_and_reduce_distances_on_chunks(
            X_start, X_end, Y_start, Y_end, thread_num,
        )

    @final
    cdef void _parallel_on_Y_init(
        self,
    ) nogil:
        cdef ITYPE_t thread_num
        PairwiseDistancesRadiusNeighborhood64._parallel_on_Y_init(self)
        self.gemm_term_computer._parallel_on_Y_init()

    @final
    cdef void _parallel_on_Y_parallel_init(
        self,
        ITYPE_t thread_num,
        ITYPE_t X_start,
        ITYPE_t X_end,
    ) nogil:
        PairwiseDistancesRadiusNeighborhood64._parallel_on_Y_parallel_init(self, thread_num, X_start, X_end)
        self.gemm_term_computer._parallel_on_Y_parallel_init(thread_num, X_start, X_end)

    @final
    cdef void _parallel_on_Y_pre_compute_and_reduce_distances_on_chunks(
        self,
        ITYPE_t X_start,
        ITYPE_t X_end,
        ITYPE_t Y_start,
        ITYPE_t Y_end,
        ITYPE_t thread_num,
    ) nogil:
        PairwiseDistancesRadiusNeighborhood64._parallel_on_Y_pre_compute_and_reduce_distances_on_chunks(
            self,
            X_start, X_end,
            Y_start, Y_end,
            thread_num,
        )
        self.gemm_term_computer._parallel_on_Y_pre_compute_and_reduce_distances_on_chunks(
            X_start, X_end, Y_start, Y_end, thread_num
        )

    @final
    cdef void compute_exact_distances(self) nogil:
        if not self.use_squared_distances:
            PairwiseDistancesRadiusNeighborhood64.compute_exact_distances(self)

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
