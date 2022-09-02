cimport numpy as cnp

from cython cimport final

from ._base cimport (
    BaseDistanceReducer64,
    _sqeuclidean_row_norms64,
)

from ._datasets_pair cimport (
    DatasetsPair64,
    DenseDenseDatasetsPair64,
)

from ._gemm_term_computer cimport GEMMTermComputer64

from ...utils._typedefs cimport ITYPE_t, DTYPE_t

import numpy as np
import warnings

from scipy.sparse import issparse
from sklearn.utils import _in_unstable_openblas_configuration
from sklearn.utils.fixes import threadpool_limits, sp_version, parse_version
from ...utils._typedefs import ITYPE, DTYPE

cnp.import_array()


def _precompute_metric_params(X, Y, metric=None, **kwds):
    """Precompute data-derived metric parameters if not provided."""
    if metric == "seuclidean" and "V" not in kwds:
        # There is a bug in scipy < 1.5 that will cause a crash if
        # X.dtype != np.double (float64). See PR #15730
        dtype = np.float64 if sp_version < parse_version("1.5") else None
        if X is Y:
            V = np.var(X, axis=0, ddof=1, dtype=dtype)
        else:
            raise ValueError(
                "The 'V' parameter is required for the seuclidean metric "
                "when Y is passed."
            )
        return {"V": V}
    if metric == "mahalanobis" and "VI" not in kwds:
        if X is Y:
            VI = np.linalg.inv(np.cov(X.T)).T
        else:
            raise ValueError(
                "The 'VI' parameter is required for the mahalanobis metric "
                "when Y is passed."
            )
        return {"VI": VI}
    return {}


cdef class PairwiseDistances64(BaseDistanceReducer64):
    """64bit implementation of PairwiseDistances."""

    @classmethod
    def compute(
        cls,
        X,
        Y,
        str metric="euclidean",
        chunk_size=None,
        dict metric_kwargs=None,
        str strategy=None,
    ):
        """Compute the pairwise-distances matrix.

        This classmethod is responsible for introspecting the arguments
        values to dispatch to the most appropriate implementation of
        :class:`PairwiseDistances64`.

        This allows decoupling the API entirely from the implementation details
        whilst maintaining RAII: all temporarily allocated datastructures necessary
        for the concrete implementation are therefore freed when this classmethod
        returns.

        No instance should directly be created outside of this class method.
        """
        if (
            metric in ("euclidean", "l2", "sqeuclidean")
            and not issparse(X)
            and not issparse(Y)
        ):
            # Specialized implementation with improved arithmetic intensity
            # and vector instructions (SIMD) by processing several vectors
            # at time to leverage a call to the BLAS GEMM routine as explained
            # in more details in the docstring.
            use_squared_distances = metric == "sqeuclidean"
            pdr = EuclideanPairwiseDistances64(
                X=X, Y=Y,
                use_squared_distances=use_squared_distances,
                chunk_size=chunk_size,
                metric_kwargs=metric_kwargs,
                strategy=strategy,
            )
        else:
            # Precompute data-derived distance metric parameters
            params = _precompute_metric_params(X, Y, metric=metric, **metric_kwargs)
            metric_kwargs.update(**params)

            # Fall back on a generic implementation that handles most scipy
            # metrics by computing the distances between 2 vectors at a time.
            pdr = PairwiseDistances64(
                datasets_pair=DatasetsPair64.get_for(X, Y, metric, metric_kwargs),
                chunk_size=chunk_size,
                metric_kwargs=metric_kwargs,
                strategy=strategy,
            )

        # Limit the number of threads in second level of nested parallelism for BLAS
        # to avoid threads over-subscription (in GEMM for instance).
        with threadpool_limits(limits=1, user_api="blas"):
            if pdr.execute_in_parallel_on_Y:
                pdr._parallel_on_Y()
            else:
                pdr._parallel_on_X()

        return pdr._finalize_results()


    def __init__(
        self,
        DatasetsPair64 datasets_pair,
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

        # Distance matrix which will be complete and returned to the caller.
        self.pairwise_distances_matrix = np.empty(
            (self.n_samples_X, self.n_samples_Y), dtype=DTYPE,
        )

    def _finalize_results(self):
        # If Y is X, then catastrophic cancellation might
        # have occurred for computations of term on the diagonal
        # which must be null. We enforce nullity of those term
        # by zeroing the diagonal.
        distance_matrix = np.asarray(self.pairwise_distances_matrix)
        if self.datasets_pair.X_is_Y:
            np.fill_diagonal(distance_matrix, 0)

        return distance_matrix

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
            DTYPE_t dist_i_j

        for i in range(X_start, X_end):
            for j in range(Y_start, Y_end):
                dist_i_j = self.datasets_pair.dist(i, j)
                self.pairwise_distances_matrix[X_start + i, Y_start + j] = dist_i_j


cdef class EuclideanPairwiseDistances64(PairwiseDistances64):
    """EuclideanDistance-specialized 64bit implementation for PairwiseDistances."""

    @classmethod
    def is_usable_for(cls, X, Y, metric) -> bool:
        return (PairwiseDistances64.is_usable_for(X, Y, metric)
                and not _in_unstable_openblas_configuration())

    def __init__(
        self,
        X,
        Y,
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
                f"usable for this case (EuclideanPairwiseDistances64) and will be ignored.",
                UserWarning,
                stacklevel=3,
            )

        super().__init__(
            # The datasets pair here is used for exact distances computations
            datasets_pair=DatasetsPair64.get_for(X, Y, metric="euclidean"),
            chunk_size=chunk_size,
            strategy=strategy,
            metric_kwargs=metric_kwargs,
        )
        # X and Y are checked by the DatasetsPair64 implemented as a DenseDenseDatasetsPair64
        cdef:
            DenseDenseDatasetsPair64 datasets_pair = <DenseDenseDatasetsPair64> self.datasets_pair
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
            self.Y_norm_squared if self.datasets_pair.X_is_Y else
            _sqeuclidean_row_norms64(datasets_pair.X, self.effective_n_threads)
        )
        self.use_squared_distances = use_squared_distances


    @final
    cdef void _parallel_on_X_parallel_init(
        self,
        ITYPE_t thread_num,
    ) nogil:
        PairwiseDistances64._parallel_on_X_parallel_init(self, thread_num)
        self.gemm_term_computer._parallel_on_X_parallel_init(thread_num)

    @final
    cdef void _parallel_on_X_init_chunk(
        self,
        ITYPE_t thread_num,
        ITYPE_t X_start,
        ITYPE_t X_end,
    ) nogil:
        PairwiseDistances64._parallel_on_X_init_chunk(self, thread_num, X_start, X_end)
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
        PairwiseDistances64._parallel_on_X_pre_compute_and_reduce_distances_on_chunks(
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
        PairwiseDistances64._parallel_on_Y_init(self)
        self.gemm_term_computer._parallel_on_Y_init()

    @final
    cdef void _parallel_on_Y_parallel_init(
        self,
        ITYPE_t thread_num,
        ITYPE_t X_start,
        ITYPE_t X_end,
    ) nogil:
        PairwiseDistances64._parallel_on_Y_parallel_init(self, thread_num, X_start, X_end)
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
        PairwiseDistances64._parallel_on_Y_pre_compute_and_reduce_distances_on_chunks(
            self,
            X_start, X_end,
            Y_start, Y_end,
            thread_num,
        )
        self.gemm_term_computer._parallel_on_Y_pre_compute_and_reduce_distances_on_chunks(
            X_start, X_end, Y_start, Y_end, thread_num
        )


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

        for i in range(n_X):
            for j in range(n_Y):
                # Using the squared euclidean distance as the rank-preserving distance:
                #
                #             ||X_c_i||² - 2 X_c_i.Y_c_j^T + ||Y_c_j||²
                #
                self.pairwise_distances_matrix[i + X_start, j + Y_start] = (
                    self.X_norm_squared[i + X_start]
                    + dist_middle_terms[i * n_Y + j]
                    + self.Y_norm_squared[j + Y_start]
                )

    def _finalize_results(self):
        distance_matrix = PairwiseDistances64._finalize_results(self)
        # Squared Euclidean distances have been used for efficiency.
        # We remap them to Euclidean distances here before finalizing
        # results.
        if not self.use_squared_distances:
            return np.sqrt(distance_matrix)
        return PairwiseDistances64._finalize_results(self)
