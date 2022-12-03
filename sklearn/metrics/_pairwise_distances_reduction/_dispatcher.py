from abc import abstractmethod

import numpy as np

from typing import List

from scipy.sparse import isspmatrix_csr

from .._dist_metrics import BOOL_METRICS, METRIC_MAPPING

from ._base import _sqeuclidean_row_norms32, _sqeuclidean_row_norms64
from ._argkmin import (
    ArgKmin64,
    ArgKmin32,
)
from ._radius_neighbors import (
    RadiusNeighbors64,
    RadiusNeighbors32,
)

from ... import get_config


def sqeuclidean_row_norms(X, num_threads):
    """Compute the squared euclidean norm of the rows of X in parallel.

    Parameters
    ----------
    X : ndarray or CSR matrix of shape (n_samples, n_features)
        Input data. Must be c-contiguous.

    num_threads : int
        The number of OpenMP threads to use.

    Returns
    -------
    sqeuclidean_row_norms : ndarray of shape (n_samples,)
        Arrays containing the squared euclidean norm of each row of X.
    """
    if X.dtype == np.float64:
        return np.asarray(_sqeuclidean_row_norms64(X, num_threads))
    if X.dtype == np.float32:
        return np.asarray(_sqeuclidean_row_norms32(X, num_threads))

    raise ValueError(
        "Only float64 or float32 datasets are supported at this time, "
        f"got: X.dtype={X.dtype}."
    )


class BaseDistancesReductionDispatcher:
    """Abstract base dispatcher for pairwise distance computation & reduction.

    Each dispatcher extending the base :class:`BaseDistancesReductionDispatcher`
    dispatcher must implement the :meth:`compute` classmethod.
    """

    @classmethod
    def valid_metrics(cls) -> List[str]:
        excluded = {
            # PyFunc cannot be supported because it necessitates interacting with
            # the CPython interpreter to call user defined functions.
            "pyfunc",
            "mahalanobis",  # is numerically unstable
            # In order to support discrete distance metrics, we need to have a
            # stable simultaneous sort which preserves the order of the indices
            # because there generally is a lot of occurrences for a given values
            # of distances in this case.
            # TODO: implement a stable simultaneous_sort.
            "hamming",
            *BOOL_METRICS,
        }
        return sorted(({"sqeuclidean"} | set(METRIC_MAPPING.keys())) - excluded)

    @classmethod
    def is_usable_for(cls, X, Y, metric) -> bool:
        """Return True if the dispatcher can be used for the
        given parameters.

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
        True if the dispatcher can be used, else False.
        """

        def is_numpy_c_ordered(X):
            return hasattr(X, "flags") and X.flags.c_contiguous

        def is_valid_sparse_matrix(X):
            return (
                isspmatrix_csr(X)
                and
                # TODO: support CSR matrices without non-zeros elements
                X.nnz > 0
                and
                # TODO: support CSR matrices with int64 indices and indptr
                # See: https://github.com/scikit-learn/scikit-learn/issues/23653
                X.indices.dtype == X.indptr.dtype == np.int32
            )

        is_usable = (
            get_config().get("enable_cython_pairwise_dist", True)
            and (is_numpy_c_ordered(X) or is_valid_sparse_matrix(X))
            and (is_numpy_c_ordered(Y) or is_valid_sparse_matrix(Y))
            and X.dtype == Y.dtype
            and X.dtype in (np.float32, np.float64)
            and metric in cls.valid_metrics()
        )

        # The other joblib-based back-end might be more efficient on fused sparse-dense
        # datasets' pairs on metric="(sq)euclidean" for some configurations because it
        # uses the Squared Euclidean matrix decomposition, i.e.:
        #
        #       ||X_c_i - Y_c_j||² = ||X_c_i||² - 2 X_c_i.Y_c_j^T + ||Y_c_j||²
        #
        # calling efficient sparse-dense routines for matrix and vectors multiplication
        # implemented in SciPy we do not use yet here.
        # See: https://github.com/scikit-learn/scikit-learn/pull/23585#issuecomment-1247996669  # noqa
        # TODO: implement specialisation for (sq)euclidean on fused sparse-dense
        # using sparse-dense routines for matrix-vector multiplications.
        # Currently, only dense-dense and sparse-sparse are optimized for
        # the Euclidean case.
        fused_sparse_dense_euclidean_case_guard = not (
            (is_valid_sparse_matrix(X) ^ is_valid_sparse_matrix(Y))  # "^" is XOR
            and isinstance(metric, str)
            and "euclidean" in metric
        )

        return is_usable and fused_sparse_dense_euclidean_case_guard

    @classmethod
    @abstractmethod
    def compute(
        cls,
        X,
        Y,
        **kwargs,
    ):
        """Compute the reduction.

        Parameters
        ----------
        X : ndarray or CSR matrix of shape (n_samples_X, n_features)
            Input data.

        Y : ndarray or CSR matrix of shape (n_samples_Y, n_features)
            Input data.

        **kwargs : additional parameters for the reduction

        Notes
        -----
        This method is an abstract class method: it has to be implemented
        for all subclasses.
        """


class ArgKmin(BaseDistancesReductionDispatcher):
    """Compute the argkmin of row vectors of X on the ones of Y.

    For each row vector of X, computes the indices of k first the rows
    vectors of Y with the smallest distances.

    ArgKmin is typically used to perform
    bruteforce k-nearest neighbors queries.

    This class is not meant to be instanciated, one should only use
    its :meth:`compute` classmethod which handles allocation and
    deallocation consistently.
    """

    @classmethod
    def compute(
        cls,
        X,
        Y,
        k,
        metric="euclidean",
        chunk_size=None,
        metric_kwargs=None,
        strategy=None,
        return_distance=False,
    ):
        """Compute the argkmin reduction.

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
                embarrassingly parallel and comes with no datastructures
                synchronisation.

              - 'parallel_on_Y' dispatches chunks of Y uniformly on threads.
                Each thread processes all the chunks of X in turn. This strategy is
                a sequence of embarrassingly parallel subtasks (the inner loop on Y
                chunks) with intermediate datastructures synchronisation at each
                iteration of the sequential outer loop on X chunks.

              - 'auto' relies on a simple heuristic to choose between
                'parallel_on_X' and 'parallel_on_Y': when `X.shape[0]` is large enough,
                'parallel_on_X' is usually the most efficient strategy.
                When `X.shape[0]` is small but `Y.shape[0]` is large, 'parallel_on_Y'
                brings more opportunity for parallelism and is therefore more efficient

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
        This classmethod inspects the arguments values to dispatch to the
        dtype-specialized implementation of :class:`ArgKmin`.

        This allows decoupling the API entirely from the implementation details
        whilst maintaining RAII: all temporarily allocated datastructures necessary
        for the concrete implementation are therefore freed when this classmethod
        returns.
        """
        if X.dtype == Y.dtype == np.float64:
            return ArgKmin64.compute(
                X=X,
                Y=Y,
                k=k,
                metric=metric,
                chunk_size=chunk_size,
                metric_kwargs=metric_kwargs,
                strategy=strategy,
                return_distance=return_distance,
            )

        if X.dtype == Y.dtype == np.float32:
            return ArgKmin32.compute(
                X=X,
                Y=Y,
                k=k,
                metric=metric,
                chunk_size=chunk_size,
                metric_kwargs=metric_kwargs,
                strategy=strategy,
                return_distance=return_distance,
            )

        raise ValueError(
            "Only float64 or float32 datasets pairs are supported at this time, "
            f"got: X.dtype={X.dtype} and Y.dtype={Y.dtype}."
        )


class RadiusNeighbors(BaseDistancesReductionDispatcher):
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
        radius,
        metric="euclidean",
        chunk_size=None,
        metric_kwargs=None,
        strategy=None,
        return_distance=False,
        sort_results=False,
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
                embarrassingly parallel and comes with no datastructures
                synchronisation.

              - 'parallel_on_Y' dispatches chunks of Y uniformly on threads.
                Each thread processes all the chunks of X in turn. This strategy is
                a sequence of embarrassingly parallel subtasks (the inner loop on Y
                chunks) with intermediate datastructures synchronisation at each
                iteration of the sequential outer loop on X chunks.

              - 'auto' relies on a simple heuristic to choose between
                'parallel_on_X' and 'parallel_on_Y': when `X.shape[0]` is large enough,
                'parallel_on_X' is usually the most efficient strategy.
                When `X.shape[0]` is small but `Y.shape[0]` is large, 'parallel_on_Y'
                brings more opportunity for parallelism and is therefore more efficient
                despite the synchronization step at each iteration of the outer loop
                on chunks of `X`.

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
        This classmethod inspects the arguments values to dispatch to the
        dtype-specialized implementation of :class:`RadiusNeighbors`.

        This allows decoupling the API entirely from the implementation details
        whilst maintaining RAII: all temporarily allocated datastructures necessary
        for the concrete implementation are therefore freed when this classmethod
        returns.
        """
        if X.dtype == Y.dtype == np.float64:
            return RadiusNeighbors64.compute(
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

        if X.dtype == Y.dtype == np.float32:
            return RadiusNeighbors32.compute(
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
            "Only float64 or float32 datasets pairs are supported at this time, "
            f"got: X.dtype={X.dtype} and Y.dtype={Y.dtype}."
        )
