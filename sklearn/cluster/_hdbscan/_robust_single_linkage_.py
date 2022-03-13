# -*- coding: utf-8 -*-
"""
Robust Single Linkage: Density based single linkage clustering.
"""
import numpy as np

from sklearn.base import BaseEstimator, ClusterMixin
from sklearn.metrics import pairwise_distances
from scipy.sparse import issparse

from joblib import Memory, cpu_count
from sklearn.utils import check_array

from ._hdbscan_linkage import mst_linkage_core, mst_linkage_core_vector, label
from ._hdbscan_boruvka import KDTreeBoruvkaAlgorithm, BallTreeBoruvkaAlgorithm
from .dist_metrics import DistanceMetric
from ._hdbscan_reachability import mutual_reachability
from ._hdbscan_tree import labelling_at_cut
from sklearn.neighbors import KDTree, BallTree

# Author: Leland McInnes <leland.mcinnes@gmail.com>
#
# License: BSD 3 clause

FAST_METRICS = KDTree.valid_metrics + BallTree.valid_metrics


def _rsl_generic(X, k=5, alpha=1.4142135623730951, metric="euclidean", **kwargs):
    distance_matrix = pairwise_distances(X, metric=metric, **kwargs)

    mutual_reachability_ = mutual_reachability(distance_matrix, k)

    min_spanning_tree = mst_linkage_core(mutual_reachability_)
    min_spanning_tree = min_spanning_tree[np.argsort(min_spanning_tree.T[2]), :]

    return label(min_spanning_tree)


def _rsl_prims_kdtree(X, k=5, alpha=1.4142135623730951, metric="euclidean", **kwargs):

    # The Cython routines used require contiguous arrays
    if not X.flags["C_CONTIGUOUS"]:
        X = np.array(X, dtype=np.double, order="C")

    dim = X.shape[0]
    k = min(dim - 1, k)

    tree = KDTree(X, metric=metric, **kwargs)

    dist_metric = DistanceMetric.get_metric(metric, **kwargs)

    core_distances = tree.query(X, k=k)[0][:, -1].copy(order="C")
    min_spanning_tree = mst_linkage_core_vector(X, core_distances, dist_metric, alpha)

    return label(min_spanning_tree)


def _rsl_prims_balltree(X, k=5, alpha=1.4142135623730951, metric="euclidean", **kwargs):

    # The Cython routines used require contiguous arrays
    if not X.flags["C_CONTIGUOUS"]:
        X = np.array(X, dtype=np.double, order="C")

    dim = X.shape[0]
    k = min(dim - 1, k)

    tree = BallTree(X, metric=metric, **kwargs)

    dist_metric = DistanceMetric.get_metric(metric, **kwargs)

    core_distances = tree.query(X, k=k)[0][:, -1].copy(order="C")
    min_spanning_tree = mst_linkage_core_vector(X, core_distances, dist_metric, alpha)

    return label(min_spanning_tree)


def _rsl_boruvka_kdtree(
    X, k=5, alpha=1.0, metric="euclidean", leaf_size=40, core_dist_n_jobs=4, **kwargs
):

    if core_dist_n_jobs < 1:
        core_dist_n_jobs = max(cpu_count() + 1 + core_dist_n_jobs, 1)

    dim = X.shape[0]
    min_samples = min(dim - 1, k)

    tree = KDTree(X, metric=metric, leaf_size=leaf_size, **kwargs)
    alg = KDTreeBoruvkaAlgorithm(
        tree, min_samples, metric=metric, alpha=alpha, leaf_size=leaf_size, **kwargs
    )
    min_spanning_tree = alg.spanning_tree()

    return label(min_spanning_tree)


def _rsl_boruvka_balltree(
    X, k=5, alpha=1.0, metric="euclidean", leaf_size=40, core_dist_n_jobs=4, **kwargs
):

    if core_dist_n_jobs < 1:
        core_dist_n_jobs = max(cpu_count() + 1 + core_dist_n_jobs, 1)

    dim = X.shape[0]
    min_samples = min(dim - 1, k)

    tree = BallTree(X, metric=metric, leaf_size=leaf_size, **kwargs)
    alg = BallTreeBoruvkaAlgorithm(
        tree, min_samples, metric=metric, alpha=alpha, leaf_size=leaf_size, **kwargs
    )
    min_spanning_tree = alg.spanning_tree()

    return label(min_spanning_tree)


def robust_single_linkage(
    X,
    cut,
    k=5,
    alpha=1.4142135623730951,
    gamma=5,
    metric="euclidean",
    algorithm="best",
    memory=None,
    leaf_size=40,
    core_dist_n_jobs=4,
    metric_params=None,
):
    """
    Perform robust single linkage clustering.

    Parameters
    ----------
    X : array or sparse (CSR) matrix of shape (n_samples, n_features), or \
            array of shape (n_samples, n_samples)
        A feature array, or array of distances between samples if
        ``metric='precomputed'``.

    cut : float
        The reachability distance value to cut the cluster heirarchy at
        to derive a flat cluster labelling.

    k : int, default=5
        Reachability distances will be computed with regard to the `k`
        nearest neighbors.

    alpha : float, default=np.sqrt(2)
        Distance scaling for reachability distance computation. Reachability
        distance is computed as

        .. math::

            \\max (core_k(a), core_k(b), 1/\\alpha d(a,b)).

    gamma : int, default=5
        Ignore any clusters in the flat clustering with size less than gamma,
        and declare points in such clusters as noise points.

    metric : str or callable, default='euclidean'
        The metric to use when calculating distance between instances in a
        feature array. If metric is a string or callable, it must be one of
        the options allowed by `metrics.pairwise.pairwise_distances` for its
        metric parameter.
        If `metric="precomputed"`, X is assumed to be a distance matrix and
        must be square.

    algorithm : str, default='best'
        Exactly which algorithm to use; hdbscan has variants specialised
        for different characteristics of the data. By default this is set
        to ``best`` which chooses the "best" algorithm given the nature of
        the data. You can force other options if you believe you know
        better. Options are:
            * ``generic``
            * ``best``
            * ``prims_kdtree``
            * ``prims_balltree``
            * ``boruvka_kdtree``
            * ``boruvka_balltree``

    memory : str, default=None
        Used to cache the output of the computation of the tree.
        By default, no caching is done. If a string is given, it is the
        path to the caching directory.

    leaf_size : int, default=40
        Leaf size for trees responsible for fast nearest
        neighbour queries.

    core_dist_n_jobs : int, default=4
        Number of parallel jobs to run in core distance computations (if
        supported by the specific algorithm). For ``core_dist_n_jobs``
        below -1, (n_cpus + 1 + core_dist_n_jobs) are used.

    metric_params : dict, default=None
        Arguments passed to the distance metric.

    Returns
    -------
    labels : ndarray, shape (n_samples, )
        Cluster labels for each point.  Noisy samples are given the label -1.

    single_linkage_tree : ndarray, shape (n_samples - 1, 4)
        The single linkage tree produced during clustering in scipy
        hierarchical clustering format
        (see http://docs.scipy.org/doc/scipy/reference/cluster.hierarchy.html).

    References
    ----------
    .. [1] Chaudhuri, K., & Dasgupta, S. (2010). Rates of convergence for the
       cluster tree. In Advances in Neural Information Processing Systems
       (pp. 343-351).
    """

    if not isinstance(k, int) or k < 1:
        raise ValueError("k must be an integer greater than zero!")

    if not isinstance(alpha, float) or alpha < 1.0:
        raise ValueError("alpha must be a float greater than or equal to 1.0!")

    if not isinstance(gamma, int) or gamma < 1:
        raise ValueError("gamma must be an integer greater than zero!")

    if not isinstance(leaf_size, int) or leaf_size < 1:
        raise ValueError("Leaf size must be at least one!")

    X = check_array(X, accept_sparse="csr")
    memory = Memory(cachedir=memory, verbose=0)

    if algorithm != "best":
        if algorithm == "generic":
            single_linkage_tree = memory.cache(_rsl_generic)(
                X, k, alpha, metric, **metric_params
            )
        elif algorithm == "prims_kdtree":
            single_linkage_tree = memory.cache(_rsl_prims_kdtree)(
                X, k, alpha, metric, **metric_params
            )
        elif algorithm == "prims_balltree":
            single_linkage_tree = memory.cache(_rsl_prims_balltree)(
                X, k, alpha, metric, **metric_params
            )
        elif algorithm == "boruvka_kdtree":
            single_linkage_tree = memory.cache(_rsl_boruvka_kdtree)(
                X, k, alpha, metric, leaf_size, core_dist_n_jobs, **metric_params
            )
        elif algorithm == "boruvka_balltree":
            single_linkage_tree = memory.cache(_rsl_boruvka_balltree)(
                X, k, alpha, metric, leaf_size, core_dist_n_jobs, **metric_params
            )
        else:
            raise TypeError("Unknown algorithm type %s specified" % algorithm)
    else:
        if issparse(X) or metric not in FAST_METRICS:
            # We can't do much with sparse matrices ...
            single_linkage_tree = memory.cache(_rsl_generic)(
                X, k, alpha, metric, **metric_params
            )
        elif metric in KDTree.valid_metrics:
            # Need heuristic to decide when to go to boruvka;
            # still debugging for now
            if X.shape[1] > 128:
                single_linkage_tree = memory.cache(_rsl_prims_kdtree)(
                    X, k, alpha, metric, **metric_params
                )
            else:
                single_linkage_tree = memory.cache(_rsl_boruvka_kdtree)(
                    X, k, alpha, metric, leaf_size, core_dist_n_jobs, **metric_params
                )
        else:  # Metric is a valid BallTree metric
            # Need heuristic to decide when to go to boruvka;
            # still debugging for now
            if X.shape[1] > 128:
                single_linkage_tree = memory.cache(_rsl_prims_kdtree)(
                    X, k, alpha, metric, **metric_params
                )
            else:
                single_linkage_tree = memory.cache(_rsl_boruvka_balltree)(
                    X, k, alpha, metric, leaf_size, core_dist_n_jobs, **metric_params
                )

    labels = labelling_at_cut(single_linkage_tree, cut, gamma)

    return labels


class RobustSingleLinkage(BaseEstimator, ClusterMixin):
    r"""
    Perform robust single linkage clustering.

    Robust single linkage is a modified version of single linkage that
    attempts to be more robust to noise. Specifically the goal is to
    more accurately approximate the level set tree of the unknown
    probability density function from which the sample data has
    been drawn.

    Parameters
    ----------
    cut : float, default=0.4
        The reachability distance value to cut the cluster heirarchy at
        to derive a flat cluster labelling.

    k : int, default=5
        Reachability distances will be computed with regard to the `k`
        nearest neighbors.

    alpha : float, default=`np.sqrt(2)`
        Distance scaling for reachability distance computation. Reachability
        distance is computed as

        .. math::

            max \\{ core_k(a), core_k(b), 1/\\alpha d(a,b) \\}.

    gamma : int, default=5
        Ignore any clusters in the flat clustering with size less than gamma,
        and declare points in such clusters as noise points.

    metric : str, or callable, default='euclidean'
        The metric to use when calculating distance between instances in a
        feature array.

        If metric is a string or callable, it must be one of
        the options allowed by `metrics.pairwise.pairwise_distances` for its
        metric parameter.

        If `metric="precomputed"`, `X` is assumed to be a distance matrix and
        must be square.

    algorithm : str, default='best'
        Exactly which algorithm to use; hdbscan has variants specialised
        for different characteristics of the data. By default this is set
        to ``best`` which chooses the "best" algorithm given the nature of
        the data. You can force other options if you believe you know
        better. Options are:
            * ``small``
            * ``small_kdtree``
            * ``large_kdtree``
            * ``large_kdtree_fastcluster``

    core_dist_n_jobs : int, default=4
        Number of parallel jobs to run in core distance computations (if
        supported by the specific algorithm). For ``core_dist_n_jobs``
        below -1, (n_cpus + 1 + core_dist_n_jobs) are used.

    metric_params : dict, default=None
        Arguments passed to the distance metric.

    Attributes
    ----------
    labels_ : ndarray, shape (n_samples, )
        Cluster labels for each point.  Noisy samples are given the label -1.

    References
    ----------
    .. [1] Chaudhuri, K., & Dasgupta, S. (2010). Rates of convergence for the
       cluster tree. In Advances in Neural Information Processing Systems
       (pp. 343-351).

    See Also
    --------


    Examples
    --------
    >>> import numpy as np
    """

    def __init__(
        self,
        cut=0.4,
        k=5,
        alpha=1.4142135623730951,
        gamma=5,
        metric="euclidean",
        algorithm="best",
        core_dist_n_jobs=4,
        metric_params=None,
    ):

        self.cut = cut
        self.k = k
        self.alpha = alpha
        self.gamma = gamma
        self.metric = metric
        self.algorithm = algorithm
        self.core_dist_n_jobs = core_dist_n_jobs
        self.metric_params = metric_params

    def fit(self, X, y=None):
        """
        Perform robust single linkage clustering on `X`.

        Assumes `X` is either a series of feature vectors or a distance matrix.

        Parameters
        ----------
        X : array or sparse (CSR) matrix of shape (n_samples, n_features), or \
                array of shape (n_samples, n_samples)
            A feature array, or array of distances between samples if
            ``metric='precomputed'``.

        y : Ignored
            Ignored.

        Returns
        -------
        self : object
            Returns self.
        """
        X = check_array(X, accept_sparse="csr")

        kwargs = self.get_params()
        kwargs["metric_params"] = self.metric_params or {}

        self.labels_ = robust_single_linkage(X, **kwargs)

        return self

    def fit_predict(self, X, y=None):
        """
        Perform clustering on X and return cluster labels.

        Assumes `X` is either a series of feature vectors or a distance matrix.

        Parameters
        ----------
        X : array or sparse (CSR) matrix of shape (n_samples, n_features), or \
                array of shape (n_samples, n_samples)
            A feature array, or array of distances between samples if
            ``metric='precomputed'``.

        y : Ignored
            Ignored.

        Returns
        -------
        y : ndarray, shape (n_samples, )
            Cluster labels.
        """

        self.fit(X)
        return self.labels_
