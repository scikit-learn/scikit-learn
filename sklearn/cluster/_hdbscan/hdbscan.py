"""
HDBSCAN: Hierarchical Density-Based Spatial Clustering
         of Applications with Noise
"""
# Author: Leland McInnes <leland.mcinnes@gmail.com>
#         Steve Astels <sastels@gmail.com>
#         John Healy <jchealy@gmail.com>
#
# License: BSD 3 clause

from numbers import Integral, Real
from pathlib import Path
from warnings import warn

import numpy as np
from joblib import Memory
from scipy.sparse import csgraph, issparse

from sklearn.base import BaseEstimator, ClusterMixin
from sklearn.metrics import pairwise_distances
from sklearn.metrics._dist_metrics import DistanceMetric
from sklearn.neighbors import BallTree, KDTree, NearestNeighbors
from sklearn.utils import check_array, gen_batches, get_chunk_n_rows
from sklearn.utils._param_validation import Interval, StrOptions, validate_params

from ._linkage import label, mst_linkage_core, mst_linkage_core_vector
from ._reachability import mutual_reachability, sparse_mutual_reachability
from ._tree import (
    compute_stability,
    condense_tree,
    get_clusters,
    labelling_at_cut,
)

FAST_METRICS = KDTree.valid_metrics + BallTree.valid_metrics
_PARAM_CONSTRAINTS = {
    "min_cluster_size": [Interval(Integral, left=2, right=None, closed="left")],
    "min_samples": [Interval(Integral, left=1, right=None, closed="left"), None],
    "cluster_selection_epsilon": [Interval(Real, left=0, right=None, closed="left")],
    "max_cluster_size": [Interval(Integral, left=0, right=None, closed="left")],
    "metric": [StrOptions(set(FAST_METRICS + ["precomputed"])), callable],
    "alpha": [Interval(Real, left=0, right=None, closed="neither")],
    "algorithm": [
        StrOptions(
            {
                "auto",
                "brute",
                "kdtree",
                "balltree",
            }
        )
    ],
    "leaf_size": [Interval(Integral, left=1, right=None, closed="left")],
    "memory": [str, None, Path],
    "n_jobs": [int],
    "cluster_selection_method": [StrOptions({"eom", "leaf"})],
    "allow_single_cluster": ["boolean"],
    "metric_params": [dict, None],
}


def _tree_to_labels(
    single_linkage_tree,
    min_cluster_size=10,
    cluster_selection_method="eom",
    allow_single_cluster=False,
    cluster_selection_epsilon=0.0,
    max_cluster_size=0,
):
    """Converts a pretrained tree and cluster size into a
    set of labels and probabilities.
    """
    condensed_tree = condense_tree(single_linkage_tree, min_cluster_size)
    stability_dict = compute_stability(condensed_tree)
    labels, probabilities = get_clusters(
        condensed_tree,
        stability_dict,
        cluster_selection_method,
        allow_single_cluster,
        cluster_selection_epsilon,
        max_cluster_size,
    )

    return (labels, probabilities, single_linkage_tree)


def _process_mst(min_spanning_tree):
    # Sort edges of the min_spanning_tree by weight
    row_order = np.argsort(min_spanning_tree.T[2])
    min_spanning_tree = min_spanning_tree[row_order, :]
    # Convert edge list into standard hierarchical clustering format
    return label(min_spanning_tree)


def _hdbscan_brute(
    X,
    min_samples=5,
    alpha=1.0,
    metric="euclidean",
    **metric_params,
):
    if metric == "precomputed":
        # Treating this case explicitly, instead of letting
        # sklearn.metrics.pairwise_distances handle it,
        # enables the usage of numpy.inf in the distance
        # matrix to indicate missing distance information.
        distance_matrix = X
    else:
        distance_matrix = pairwise_distances(X, metric=metric, **metric_params)

    if issparse(distance_matrix):
        return _hdbscan_sparse_distance_matrix(
            distance_matrix,
            min_samples,
            alpha,
            **metric_params,
        )

    mutual_reachability_ = mutual_reachability(distance_matrix, min_samples, alpha)

    min_spanning_tree = mst_linkage_core(mutual_reachability_)

    # Warn if the MST couldn't be constructed around the missing distances
    if np.isinf(min_spanning_tree.T[2]).any():
        warn(
            "The minimum spanning tree contains edge weights with value "
            "infinity. Potentially, you are missing too many distances "
            "in the initial distance matrix for the given neighborhood "
            "size.",
            UserWarning,
        )

    return _process_mst(min_spanning_tree)


def _hdbscan_sparse_distance_matrix(
    X,
    min_samples=5,
    alpha=1.0,
    **metric_params,
):
    assert issparse(X)

    # Compute sparse mutual reachability graph
    # if max_dist > 0, max distance to use when the reachability is infinite
    max_dist = metric_params.get("max_dist", 0.0)
    mutual_reachability_ = sparse_mutual_reachability(
        X.tolil(), min_points=min_samples, max_dist=max_dist, alpha=alpha
    )
    # Check connected component on mutual reachability
    # If more than one component, it means that even if the distance matrix X
    # has one component, there exists with less than `min_samples` neighbors
    if (
        csgraph.connected_components(
            mutual_reachability_, directed=False, return_labels=False
        )
        > 1
    ):
        raise ValueError(
            "There exists points with less than %s neighbors. "
            "Ensure your distance matrix has non zeros values for "
            "at least `min_sample`=%s neighbors for each points (i.e. K-nn graph), "
            "or specify a `max_dist` to use when distances are missing."
            % (min_samples, min_samples)
        )

    # Compute the minimum spanning tree for the sparse graph
    sparse_min_spanning_tree = csgraph.minimum_spanning_tree(mutual_reachability_)

    # Convert the graph to scipy cluster array format
    nonzeros = sparse_min_spanning_tree.nonzero()
    nonzero_vals = sparse_min_spanning_tree[nonzeros]
    min_spanning_tree = np.vstack(nonzeros + (nonzero_vals,)).T

    # Sort edges of the min_spanning_tree by weight
    min_spanning_tree = min_spanning_tree[np.argsort(min_spanning_tree.T[2]), :][0]

    # Convert edge list into standard hierarchical clustering format
    single_linkage_tree = label(min_spanning_tree)

    return single_linkage_tree


def _hdbscan_prims(
    X,
    algo,
    min_samples=5,
    alpha=1.0,
    metric="euclidean",
    leaf_size=40,
    n_jobs=4,
    **metric_params,
):
    # The Cython routines used require contiguous arrays
    if not X.flags["C_CONTIGUOUS"]:
        X = np.array(X, order="C")

    # Get distance to kth nearest neighbour
    nbrs = NearestNeighbors(
        n_neighbors=min_samples,
        algorithm=algo,
        leaf_size=leaf_size,
        metric=metric,
        metric_params=metric_params,
        n_jobs=n_jobs,
        p=None,
    ).fit(X)

    n_samples = X.shape[0]
    core_distances = np.empty(n_samples)
    core_distances.fill(np.nan)

    chunk_n_rows = get_chunk_n_rows(row_bytes=16 * min_samples, max_n_rows=n_samples)
    slices = gen_batches(n_samples, chunk_n_rows)
    for sl in slices:
        core_distances[sl] = nbrs.kneighbors(X[sl], min_samples)[0][:, -1]

    dist_metric = DistanceMetric.get_metric(metric, **metric_params)

    # Mutual reachability distance is implicit in mst_linkage_core_vector
    min_spanning_tree = mst_linkage_core_vector(X, core_distances, dist_metric, alpha)

    return _process_mst(min_spanning_tree)


def remap_single_linkage_tree(tree, internal_to_raw, outliers):
    """
    Takes an internal single_linkage_tree structure and adds back in a set of points
    that were initially detected as non-finite and returns that new tree.
    These points will all be merged into the final node at np.inf distance and
    considered noise points.

    Parameters
    ----------
    tree: single_linkage_tree
    internal_to_raw: dict
        a mapping from internal integer index to the raw integer index
    finite_index: ndarray
        Boolean array of which entries in the raw data were finite
    """
    finite_count = len(internal_to_raw)

    outlier_count = len(outliers)
    for i, (left, right, distance, size) in enumerate(tree):
        if left < finite_count:
            tree[i, 0] = internal_to_raw[left]
        else:
            tree[i, 0] = left + outlier_count
        if right < finite_count:
            tree[i, 1] = internal_to_raw[right]
        else:
            tree[i, 1] = right + outlier_count

    outlier_tree = np.zeros((len(outliers), 4))
    last_cluster_id = tree[tree.shape[0] - 1][0:2].max()
    last_cluster_size = tree[tree.shape[0] - 1][3]
    for i, outlier in enumerate(outliers):
        outlier_tree[i] = (outlier, last_cluster_id + 1, np.inf, last_cluster_size + 1)
        last_cluster_id += 1
        last_cluster_size += 1
    tree = np.vstack([tree, outlier_tree])
    return tree


def get_finite_row_indices(matrix):
    """
    Returns the indices of the purely finite rows of a
    sparse matrix or dense ndarray
    """
    if issparse(matrix):
        row_indices = np.array(
            [i for i, row in enumerate(matrix.tolil().data) if np.all(np.isfinite(row))]
        )
    else:
        row_indices = np.where(np.isfinite(matrix).sum(axis=1) == matrix.shape[1])[0]
    return row_indices


@validate_params(
    {
        **_PARAM_CONSTRAINTS,
        "X": ["array-like", "sparse matrix"],
    }
)
def hdbscan(
    X,
    min_cluster_size=5,
    min_samples=None,
    alpha=1.0,
    cluster_selection_epsilon=0.0,
    max_cluster_size=0,
    metric="euclidean",
    leaf_size=40,
    algorithm="auto",
    memory=None,
    n_jobs=4,
    cluster_selection_method="eom",
    allow_single_cluster=False,
    metric_params=None,
):
    """Perform HDBSCAN clustering from a vector array or distance matrix.

    Parameters
    ----------
    X : array or sparse (CSR) matrix of shape (n_samples, n_features), or \
            array of shape (n_samples, n_samples)
        A feature array, or array of distances between samples if
        `metric='precomputed'`.

    min_cluster_size : int, default=5
        The minimum number of samples in a group for that group to be
        considered a cluster; groupings smaller than this size will be left
        as noise.

    min_samples : int, default=None
        The number of samples in a neighborhood for a point
        to be considered as a core point. This includes the point itself.
        defaults to the `min_cluster_size`.

    alpha : float, default=1.0
        A distance scaling parameter as used in robust single linkage.
        See [2]_ for more information.

    cluster_selection_epsilon : float, default=0.0
        A distance threshold. Clusters below this value will be merged.
        See [3]_ for more information.

    max_cluster_size : int, default=0
        A limit to the size of clusters returned by the `eom` cluster selection
        algorithm. Has no effect if `cluster_selection_method=leaf`. Can be
        overridden in rare cases by a high value for
        `cluster_selection_epsilon`.

    metric : str or callable, default='minkowski'
        The metric to use when calculating distance between instances in a
        feature array.

        - If metric is a string or callable, it must be one of
          the options allowed by :func:`metrics.pairwise.pairwise_distances`
          for its metric parameter.

        - If metric is "precomputed", `X` is assumed to be a distance matrix and
          must be square.

    leaf_size : int, default=40
        Leaf size for trees responsible for fast nearest neighbour queries. A
        large dataset size and small leaf_size may induce excessive memory
        usage. If you are running out of memory consider increasing the
        `leaf_size` parameter.

    algorithm : str, default='auto'
        Exactly which algorithm to use; hdbscan has variants specialised
        for different characteristics of the data. By default this is set
        to `'auto'` which attempts to use a `KDTree` tree if possible,
        otherwise it uses a `BallTree` tree.

        If `X` is sparse or `metric` is invalid for both `KDTree` and
        `BallTree`, then it resolves to use the `brute` algorithm.

        Available algorithms:
        - `'brute'`
        - `'kdtree'`
        - `'balltree'`

    memory : str, default=None
        Used to cache the output of the computation of the tree.
        By default, no caching is done. If a string is given, it is the
        path to the caching directory.

    n_jobs : int, default=4
        Number of parallel jobs to run in core distance computations (if
        supported by the specific algorithm). For `n_jobs<0`,
        `(n_cpus + n_jobs + 1)` are used.

    cluster_selection_method : str, default='eom'
        The method used to select clusters from the condensed tree. The
        standard approach for HDBSCAN* is to use an Excess of Mass algorithm
        to find the most persistent clusters. Alternatively you can instead
        select the clusters at the leaves of the tree -- this provides the
        most fine grained and homogeneous clusters. Options are:
        - `eom`
        - `leaf`

    allow_single_cluster : bool, default=False
        By default HDBSCAN* will not produce a single cluster. Setting this to
        `True` will allow single cluster results in the case that you feel this
        is a valid result for your dataset.

    metric_params : dict, default=None
        Arguments passed to the distance metric.

    Returns
    -------
    labels : ndarray, shape (n_samples, )
        Cluster labels for each point.  Noisy samples are given the label -1.

    probabilities : ndarray, shape (n_samples, )
        Cluster membership strengths for each point. Noisy samples are assigned
        0.

    single_linkage_tree : ndarray, shape (n_samples - 1, 4)
        The single linkage tree produced during clustering in scipy
        hierarchical clustering format
        (see http://docs.scipy.org/doc/scipy/reference/cluster.hierarchy.html).

    References
    ----------

    .. [1] Campello, R. J., Moulavi, D., & Sander, J. (2013, April).
       Density-based clustering based on hierarchical density estimates.
       In Pacific-Asia Conference on Knowledge Discovery and Data Mining
       (pp. 160-172). Springer Berlin Heidelberg.

    .. [2] Chaudhuri, K., & Dasgupta, S. (2010). Rates of convergence for the
       cluster tree. In Advances in Neural Information Processing Systems
       (pp. 343-351).

    .. [3] Malzer, C., & Baum, M. (2019). A Hybrid Approach To Hierarchical
       Density-based Cluster Selection. arxiv preprint 1911.02282.
    """
    if min_samples is None:
        min_samples = min_cluster_size

    # Checks input and converts to an nd-array where possible
    if metric != "precomputed" or issparse(X):
        X = check_array(X, accept_sparse="csr", force_all_finite=False)
    elif isinstance(X, np.ndarray):
        # Only non-sparse, precomputed distance matrices are handled here
        # and thereby allowed to contain numpy.inf for missing distances

        # Perform check_array(X) after removing infinite values (numpy.inf)
        # from the given distance matrix.
        tmp = X.copy()
        tmp[np.isinf(tmp)] = 1
        check_array(tmp)

    memory = Memory(location=memory, verbose=0)

    metric_params = metric_params or {}
    func = None
    kwargs = dict(
        X=X,
        algo="kd_tree",
        min_samples=min_samples,
        alpha=alpha,
        metric=metric,
        leaf_size=leaf_size,
        n_jobs=n_jobs,
        **metric_params,
    )
    if "kdtree" in algorithm and metric not in KDTree.valid_metrics:
        raise ValueError(
            f"{metric} is not a valid metric for a KDTree-based algorithm. Please"
            " select a different metric."
        )
    elif "balltree" in algorithm and metric not in BallTree.valid_metrics:
        raise ValueError(
            f"{metric} is not a valid metric for a BallTree-based algorithm. Please"
            " select a different metric."
        )

    if algorithm != "auto":
        if metric != "precomputed" and issparse(X) and algorithm != "brute":
            raise ValueError("Sparse data matrices only support algorithm `brute`.")

        if algorithm == "brute":
            func = _hdbscan_brute
            for key in ("algo", "leaf_size", "n_jobs"):
                kwargs.pop(key, None)
        elif algorithm == "kdtree":
            func = _hdbscan_prims
        elif algorithm == "balltree":
            func = _hdbscan_prims
            kwargs["algo"] = "ball_tree"
    else:
        if issparse(X) or metric not in FAST_METRICS:
            # We can't do much with sparse matrices ...
            func = _hdbscan_brute
            for key in ("algo", "leaf_size", "n_jobs"):
                kwargs.pop(key, None)
        elif metric in KDTree.valid_metrics:
            func = _hdbscan_prims
        else:  # Metric is a valid BallTree metric
            func = _hdbscan_prims
            kwargs["algo"] = "ball_tree"

    single_linkage_tree = memory.cache(func)(**kwargs)

    return _tree_to_labels(
        single_linkage_tree,
        min_cluster_size,
        cluster_selection_method,
        allow_single_cluster,
        cluster_selection_epsilon,
        max_cluster_size,
    )


class HDBSCAN(ClusterMixin, BaseEstimator):
    """Perform HDBSCAN clustering from vector array or distance matrix.

    HDBSCAN - Hierarchical Density-Based Spatial Clustering of Applications
    with Noise. Performs DBSCAN over varying epsilon values and integrates
    the result to find a clustering that gives the best stability over epsilon.
    This allows HDBSCAN to find clusters of varying densities (unlike DBSCAN),
    and be more robust to parameter selection.

    Parameters
    ----------
    min_cluster_size : int, default=5
        The minimum number of samples in a group for that group to be
        considered a cluster; groupings smaller than this size will be left
        as noise.

    min_samples : int, default=None
        The number of samples in a neighborhood for a point
        to be considered as a core point. This includes the point itself.
        defaults to the `min_cluster_size`.

    cluster_selection_epsilon : float, default=0.0
        A distance threshold. Clusters below this value will be merged.
        See [5]_ for more information.

    max_cluster_size : int, default=0
        A limit to the size of clusters returned by the `eom` cluster selection
        algorithm. Has no effect if `cluster_selection_method=leaf`. Can be
        overridden in rare cases by a high value for
        `cluster_selection_epsilon`.

    metric : str or callable, default='euclidean'
        The metric to use when calculating distance between instances in a
        feature array.

        - If metric is a string or callable, it must be one of
          the options allowed by `metrics.pairwise.pairwise_distances` for its
          metric parameter.

        - If metric is "precomputed", X is assumed to be a distance matrix and
          must be square.

    alpha : float, default=1.0
        A distance scaling parameter as used in robust single linkage.
        See [3]_ for more information.

    algorithm : str, default='auto'
        Exactly which algorithm to use; hdbscan has variants specialised
        for different characteristics of the data. By default this is set
        to `'auto'` which attempts to use a `KDTree` tree if possible,
        otherwise it uses a `BallTree` tree.

        If the `X` passed during `fit` is sparse or `metric` is invalid for
        both `KDTree` and `BallTree`, then it resolves to use the `brute`
        algorithm.

        Available algorithms:
        - `'brute'`
        - `'kdtree'`
        - `'balltree'`

    leaf_size : int, default=40
        Leaf size for trees responsible for fast nearest neighbour queries. A
        large dataset size and small leaf_size may induce excessive memory
        usage. If you are running out of memory consider increasing the
        `leaf_size` parameter. Ignored for `algorithm=brute`.

    memory : str, default=None
        Used to cache the output of the computation of the tree.
        By default, no caching is done. If a string is given, it is the
        path to the caching directory.

    n_jobs : int, default=4
        Number of parallel jobs to run in core distance computations (if
        supported by the specific algorithm). For `n_jobs<0`,
        `(n_cpus + n_jobs + 1)` are used.

    cluster_selection_method : str, default='eom'
        The method used to select clusters from the condensed tree. The
        standard approach for HDBSCAN* is to use an Excess of Mass algorithm
        to find the most persistent clusters. Alternatively you can instead
        select the clusters at the leaves of the tree -- this provides the
        most fine grained and homogeneous clusters. Options are:
        - `eom`
        - `leaf`

    allow_single_cluster : bool, default=False
        By default HDBSCAN* will not produce a single cluster, setting this
        to True will override this and allow single cluster results in
        the case that you feel this is a valid result for your dataset.

    metric_params : dict, default=None
        Arguments passed to the distance metric.

    Attributes
    ----------
    labels_ : ndarray, shape (n_samples, )
        Cluster labels for each point in the dataset given to fit().
        Noisy samples are given the label -1.

    probabilities_ : ndarray, shape (n_samples, )
        The strength with which each sample is a member of its assigned
        cluster. Noise points have probability zero; points in clusters
        have values assigned proportional to the degree that they
        persist as part of the cluster.

    n_features_in_ : int
        Number of features seen during :term:`fit`.

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Defined only when `X`
        has feature names that are all strings.

    See Also
    --------
    DBSCAN : Density-Based Spatial Clustering of Applications
        with Noise.
    OPTICS : Ordering Points To Identify the Clustering Structure.
    BIRCH : Memory-efficient, online-learning algorithm.

    References
    ----------

    .. [1] Campello, R. J., Moulavi, D., & Sander, J. (2013, April).
       Density-based clustering based on hierarchical density estimates.
       In Pacific-Asia Conference on Knowledge Discovery and Data Mining
       (pp. 160-172). Springer Berlin Heidelberg.

    .. [2] Campello, R. J., Moulavi, D., Zimek, A., & Sander, J. (2015).
       Hierarchical density estimates for data clustering, visualization,
       and outlier detection. ACM Transactions on Knowledge Discovery
       from Data (TKDD), 10(1), 5.

    .. [3] Chaudhuri, K., & Dasgupta, S. (2010). Rates of convergence for the
       cluster tree. In Advances in Neural Information Processing Systems
       (pp. 343-351).

    .. [4] Moulavi, D., Jaskowiak, P.A., Campello, R.J., Zimek, A. and
       Sander, J., 2014. Density-Based Clustering Validation. In SDM
       (pp. 839-847).

    .. [5] Malzer, C., & Baum, M. (2019). A Hybrid Approach To Hierarchical
           Density-based Cluster Selection. arxiv preprint 1911.02282.

    Examples
    --------
    >>> from sklearn.cluster import HDBSCAN
    >>> from sklearn.datasets import load_digits
    >>> X, _ = load_digits(return_X_y=True)
    >>> hdb = HDBSCAN(min_cluster_size=20)
    >>> hdb.fit(X)
    HDBSCAN(min_cluster_size=20)
    >>> hdb.labels_
    array([ 2,  6, -1, ..., -1, -1, -1])
    """

    _parameter_constraints = _PARAM_CONSTRAINTS

    def __init__(
        self,
        min_cluster_size=5,
        min_samples=None,
        cluster_selection_epsilon=0.0,
        max_cluster_size=0,
        metric="euclidean",
        alpha=1.0,
        algorithm="auto",
        leaf_size=40,
        memory=None,
        n_jobs=4,
        cluster_selection_method="eom",
        allow_single_cluster=False,
        metric_params=None,
    ):
        self.min_cluster_size = min_cluster_size
        self.min_samples = min_samples
        self.alpha = alpha
        self.max_cluster_size = max_cluster_size
        self.cluster_selection_epsilon = cluster_selection_epsilon
        self.metric = metric
        self.algorithm = algorithm
        self.leaf_size = leaf_size
        self.memory = memory
        self.n_jobs = n_jobs
        self.cluster_selection_method = cluster_selection_method
        self.allow_single_cluster = allow_single_cluster
        self.metric_params = metric_params

    def fit(self, X, y=None):
        """Perform HDBSCAN clustering from features or distance matrix.

        Parameters
        ----------
        X : array or sparse (CSR) matrix of shape (n_samples, n_features), or \
                array of shape (n_samples, n_samples)
            A feature array, or array of distances between samples if
            `metric='precomputed'`.

        y : Ignored
            Ignored.

        Returns
        -------
        self : object
            Returns self.
        """
        self._validate_params()
        metric_params = self.metric_params or {}
        if self.metric != "precomputed":
            # Non-precomputed matrices may contain non-finite values.
            # Rows with these values
            X = self._validate_data(X, force_all_finite=False, accept_sparse="csr")
            self._raw_data = X

            self._all_finite = (
                np.all(np.isfinite(X.data)) if issparse(X) else np.all(np.isfinite(X))
            )

            if not self._all_finite:
                # Pass only the purely finite indices into hdbscan
                # We will later assign all non-finite points to the
                # background-1 cluster
                finite_index = get_finite_row_indices(X)
                X = X[finite_index]
                internal_to_raw = {x: y for x, y in enumerate(finite_index)}
                outliers = list(set(range(X.shape[0])) - set(finite_index))
        elif issparse(X):
            # Handle sparse precomputed distance matrices separately
            X = self._validate_data(X, accept_sparse="csr")
        else:
            # Only non-sparse, precomputed distance matrices are allowed
            #   to have numpy.inf values indicating missing distances
            X = self._validate_data(X, force_all_finite="allow-nan")

        self.n_features_in_ = X.shape[1]
        kwargs = self.get_params()
        # prediction data only applies to the persistent model, so remove
        # it from the keyword args we pass on the the function
        kwargs["metric_params"] = metric_params

        (
            self.labels_,
            self.probabilities_,
            self._single_linkage_tree_,
        ) = hdbscan(X, **kwargs)

        if self.metric != "precomputed" and not self._all_finite:
            # remap indices to align with original data in the case of
            # non-finite entries.
            self._single_linkage_tree_ = remap_single_linkage_tree(
                self._single_linkage_tree_, internal_to_raw, outliers
            )
            new_labels = np.full(self._raw_data.shape[0], -1)
            new_labels[finite_index] = self.labels_
            self.labels_ = new_labels

            new_probabilities = np.zeros(self._raw_data.shape[0])
            new_probabilities[finite_index] = self.probabilities_
            self.probabilities_ = new_probabilities

        return self

    def fit_predict(self, X, y=None):
        """Perform clustering on X and return cluster labels.

        Parameters
        ----------
        X : array or sparse (CSR) matrix of shape (n_samples, n_features), or \
                array of shape (n_samples, n_samples)
            A feature array, or array of distances between samples if
            `metric='precomputed'`.

        y : Ignored
            Ignored.

        Returns
        -------
        y : ndarray, shape (n_samples, )
            Cluster labels.
        """
        self.fit(X)
        return self.labels_

    def weighted_cluster_centroid(self, cluster_id):
        """
        Provide an approximate representative point for a given cluster.

        Note that this technique assumes a euclidean metric for speed of
        computation. For more general metrics use the `weighted_cluster_medoid`
        method which is slower, but can work with more general metrics.

        Parameters
        ----------
        cluster_id : int
            The id of the cluster to compute a centroid for.

        Returns
        -------
        centroid : array of shape (n_features,)
            A representative centroid for cluster `cluster_id`.
        """
        if not hasattr(self, "labels_"):
            raise AttributeError("Model has not been fit to data")

        if cluster_id == -1:
            raise ValueError(
                "Cannot calculate weighted centroid for -1 cluster "
                "since it is a noise cluster"
            )

        mask = self.labels_ == cluster_id
        cluster_data = self._raw_data[mask]
        cluster_membership_strengths = self.probabilities_[mask]

        return np.average(cluster_data, weights=cluster_membership_strengths, axis=0)

    def weighted_cluster_medoid(self, cluster_id):
        """
        Provide an approximate representative point for a given cluster.

        Note that this technique can be very slow and memory intensive for
        large clusters. For faster results use the `weighted_cluster_centroid`
        method which is faster, but assumes a euclidean metric.

        Parameters
        ----------
        cluster_id : int
            The id of the cluster to compute a medoid for.

        Returns
        -------
        centroid : array of shape (n_features,)
            A representative medoid for cluster `cluster_id`.
        """
        if not hasattr(self, "labels_"):
            raise AttributeError("Model has not been fit to data")

        if cluster_id == -1:
            raise ValueError(
                "Cannot calculate weighted centroid for -1 cluster "
                "since it is a noise cluster"
            )

        mask = self.labels_ == cluster_id
        cluster_data = self._raw_data[mask]
        cluster_membership_strengths = self.probabilities_[mask]
        metric_params = self.metric_params or {}

        dist_mat = pairwise_distances(cluster_data, metric=self.metric, **metric_params)

        dist_mat = dist_mat * cluster_membership_strengths
        medoid_index = np.argmin(dist_mat.sum(axis=1))
        return cluster_data[medoid_index]

    def dbscan_clustering(self, cut_distance, min_cluster_size=5):
        """
        Return clustering given by DBSCAN without border points.

        Return clustering that would be equivalent to running DBSCAN* for a
        particular cut_distance (or epsilon) DBSCAN* can be thought of as
        DBSCAN without the border points.  As such these results may differ
        slightly from `cluster.DBSCAN` due to the difference in implementation
        over the non-core points.

        This can also be thought of as a flat clustering derived from constant
        height cut through the single linkage tree.

        This represents the result of selecting a cut value for robust single linkage
        clustering. The `min_cluster_size` allows the flat clustering to declare noise
        points (and cluster smaller than `min_cluster_size`).

        Parameters
        ----------

        cut_distance : float
            The mutual reachability distance cut value to use to generate a
            flat clustering.

        min_cluster_size : int, optional
            Clusters smaller than this value with be called 'noise' and remain
            unclustered in the resulting flat clustering.

        Returns
        -------

        labels : array [n_samples]
            An array of cluster labels, one per datapoint. Unclustered points
            are assigned the label -1.
        """
        return labelling_at_cut(
            self._single_linkage_tree_, cut_distance, min_cluster_size
        )

    def _more_tags(self):
        return {"allow_nan": True}
