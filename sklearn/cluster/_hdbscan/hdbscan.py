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
from warnings import warn

import numpy as np
from scipy.sparse import csgraph, issparse

from sklearn.base import BaseEstimator, ClusterMixin
from sklearn.metrics import pairwise_distances
from sklearn.metrics._dist_metrics import DistanceMetric
from sklearn.neighbors import BallTree, KDTree, NearestNeighbors
from sklearn.utils._param_validation import Interval, StrOptions

from ._linkage import label, mst_linkage_core, mst_linkage_core_vector
from ._reachability import mutual_reachability, sparse_mutual_reachability
from ._tree import compute_stability, condense_tree, get_clusters, labelling_at_cut

FAST_METRICS = KDTree.valid_metrics + BallTree.valid_metrics


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
    n_jobs=None,
    **metric_params,
):
    if metric == "precomputed":
        # Treating this case explicitly, instead of letting
        # sklearn.metrics.pairwise_distances handle it,
        # enables the usage of numpy.inf in the distance
        # matrix to indicate missing distance information.
        distance_matrix = X
    else:
        distance_matrix = pairwise_distances(
            X, metric=metric, n_jobs=n_jobs, **metric_params
        )

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
            "There exists points with fewer than %s neighbors. "
            "Ensure your distance matrix has non-zero values for "
            "at least `min_sample`=%s neighbors for each points (i.e. K-nn graph), "
            "or specify a `max_dist` in `metric_params` to use when distances "
            "are missing." % (min_samples, min_samples)
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
    n_jobs=None,
    **metric_params,
):
    # The Cython routines used require contiguous arrays
    X = np.asarray(X, order="C")

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

    neighbors_distances, _ = nbrs.kneighbors(X, min_samples, return_distance=True)
    core_distances = np.ascontiguousarray(neighbors_distances[:, -1])
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


class HDBSCAN(ClusterMixin, BaseEstimator):
    """Cluster data using hierarchical density-based clustering.

    HDBSCAN - Hierarchical Density-Based Spatial Clustering of Applications
    with Noise. Performs DBSCAN over varying epsilon values and integrates
    the result to find a clustering that gives the best stability over epsilon.
    This allows HDBSCAN to find clusters of varying densities (unlike DBSCAN),
    and be more robust to parameter selection.

    ..versionadded:: 1.2

    Parameters
    ----------
    min_cluster_size : int, default=5
        The minimum number of samples in a group for that group to be
        considered a cluster; groupings smaller than this size will be left
        as noise.

    min_samples : int, default=None
        The number of samples in a neighborhood for a point
        to be considered as a core point. This includes the point itself.
        When `None`, defaults to `min_cluster_size`.

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
          the options allowed by :func:`~sklearn.metrics.pairwise.pairwise_distances` for its
          metric parameter.

        - If metric is "precomputed", X is assumed to be a distance matrix and
          must be square.

    alpha : float, default=1.0
        A distance scaling parameter as used in robust single linkage.
        See [3]_ for more information.

    algorithm : {"auto", "brute", "kdtree", "balltree"}, default="auto"
        Exactly which algorithm to use; hdbscan has variants specialised
        for different characteristics of the data. By default this is set
        to `'auto'` which attempts to use a :class:`~sklearn.neighbors.KDTree` tree if possible,
        otherwise it uses a :class:`~sklearn.neighbors.BallTree` tree.

        If the `X` passed during `fit` is sparse or `metric` is invalid for
        both :class:`~sklearn.neighbors.KDTree` and :class:`~sklearn.neighbors.BallTree`, then it resolves to use the `"brute"`
        algorithm.

        Available algorithms:
        - `'brute'`
        - `'kdtree'`
        - `'balltree'`

    leaf_size : int, default=40
        Leaf size for trees responsible for fast nearest neighbour queries when a KDTree or a BallTree are used as algorithms. A
        large dataset size and small `leaf_size` may induce excessive memory
        usage. If you are running out of memory consider increasing the
        `leaf_size` parameter. Ignored for `algorithm=brute`.

    n_jobs : int, default=None
        Number of jobs to run in parallel to calculate distances.
        `None` means 1 unless in a :obj:`joblib.parallel_backend` context.
        `-1` means using all processors. See :term:`Glossary <n_jobs>`
        for more details.

    cluster_selection_method : {"eom", "leaf"}, default="eom"
        The method used to select clusters from the condensed tree. The
        standard approach for HDBSCAN* is to use an Excess of Mass algorithm
        to find the most persistent clusters. Alternatively you can instead
        select the clusters at the leaves of the tree -- this provides the
        most fine grained and homogeneous clusters.

    allow_single_cluster : bool, default=False
        By default HDBSCAN* will not produce a single cluster, setting this
        to True will override this and allow single cluster results in
        the case that you feel this is a valid result for your dataset.

    store_centers : str, default=None
        Which, if any, cluster centers to compute and store. The options are:
        - `None` which does not compute nor store any centers.
        - `"centroid"` which calculates the center by taking the weighted
          average of their positions. Note that the algorithm uses the
          euclidean metric and does not guarantee that the output will be
          an observed data point.
        - `"medoid"` which calculates the center by taking the point in the
          fitted data which minimizes the distance to all other points in
          the cluster. This is slower than "centroid" since it requires
          computing additional pairwise distances between points of the
          same cluster but guarantees the output is an observed data point.
          The medoid is also well-defined for arbitrary metrics, and does not
          depend on a euclidean metric.
        - `"both"`which computes and stores both forms of centers.

    metric_params : dict, default=None
        Arguments passed to the distance metric.

    Attributes
    ----------
    labels_ : ndarray of shape (n_samples,)
        Cluster labels for each point in the dataset given to :term:`fit`.
        Noisy samples are given the label -1.

    probabilities_ : ndarray of shape (n_samples,)
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

    _parameter_constraints = {
        "min_cluster_size": [Interval(Integral, left=2, right=None, closed="left")],
        "min_samples": [Interval(Integral, left=1, right=None, closed="left"), None],
        "cluster_selection_epsilon": [
            Interval(Real, left=0, right=None, closed="left")
        ],
        "max_cluster_size": [Interval(Integral, left=0, right=None, closed="left")],
        "metric": [StrOptions(set(FAST_METRICS + {"precomputed"})), callable],
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
        "n_jobs": [Integral, None],
        "cluster_selection_method": [StrOptions({"eom", "leaf"})],
        "allow_single_cluster": ["boolean"],
        "store_centers": [None, StrOptions({"centroid", "medoid", "both"})],
        "metric_params": [dict, None],
    }

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
        n_jobs=4,
        cluster_selection_method="eom",
        allow_single_cluster=False,
        store_centers=None,
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
        self.n_jobs = n_jobs
        self.cluster_selection_method = cluster_selection_method
        self.allow_single_cluster = allow_single_cluster
        self.store_centers = store_centers
        self.metric_params = metric_params

    def fit(self, X, y=None):
        """Find clusters based on hierarchical density-based clustering.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features), or \
                ndarray of shape (n_samples, n_samples)
            A feature array, or array of distances between samples if
            `metric='precomputed'`.

        y : None
            Ignored.

        Returns
        -------
        self : object
            Returns self.
        """
        self._validate_params()
        self._metric_params = self.metric_params or {}
        if self.metric != "precomputed":
            # Non-precomputed matrices may contain non-finite values.
            X = self._validate_data(
                X, accept_sparse="csr", force_all_finite=False, dtype=np.float64
            )
            self._raw_data = X

            all_finite = (
                np.all(np.isfinite(X.data)) if issparse(X) else np.all(np.isfinite(X))
            )

            if not all_finite:
                # Pass only the purely finite indices into hdbscan
                # We will later assign all non-finite points to the
                # background-1 cluster
                finite_index = get_finite_row_indices(X)
                X = X[finite_index]
                internal_to_raw = {x: y for x, y in enumerate(finite_index)}
                outliers = list(set(range(X.shape[0])) - set(finite_index))
        elif issparse(X):
            # Handle sparse precomputed distance matrices separately
            X = self._validate_data(
                X,
                accept_sparse="csr",
                dtype=np.float64,
            )
        else:
            # Only non-sparse, precomputed distance matrices are handled here
            # and thereby allowed to contain numpy.inf for missing distances

            # Perform data validation after removing infinite values (numpy.inf)
            # from the given distance matrix.
            X = self._validate_data(X, force_all_finite=False, dtype=np.float64)
            if np.isnan(X).any():
                # TODO: Support np.nan in Cython implementation for sparse
                # HDBSCAN
                raise ValueError("np.nan values found in precomputed-sparse")
        self.n_features_in_ = X.shape[1]
        self._min_samples = (
            self.min_cluster_size if self.min_samples is None else self.min_samples
        )

        func = None
        kwargs = dict(
            X=X,
            algo="kd_tree",
            min_samples=self._min_samples,
            alpha=self.alpha,
            metric=self.metric,
            leaf_size=self.leaf_size,
            n_jobs=self.n_jobs,
            **self._metric_params,
        )
        if "kdtree" in self.algorithm and self.metric not in KDTree.valid_metrics:
            raise ValueError(
                f"{self.metric} is not a valid metric for a KDTree-based algorithm."
                " Please select a different metric."
            )
        elif "balltree" in self.algorithm and self.metric not in BallTree.valid_metrics:
            raise ValueError(
                f"{self.metric} is not a valid metric for a BallTree-based algorithm."
                " Please select a different metric."
            )

        if self.algorithm != "auto":
            if (
                self.metric != "precomputed"
                and issparse(X)
                and self.algorithm != "brute"
            ):
                raise ValueError("Sparse data matrices only support algorithm `brute`.")

            if self.algorithm == "brute":
                func = _hdbscan_brute
                for key in ("algo", "leaf_size"):
                    kwargs.pop(key, None)
            elif self.algorithm == "kdtree":
                func = _hdbscan_prims
            elif self.algorithm == "balltree":
                func = _hdbscan_prims
                kwargs["algo"] = "ball_tree"
        else:
            if issparse(X) or self.metric not in FAST_METRICS:
                # We can't do much with sparse matrices ...
                func = _hdbscan_brute
                for key in ("algo", "leaf_size"):
                    kwargs.pop(key, None)
            elif self.metric in KDTree.valid_metrics:
                # TODO: Benchmark KD vs Ball Tree efficacy
                func = _hdbscan_prims
            else:
                # Metric is a valid BallTree metric
                func = _hdbscan_prims
                kwargs["algo"] = "ball_tree"

        single_linkage_tree = func(**kwargs)

        (
            self.labels_,
            self.probabilities_,
            self._single_linkage_tree_,
        ) = _tree_to_labels(
            single_linkage_tree,
            self.min_cluster_size,
            self.cluster_selection_method,
            self.allow_single_cluster,
            self.cluster_selection_epsilon,
            self.max_cluster_size,
        )
        if self.metric != "precomputed" and not all_finite:
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

        if self.store_centers:
            self._weighted_cluster_center(X)
        return self

    def fit_predict(self, X, y=None):
        """Cluster X and return the associated cluster labels.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features), or \
                ndarray of shape (n_samples, n_samples)
            A feature array, or array of distances between samples if
            `metric='precomputed'`.

        y : None
            Ignored.

        Returns
        -------
        y : ndarray of shape (n_samples,)
            Cluster labels.
        """
        self.fit(X)
        return self.labels_

    def _weighted_cluster_center(self, X):
        n_clusters = len(set(self.labels_))
        mask = np.empty((X.shape[0],), dtype=np.bool_)
        make_centroids = self.store_centers in ("centroid", "both")
        make_medoids = self.store_centers in ("medoid", "both")

        if make_centroids:
            self.centroids_ = np.empty((n_clusters, X.shape[1]), dtype=np.float64)
        if make_medoids:
            self.medoids_ = np.empty((n_clusters, X.shape[1]), dtype=np.float64)

        # Need to handle iteratively seen each cluster may have a different
        # number of samples, hence we can't create a homogenous 3D array.
        for idx in range(n_clusters):
            mask = self.labels_ == idx
            data = X[mask]
            strength = self.probabilities_[mask]
            if make_centroids:
                self.centroids_[idx] = np.average(data, weights=strength, axis=0)
            if make_medoids:
                # TODO: Implement weighted argmin PWD backend
                dist_mat = pairwise_distances(
                    data, metric=self.metric, **self._metric_params
                )
                dist_mat = dist_mat * strength
                medoid_index = np.argmin(dist_mat.sum(axis=1))
                self.medoids_[idx] = data[medoid_index]
        return

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

        min_cluster_size : int, default=5
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
        return {"allow_nan": self.metric != "precomputed"}
