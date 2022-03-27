"""
HDBSCAN: Hierarchical Density-Based Spatial Clustering
         of Applications with Noise
"""
# Author: Leland McInnes <leland.mcinnes@gmail.com>
#         Steve Astels <sastels@gmail.com>
#         John Healy <jchealy@gmail.com>
#
# License: BSD 3 clause

import numpy as np

from sklearn.base import BaseEstimator, ClusterMixin
from sklearn.metrics import pairwise_distances
from scipy.sparse import issparse
from sklearn.neighbors import KDTree, BallTree
from joblib import Memory
from warnings import warn
from sklearn.utils import check_array
from joblib.parallel import cpu_count

from scipy.sparse import csgraph

from ._hdbscan_linkage import (
    mst_linkage_core,
    mst_linkage_core_vector,
    label,
)
from ._hdbscan_tree import (
    condense_tree,
    compute_stability,
    get_clusters,
    labelling_at_cut,
)
from ._hdbscan_reachability import mutual_reachability, sparse_mutual_reachability

from ._hdbscan_boruvka import KDTreeBoruvkaAlgorithm, BallTreeBoruvkaAlgorithm
from sklearn.metrics._dist_metrics import DistanceMetric

FAST_METRICS = KDTree.valid_metrics + BallTree.valid_metrics + ["cosine", "arccos"]


def _tree_to_labels(
    single_linkage_tree,
    min_cluster_size=10,
    cluster_selection_method="eom",
    allow_single_cluster=False,
    match_reference_implementation=False,
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
        match_reference_implementation,
        cluster_selection_epsilon,
        max_cluster_size,
    )

    return (labels, probabilities, single_linkage_tree)


def _hdbscan_generic(
    X,
    min_samples=5,
    alpha=1.0,
    metric="minkowski",
    p=2,
    leaf_size=None,
    **kwargs,
):
    if metric == "minkowski":
        distance_matrix = pairwise_distances(X, metric=metric, p=p)
    elif metric == "arccos":
        distance_matrix = pairwise_distances(X, metric="cosine", **kwargs)
    elif metric == "precomputed":
        # Treating this case explicitly, instead of letting
        #   sklearn.metrics.pairwise_distances handle it,
        #   enables the usage of numpy.inf in the distance
        #   matrix to indicate missing distance information.
        # TODO: Check if copying is necessary
        distance_matrix = X.copy()
    else:
        distance_matrix = pairwise_distances(X, metric=metric, **kwargs)

    if issparse(distance_matrix):
        # raise TypeError('Sparse distance matrices not yet supported')
        return _hdbscan_sparse_distance_matrix(
            distance_matrix,
            min_samples,
            alpha,
            **kwargs,
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

    # Sort edges of the min_spanning_tree by weight
    min_spanning_tree = min_spanning_tree[np.argsort(min_spanning_tree.T[2]), :]

    # Convert edge list into standard hierarchical clustering format
    single_linkage_tree = label(min_spanning_tree)

    return single_linkage_tree


def _hdbscan_sparse_distance_matrix(
    X,
    min_samples=5,
    alpha=1.0,
    **kwargs,
):
    assert issparse(X)
    # Check for connected component on X
    if csgraph.connected_components(X, directed=False, return_labels=False) > 1:
        raise ValueError(
            "Sparse distance matrix has multiple connected "
            "components!\nThat is, there exist groups of points "
            "that are completely disjoint -- there are no distance "
            "relations connecting them\n"
            "Run hdbscan on each component."
        )

    lil_matrix = X.tolil()

    # Compute sparse mutual reachability graph
    # if max_dist > 0, max distance to use when the reachability is infinite
    max_dist = kwargs.get("max_dist", 0.0)
    mutual_reachability_ = sparse_mutual_reachability(
        lil_matrix, min_points=min_samples, max_dist=max_dist, alpha=alpha
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


def _hdbscan_prims_kdtree(
    X,
    min_samples=5,
    alpha=1.0,
    metric="minkowski",
    leaf_size=40,
    **kwargs,
):
    if X.dtype != np.float64:
        X = X.astype(np.float64)

    # The Cython routines used require contiguous arrays
    if not X.flags["C_CONTIGUOUS"]:
        X = np.array(X, dtype=np.double, order="C")

    tree = KDTree(X, metric=metric, leaf_size=leaf_size, **kwargs)

    # TO DO: Deal with p for minkowski appropriately
    dist_metric = DistanceMetric.get_metric(metric, **kwargs)

    # Get distance to kth nearest neighbour
    core_distances = tree.query(
        X, k=min_samples + 1, dualtree=True, breadth_first=True
    )[0][:, -1].copy(order="C")

    # Mutual reachability distance is implicit in mst_linkage_core_vector
    min_spanning_tree = mst_linkage_core_vector(X, core_distances, dist_metric, alpha)

    # Sort edges of the min_spanning_tree by weight
    min_spanning_tree = min_spanning_tree[np.argsort(min_spanning_tree.T[2]), :]

    # Convert edge list into standard hierarchical clustering format
    single_linkage_tree = label(min_spanning_tree)

    return single_linkage_tree


def _hdbscan_prims_balltree(
    X,
    min_samples=5,
    alpha=1.0,
    metric="minkowski",
    leaf_size=40,
    **kwargs,
):
    if X.dtype != np.float64:
        X = X.astype(np.float64)

    # The Cython routines used require contiguous arrays
    if not X.flags["C_CONTIGUOUS"]:
        X = np.array(X, dtype=np.double, order="C")

    tree = BallTree(X, metric=metric, leaf_size=leaf_size, **kwargs)

    dist_metric = DistanceMetric.get_metric(metric, **kwargs)

    # Get distance to kth nearest neighbour
    core_distances = tree.query(
        X, k=min_samples + 1, dualtree=True, breadth_first=True
    )[0][:, -1].copy(order="C")

    # Mutual reachability distance is implicit in mst_linkage_core_vector
    min_spanning_tree = mst_linkage_core_vector(X, core_distances, dist_metric, alpha)
    # Sort edges of the min_spanning_tree by weight
    min_spanning_tree = min_spanning_tree[np.argsort(min_spanning_tree.T[2]), :]
    # Convert edge list into standard hierarchical clustering format
    single_linkage_tree = label(min_spanning_tree)

    return single_linkage_tree


def _hdbscan_boruvka_kdtree(
    X,
    min_samples=5,
    metric="minkowski",
    leaf_size=40,
    approx_min_span_tree=True,
    n_jobs=4,
    **kwargs,
):
    if leaf_size < 3:
        leaf_size = 3

    if n_jobs < 1:
        n_jobs = max(cpu_count() + 1 + n_jobs, 1)

    if X.dtype != np.float64:
        X = X.astype(np.float64)

    tree = KDTree(X, metric=metric, leaf_size=leaf_size, **kwargs)

    n_samples = X.shape[0]
    if min_samples + 1 > n_samples:
        raise ValueError(
            "Expected min_samples + 1 <= n_samples, "
            f" but {min_samples+1=}, {n_samples=}"
        )

    alg = KDTreeBoruvkaAlgorithm(
        tree,
        min_samples,
        metric=metric,
        leaf_size=leaf_size // 3,
        approx_min_span_tree=approx_min_span_tree,
        n_jobs=n_jobs,
        **kwargs,
    )
    min_spanning_tree = alg.spanning_tree()
    # Sort edges of the min_spanning_tree by weight
    row_order = np.argsort(min_spanning_tree.T[2])
    min_spanning_tree = min_spanning_tree[row_order, :]
    # Convert edge list into standard hierarchical clustering format
    single_linkage_tree = label(min_spanning_tree)

    return single_linkage_tree


def _hdbscan_boruvka_balltree(
    X,
    min_samples=5,
    metric="minkowski",
    leaf_size=40,
    approx_min_span_tree=True,
    n_jobs=4,
    **kwargs,
):
    if leaf_size < 3:
        leaf_size = 3

    if n_jobs < 1:
        n_jobs = max(cpu_count() + 1 + n_jobs, 1)

    if X.dtype != np.float64:
        X = X.astype(np.float64)

    tree = BallTree(X, metric=metric, leaf_size=leaf_size, **kwargs)
    alg = BallTreeBoruvkaAlgorithm(
        tree,
        min_samples,
        metric=metric,
        leaf_size=leaf_size // 3,
        approx_min_span_tree=approx_min_span_tree,
        n_jobs=n_jobs,
        **kwargs,
    )
    min_spanning_tree = alg.spanning_tree()
    # Sort edges of the min_spanning_tree by weight
    min_spanning_tree = min_spanning_tree[np.argsort(min_spanning_tree.T[2]), :]
    # Convert edge list into standard hierarchical clustering format
    single_linkage_tree = label(min_spanning_tree)

    return single_linkage_tree


def remap_condensed_tree(tree, internal_to_raw, outliers):
    """
    Takes an internal condensed_tree structure and adds back in a set of points
    that were initially detected as non-finite and returns that new tree.
    These points will all be split off from the maximal node at lambda zero and
    considered noise points.

    Parameters
    ----------
    tree: condensed_tree
    internal_to_raw: dict
        a mapping from internal integer index to the raw integer index
    finite_index: ndarray
        Boolean array of which entries in the raw data were finite
    """
    finite_count = len(internal_to_raw)

    outlier_count = len(outliers)
    for i, (parent, child, lambda_val, child_size) in enumerate(tree):
        if child < finite_count:
            child = internal_to_raw[child]
        else:
            child = child + outlier_count
        tree[i] = (parent + outlier_count, child, lambda_val, child_size)

    outlier_list = []
    root = tree[0][0]  # Should I check to be sure this is the minimal lambda?
    for outlier in outliers:
        outlier_list.append((root, outlier, 0, 1))

    outlier_tree = np.array(
        outlier_list,
        dtype=[
            ("parent", np.intp),
            ("child", np.intp),
            ("lambda_val", float),
            ("child_size", np.intp),
        ],
    )
    tree = np.append(outlier_tree, tree)
    return tree


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


def hdbscan(
    X,
    min_cluster_size=5,
    min_samples=None,
    alpha=1.0,
    cluster_selection_epsilon=0.0,
    max_cluster_size=0,
    metric="minkowski",
    leaf_size=40,
    algorithm="best",
    memory=None,
    approx_min_span_tree=True,
    n_jobs=4,
    cluster_selection_method="eom",
    allow_single_cluster=False,
    match_reference_implementation=False,
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
        defaults to the min_cluster_size.

    alpha : float, default=1.0
        A distance scaling parameter as used in robust single linkage.
        See [2]_ for more information.

    cluster_selection_epsilon : float, default=0.0
        A distance threshold. Clusters below this value will be merged.
        See [3]_ for more information.

    max_cluster_size : int, default=0
        A limit to the size of clusters returned by the eom algorithm.
        Has no effect when using leaf clustering (where clusters are
        usually small regardless) and can also be overridden in rare
        cases by a high value for cluster_selection_epsilon.

    metric : str or callable, default='minkowski'
        The metric to use when calculating distance between instances in a
        feature array.

        * If metric is a string or callable, it must be one of
          the options allowed by `metrics.pairwise.pairwise_distances` for its
          metric parameter.

        * If metric is "precomputed", X is assumed to be a distance matrix and
          must be square.

    leaf_size : int, default=40
        Leaf size for trees responsible for fast nearest
        neighbour queries.

    algorithm : str, default='best'
        Exactly which algorithm to use; hdbscan has variants specialised
        for different characteristics of the data. By default this is set
        to `best` which chooses the "best" algorithm given the nature of
        the data. You can force other options if you believe you know
        better. Options are:
        - `best`
        - `generic`
        - `prims_kdtree`
        - `prims_balltree`
        - `boruvka_kdtree`
        - `boruvka_balltree`

    memory : str, default=None
        Used to cache the output of the computation of the tree.
        By default, no caching is done. If a string is given, it is the
        path to the caching directory.

    approx_min_span_tree : bool, default=True
        Whether to accept an only approximate minimum spanning tree.
        For some algorithms this can provide a significant speedup, but
        the resulting clustering may be of marginally lower quality.
        If you are willing to sacrifice speed for correctness you may want
        to explore this; in general this should be left at the default True.

    n_jobs : int, default=4
        Number of parallel jobs to run in core distance computations (if
        supported by the specific algorithm). For `n_jobs`
        below -1, (n_cpus + 1 + n_jobs) are used.

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
        to t=True will override this and allow single cluster results in
        the case that you feel this is a valid result for your dataset.

    match_reference_implementation : bool, default=False
        There exist some interpretational differences between this
        HDBSCAN* implementation and the original authors reference
        implementation in Java. This can result in very minor differences
        in clustering results. Setting this flag to True will, at a some
        performance cost, ensure that the clustering results match the
        reference implementation.

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

    if type(min_samples) is not int or type(min_cluster_size) is not int:
        raise ValueError("Min samples and min cluster size must be integers!")

    if min_samples <= 0 or min_cluster_size <= 0:
        raise ValueError("Min samples and Min cluster size must be positive integers")

    if min_cluster_size == 1:
        raise ValueError("Min cluster size must be greater than one")

    if type(cluster_selection_epsilon) is int:
        cluster_selection_epsilon = float(cluster_selection_epsilon)

    if type(cluster_selection_epsilon) is not float or cluster_selection_epsilon < 0.0:
        raise ValueError("Epsilon must be a float value greater than or equal to 0!")

    if not isinstance(alpha, float) or alpha <= 0.0:
        raise ValueError("Alpha must be a positive float value greater than 0!")

    if leaf_size < 1:
        raise ValueError("Leaf size must be greater than 0!")

    if match_reference_implementation:
        min_samples = min_samples - 1
        min_cluster_size = min_cluster_size + 1
        approx_min_span_tree = False

    if cluster_selection_method not in ("eom", "leaf"):
        raise ValueError(
            'Invalid Cluster Selection Method: %s\nShould be one of: "eom", "leaf"\n'
        )

    # Checks input and converts to an nd-array where possible
    if metric != "precomputed" or issparse(X):
        X = check_array(X, accept_sparse="csr", force_all_finite=False)
    else:
        # Only non-sparse, precomputed distance matrices are handled here
        # and thereby allowed to contain numpy.inf for missing distances

        # Perform check_array(X) after removing infinite values (numpy.inf)
        # from the given distance matrix.
        tmp = X.copy()
        tmp[np.isinf(tmp)] = 1
        check_array(tmp)

    # Python 2 and 3 compliant string_type checking
    memory = Memory(location=memory, verbose=0)

    size = X.shape[0]
    min_samples = min(size - 1, min_samples)
    if min_samples == 0:
        min_samples = 1

    metric_params = metric_params or {}
    if algorithm != "best":
        if metric != "precomputed" and issparse(X) and algorithm != "generic":
            raise ValueError("Sparse data matrices only support algorithm 'generic'.")

        if algorithm == "generic":
            single_linkage_tree = memory.cache(_hdbscan_generic)(
                X,
                min_samples,
                alpha,
                metric,
                leaf_size,
                **metric_params,
            )
        elif algorithm == "prims_kdtree":
            if metric not in KDTree.valid_metrics:
                raise ValueError("Cannot use Prim's with KDTree for this metric!")
            single_linkage_tree = memory.cache(_hdbscan_prims_kdtree)(
                X,
                min_samples,
                alpha,
                metric,
                leaf_size,
                **metric_params,
            )
        elif algorithm == "prims_balltree":
            if metric not in BallTree.valid_metrics:
                raise ValueError("Cannot use Prim's with BallTree for this metric!")
            single_linkage_tree = memory.cache(_hdbscan_prims_balltree)(
                X,
                min_samples,
                alpha,
                metric,
                leaf_size,
                **metric_params,
            )
        elif algorithm == "boruvka_kdtree":
            if metric not in KDTree.valid_metrics:
                raise ValueError("Cannot use Boruvka with KDTree for this metric!")
            single_linkage_tree = memory.cache(_hdbscan_boruvka_kdtree)(
                X,
                min_samples,
                metric,
                leaf_size,
                approx_min_span_tree,
                n_jobs,
                **metric_params,
            )
        elif algorithm == "boruvka_balltree":
            if metric not in BallTree.valid_metrics:
                raise ValueError("Cannot use Boruvka with BallTree for this metric!")
            if (X.shape[0] // leaf_size) > 16000:
                warn(
                    "A large dataset size and small leaf_size may induce excessive "
                    "memory usage. If you are running out of memory consider "
                    "increasing the `leaf_size` parameter."
                )
            single_linkage_tree = memory.cache(_hdbscan_boruvka_balltree)(
                X,
                min_samples,
                metric,
                leaf_size,
                approx_min_span_tree,
                n_jobs,
                **metric_params,
            )
        else:
            raise TypeError("Unknown algorithm type %s specified" % algorithm)
    else:

        if issparse(X) or metric not in FAST_METRICS:
            # We can't do much with sparse matrices ...
            single_linkage_tree = memory.cache(_hdbscan_generic)(
                X,
                min_samples,
                alpha,
                metric,
                leaf_size,
                **metric_params,
            )
        elif metric in KDTree.valid_metrics:
            # TO DO: Need heuristic to decide when to go to boruvka;
            # still debugging for now
            if X.shape[1] > 60:
                single_linkage_tree = memory.cache(_hdbscan_prims_kdtree)(
                    X,
                    min_samples,
                    alpha,
                    metric,
                    leaf_size,
                    **metric_params,
                )
            else:
                single_linkage_tree = memory.cache(_hdbscan_boruvka_kdtree)(
                    X,
                    min_samples,
                    metric,
                    leaf_size,
                    approx_min_span_tree,
                    n_jobs,
                    **metric_params,
                )
        else:  # Metric is a valid BallTree metric
            # TO DO: Need heuristic to decide when to go to boruvka;
            # still debugging for now
            if X.shape[1] > 60:
                single_linkage_tree = memory.cache(_hdbscan_prims_balltree)(
                    X,
                    min_samples,
                    alpha,
                    metric,
                    leaf_size,
                    **metric_params,
                )
            else:
                single_linkage_tree = memory.cache(_hdbscan_boruvka_balltree)(
                    X,
                    min_samples,
                    metric,
                    leaf_size,
                    approx_min_span_tree,
                    n_jobs,
                    **metric_params,
                )

    return _tree_to_labels(
        single_linkage_tree,
        min_cluster_size,
        cluster_selection_method,
        allow_single_cluster,
        match_reference_implementation,
        cluster_selection_epsilon,
        max_cluster_size,
    )


# Inherits from sklearn
class HDBSCAN(BaseEstimator, ClusterMixin):
    """Perform HDBSCAN clustering from vector array or distance matrix.

    HDBSCAN - Hierarchical Density-Based Spatial Clustering of Applications
    with Noise. Performs DBSCAN over varying epsilon values and integrates
    the result to find a clustering that gives the best stability over epsilon.
    This allows HDBSCAN to find clusters of varying densities (unlike DBSCAN),
    and be more robust to parameter selection.

    Parameters
    ----------
    min_cluster_size : int, default=5
        The minimum size of clusters; single linkage splits that contain
        fewer points than this will be considered points "falling out" of a
        cluster rather than a cluster splitting into two new clusters.

    min_samples : int, default=None
        The number of samples in a neighbourhood for a point to be
        considered a core point.

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

    algorithm : str, default='best'
        Exactly which algorithm to use; hdbscan has variants specialised
        for different characteristics of the data. By default this is set
        to `best` which chooses the "best" algorithm given the nature of
        the data. You can force other options if you believe you know
        better. Options are:
        - `best`
        - `generic`
        - `prims_kdtree`
        - `prims_balltree`
        - `boruvka_kdtree`
        - `boruvka_balltree`

    leaf_size : int, default=40
        If using a space tree algorithm (`KDTree`, or `BallTree`) the number
        of points ina leaf node of the tree. This does not alter the
        resulting clustering, but may have an effect on the runtime
        of the algorithm.

    memory : str, default=None
        Used to cache the output of the computation of the tree.
        By default, no caching is done. If a string is given, it is the
        path to the caching directory.

    approx_min_span_tree : bool, default=True
        Whether to accept an only approximate minimum spanning tree.
        For some algorithms this can provide a significant speedup, but
        the resulting clustering may be of marginally lower quality.
        If you are willing to sacrifice speed for correctness you may want
        to explore this; in general this should be left at the default `True`.

    n_jobs : int, default=4
        Number of parallel jobs to run in core distance computations (if
        supported by the specific algorithm). For `n_jobs`
        below -1, (n_cpus + 1 + n_jobs) are used.

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

    match_reference_implementation : bool, default=False
        There exist some interpretational differences between this
        HDBSCAN* implementation and the original authors reference
        implementation in Java. This can result in very minor differences
        in clustering results. Setting this flag to True will, at a some
        performance cost, ensure that the clustering results match the
        reference implementation.

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

    def _more_tags(self):
        return {"allow_nan": True}

    def __init__(
        self,
        min_cluster_size=5,
        min_samples=None,
        cluster_selection_epsilon=0.0,
        max_cluster_size=0,
        metric="euclidean",
        alpha=1.0,
        algorithm="best",
        leaf_size=40,
        memory=None,
        approx_min_span_tree=True,
        n_jobs=4,
        cluster_selection_method="eom",
        allow_single_cluster=False,
        match_reference_implementation=False,
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
        self.approx_min_span_tree = approx_min_span_tree
        self.n_jobs = n_jobs
        self.cluster_selection_method = cluster_selection_method
        self.allow_single_cluster = allow_single_cluster
        self.match_reference_implementation = match_reference_implementation
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
        metric_params = self.metric_params or {}
        if self.metric != "precomputed":
            # Non-precomputed matrices may contain non-finite values.
            # Rows with these values
            X = self._validate_data(X, force_all_finite=False, accept_sparse="csr")
            self._raw_data = X

            self._all_finite = (
                np.alltrue(np.isfinite(X.tocoo().data))
                if issparse(X)
                else np.alltrue(np.isfinite(X))
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
