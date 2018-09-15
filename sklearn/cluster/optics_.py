# -*- coding: utf-8 -*-
"""Ordering Points To Identify the Clustering Structure (OPTICS)

These routines execute the OPTICS algorithm, and implement various
cluster extraction methods of the ordered list.

Authors: Shane Grigsby <refuge@rocktalus.com>
         Amy X. Zhang <axz@mit.edu>
         Adrin Jalali <adrin.jalali@gmail.com>
License: BSD 3 clause
"""

from __future__ import division
import warnings
import numpy as np

from ..utils import check_array
from ..utils.validation import check_is_fitted
from ..neighbors import NearestNeighbors
from ..base import BaseEstimator, ClusterMixin
from ..metrics import pairwise_distances
from ._optics_inner import quick_scan


def optics(X, min_samples=5, max_eps=np.inf, metric='euclidean',
           p=2, metric_params=None, maxima_ratio=.75,
           rejection_ratio=.7, similarity_threshold=0.4,
           significant_min=.003, min_cluster_size=.005,
           min_maxima_ratio=0.001, algorithm='ball_tree',
           leaf_size=30, n_jobs=None):
    """Perform OPTICS clustering from vector array

    OPTICS: Ordering Points To Identify the Clustering Structure
    Equivalent to DBSCAN, finds core sample of high density and expands
    clusters from them. Unlike DBSCAN, keeps cluster hierarchy for a variable
    neighborhood radius. Optimized for usage on large point datasets.

    Read more in the :ref:`User Guide <optics>`.

    Parameters
    ----------
    X : array, shape (n_samples, n_features)
        The data.

    min_samples : int (default=5)
        The number of samples in a neighborhood for a point to be considered
        as a core point.

    max_eps : float, optional (default=np.inf)
        The maximum distance between two samples for them to be considered
        as in the same neighborhood. Default value of "np.inf" will identify
        clusters across all scales; reducing `max_eps` will result in
        shorter run times.

    metric : string or callable, optional (default='euclidean')
        The distance metric to use for neighborhood lookups. Default is
        "euclidean". Other options include "minkowski", "manhattan",
        "chebyshev", "haversine", "seuclidean", "hamming", "canberra",
        and "braycurtis". The "wminkowski" and "mahalanobis" metrics are
        also valid with an additional argument.

    p : integer, optional (default=2)
        Parameter for the Minkowski metric from
        :class:`sklearn.metrics.pairwise_distances`. When p = 1, this is
        equivalent to using manhattan_distance (l1), and euclidean_distance
        (l2) for p = 2. For arbitrary p, minkowski_distance (l_p) is used.

    metric_params : dict, optional (default=None)
        Additional keyword arguments for the metric function.

    maxima_ratio : float, optional (default=.75)
        The maximum ratio we allow of average height of clusters on the
        right and left to the local maxima in question. The higher the
        ratio, the more generous the algorithm is to preserving local
        minima, and the more cuts the resulting tree will have.

    rejection_ratio : float, optional (default=.7)
        Adjusts the fitness of the clustering. When the maxima_ratio is
        exceeded, determine which of the clusters to the left and right to
        reject based on rejection_ratio. Higher values will result in points
        being more readily classified as noise; conversely, lower values will
        result in more points being clustered.

    similarity_threshold : float, optional (default=.4)
        Used to check if nodes can be moved up one level, that is, if the
        new cluster created is too "similar" to its parent, given the
        similarity threshold. Similarity can be determined by 1) the size
        of the new cluster relative to the size of the parent node or
        2) the average of the reachability values of the new cluster
        relative to the average of the reachability values of the parent
        node. A lower value for the similarity threshold means less levels
        in the tree.

    significant_min : float, optional (default=.003)
        Sets a lower threshold on how small a significant maxima can be.

    min_cluster_size : int > 1 or float between 0 and 1 (default=0.005)
        Minimum number of samples in an OPTICS cluster, expressed as an
        absolute number or a fraction of the number of samples (rounded
        to be at least 2).

    min_maxima_ratio : float, optional (default=.001)
        Used to determine neighborhood size for minimum cluster membership.
        Each local maxima should be a largest value in a neighborhood
        of the `size min_maxima_ratio * len(X)` from left and right.

    algorithm : {'auto', 'ball_tree', 'kd_tree', 'brute'}, optional
        Algorithm used to compute the nearest neighbors:

        - 'ball_tree' will use :class:`BallTree` (default)
        - 'kd_tree' will use :class:`KDTree`
        - 'brute' will use a brute-force search.
        - 'auto' will attempt to decide the most appropriate algorithm
          based on the values passed to :meth:`fit` method.

        Note: fitting on sparse input will override the setting of
        this parameter, using brute force.

    leaf_size : int, optional (default=30)
        Leaf size passed to :class:`BallTree` or :class:`KDTree`. This can
        affect the speed of the construction and query, as well as the memory
        required to store the tree. The optimal value depends on the
        nature of the problem.

    n_jobs : int or None, optional (default=None)
        The number of parallel jobs to run for neighbors search.
        ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
        ``-1`` means using all processors. See :term:`Glossary <n_jobs>`
        for more details.

    Returns
    -------
    core_sample_indices_ : array, shape (n_core_samples,)
        The indices of the core samples.

    labels_ : array, shape (n_samples,)
        The estimated labels.

    See also
    --------
    OPTICS
        An estimator interface for this clustering algorithm.
    dbscan
        A similar clustering for a specified neighborhood radius (eps).
        Our implementation is optimized for runtime.

    References
    ----------
    Ankerst, Mihael, Markus M. Breunig, Hans-Peter Kriegel, and Jörg Sander.
    "OPTICS: ordering points to identify the clustering structure." ACM SIGMOD
    Record 28, no. 2 (1999): 49-60.
    """

    clust = OPTICS(min_samples, max_eps, metric, p, metric_params,
                   maxima_ratio, rejection_ratio,
                   similarity_threshold, significant_min,
                   min_cluster_size, min_maxima_ratio,
                   algorithm, leaf_size, n_jobs)
    clust.fit(X)
    return clust.core_sample_indices_, clust.labels_


class OPTICS(BaseEstimator, ClusterMixin):
    """Estimate clustering structure from vector array

    OPTICS: Ordering Points To Identify the Clustering Structure
    Equivalent to DBSCAN, finds core sample of high density and expands
    clusters from them. Unlike DBSCAN, keeps cluster hierarchy for a variable
    neighborhood radius. Optimized for usage on large point datasets.

    Read more in the :ref:`User Guide <optics>`.

    Parameters
    ----------
    min_samples : int (default=5)
        The number of samples in a neighborhood for a point to be considered
        as a core point.

    max_eps : float, optional (default=np.inf)
        The maximum distance between two samples for them to be considered
        as in the same neighborhood. Default value of "np.inf" will identify
        clusters across all scales; reducing `max_eps` will result in
        shorter run times.

    metric : string or callable, optional (default='euclidean')
        The distance metric to use for neighborhood lookups. Default is
        "euclidean". Other options include "minkowski", "manhattan",
        "chebyshev", "haversine", "seuclidean", "hamming", "canberra",
        and "braycurtis". The "wminkowski" and "mahalanobis" metrics are
        also valid with an additional argument.

    p : integer, optional (default=2)
        Parameter for the Minkowski metric from
        :class:`sklearn.metrics.pairwise_distances`. When p = 1, this is
        equivalent to using manhattan_distance (l1), and euclidean_distance
        (l2) for p = 2. For arbitrary p, minkowski_distance (l_p) is used.

    metric_params : dict, optional (default=None)
        Additional keyword arguments for the metric function.

    maxima_ratio : float, optional (default=.75)
        The maximum ratio we allow of average height of clusters on the
        right and left to the local maxima in question. The higher the
        ratio, the more generous the algorithm is to preserving local
        minima, and the more cuts the resulting tree will have.

    rejection_ratio : float, optional (default=.7)
        Adjusts the fitness of the clustering. When the maxima_ratio is
        exceeded, determine which of the clusters to the left and right to
        reject based on rejection_ratio. Higher values will result in points
        being more readily classified as noise; conversely, lower values will
        result in more points being clustered.

    similarity_threshold : float, optional (default=.4)
        Used to check if nodes can be moved up one level, that is, if the
        new cluster created is too "similar" to its parent, given the
        similarity threshold. Similarity can be determined by 1) the size
        of the new cluster relative to the size of the parent node or
        2) the average of the reachability values of the new cluster
        relative to the average of the reachability values of the parent
        node. A lower value for the similarity threshold means less levels
        in the tree.

    significant_min : float, optional (default=.003)
        Sets a lower threshold on how small a significant maxima can be.

    min_cluster_size : int > 1 or float between 0 and 1 (default=0.005)
        Minimum number of samples in an OPTICS cluster, expressed as an
        absolute number or a fraction of the number of samples (rounded
        to be at least 2).

    min_maxima_ratio : float, optional (default=.001)
        Used to determine neighborhood size for minimum cluster membership.
        Each local maxima should be a largest value in a neighborhood
        of the `size min_maxima_ratio * len(X)` from left and right.

    algorithm : {'auto', 'ball_tree', 'kd_tree', 'brute'}, optional
        Algorithm used to compute the nearest neighbors:

        - 'ball_tree' will use :class:`BallTree` (default)
        - 'kd_tree' will use :class:`KDTree`
        - 'brute' will use a brute-force search.
        - 'auto' will attempt to decide the most appropriate algorithm
          based on the values passed to :meth:`fit` method.

        Note: fitting on sparse input will override the setting of
        this parameter, using brute force.

    leaf_size : int, optional (default=30)
        Leaf size passed to :class:`BallTree` or :class:`KDTree`. This can
        affect the speed of the construction and query, as well as the memory
        required to store the tree. The optimal value depends on the
        nature of the problem.

    n_jobs : int or None, optional (default=None)
        The number of parallel jobs to run for neighbors search.
        ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
        ``-1`` means using all processors. See :term:`Glossary <n_jobs>`
        for more details.

    Attributes
    ----------
    core_sample_indices_ : array, shape (n_core_samples,)
        Indices of core samples.

    labels_ : array, shape (n_samples,)
        Cluster labels for each point in the dataset given to fit().
        Noisy samples are given the label -1.

    reachability_ : array, shape (n_samples,)
        Reachability distances per sample.

    ordering_ : array, shape (n_samples,)
        The cluster ordered list of sample indices

    core_distances_ : array, shape (n_samples,)
        Distance at which each sample becomes a core point.
        Points which will never be core have a distance of inf.

    See also
    --------

    DBSCAN
        A similar clustering for a specified neighborhood radius (eps).
        Our implementation is optimized for runtime.

    References
    ----------
    Ankerst, Mihael, Markus M. Breunig, Hans-Peter Kriegel, and Jörg Sander.
    "OPTICS: ordering points to identify the clustering structure." ACM SIGMOD
    Record 28, no. 2 (1999): 49-60.
    """

    def __init__(self, min_samples=5, max_eps=np.inf, metric='euclidean',
                 p=2, metric_params=None, maxima_ratio=.75,
                 rejection_ratio=.7, similarity_threshold=0.4,
                 significant_min=.003, min_cluster_size=.005,
                 min_maxima_ratio=0.001, algorithm='ball_tree',
                 leaf_size=30, n_jobs=None):

        self.max_eps = max_eps
        self.min_samples = min_samples
        self.maxima_ratio = maxima_ratio
        self.rejection_ratio = rejection_ratio
        self.similarity_threshold = similarity_threshold
        self.significant_min = significant_min
        self.min_cluster_size = min_cluster_size
        self.min_maxima_ratio = min_maxima_ratio
        self.algorithm = algorithm
        self.metric = metric
        self.metric_params = metric_params
        self.p = p
        self.leaf_size = leaf_size
        self.n_jobs = n_jobs

    def fit(self, X, y=None):
        """Perform OPTICS clustering

        Extracts an ordered list of points and reachability distances, and
        performs initial clustering using `max_eps` distance specified at
        OPTICS object instantiation.

        Parameters
        ----------
        X : array, shape (n_samples, n_features)
            The data.

        y : ignored

        Returns
        -------
        self : instance of OPTICS
            The instance.
        """
        X = check_array(X, dtype=np.float)

        n_samples = len(X)

        if self.min_samples > n_samples:
            raise ValueError("Number of training samples (n_samples=%d) must "
                             "be greater than min_samples (min_samples=%d) "
                             "used for clustering." %
                             (n_samples, self.min_samples))

        if self.min_cluster_size <= 0 or (self.min_cluster_size !=
                                          int(self.min_cluster_size)
                                          and self.min_cluster_size > 1):
            raise ValueError('min_cluster_size must be a positive integer or '
                             'a float between 0 and 1. Got %r' %
                             self.min_cluster_size)
        elif self.min_cluster_size > n_samples:
            raise ValueError('min_cluster_size must be no greater than the '
                             'number of samples (%d). Got %d' %
                             (n_samples, self.min_cluster_size))

        # Start all points as 'unprocessed' ##
        self.reachability_ = np.empty(n_samples)
        self.reachability_.fill(np.inf)
        self.core_distances_ = np.empty(n_samples)
        self.core_distances_.fill(np.nan)
        # Start all points as noise ##
        self.labels_ = np.full(n_samples, -1, dtype=int)

        nbrs = NearestNeighbors(n_neighbors=self.min_samples,
                                algorithm=self.algorithm,
                                leaf_size=self.leaf_size, metric=self.metric,
                                metric_params=self.metric_params, p=self.p,
                                n_jobs=self.n_jobs)

        nbrs.fit(X)
        self.core_distances_[:] = nbrs.kneighbors(X,
                                                  self.min_samples)[0][:, -1]

        self.ordering_ = self._calculate_optics_order(X, nbrs)

        indices_, self.labels_ = _extract_optics(self.ordering_,
                                                 self.reachability_,
                                                 self.maxima_ratio,
                                                 self.rejection_ratio,
                                                 self.similarity_threshold,
                                                 self.significant_min,
                                                 self.min_cluster_size,
                                                 self.min_maxima_ratio)
        self.core_sample_indices_ = indices_
        return self

    # OPTICS helper functions

    def _calculate_optics_order(self, X, nbrs):
        # Main OPTICS loop. Not parallelizable. The order that entries are
        # written to the 'ordering_' list is important!
        processed = np.zeros(X.shape[0], dtype=bool)
        ordering = np.zeros(X.shape[0], dtype=int)
        ordering_idx = 0
        for point in range(X.shape[0]):
            if processed[point]:
                continue
            if self.core_distances_[point] <= self.max_eps:
                while not processed[point]:
                    processed[point] = True
                    ordering[ordering_idx] = point
                    ordering_idx += 1
                    point = self._set_reach_dist(point, processed, X, nbrs)
            else:  # For very noisy points
                ordering[ordering_idx] = point
                ordering_idx += 1
                processed[point] = True
        return ordering

    def _set_reach_dist(self, point_index, processed, X, nbrs):
        P = X[point_index:point_index + 1]
        indices = nbrs.radius_neighbors(P, radius=self.max_eps,
                                        return_distance=False)[0]

        # Getting indices of neighbors that have not been processed
        unproc = np.compress((~np.take(processed, indices)).ravel(),
                             indices, axis=0)
        # Keep n_jobs = 1 in the following lines...please
        if not unproc.size:
            # Everything is already processed. Return to main loop
            return point_index

        dists = pairwise_distances(P, np.take(X, unproc, axis=0),
                                   self.metric, n_jobs=1).ravel()

        rdists = np.maximum(dists, self.core_distances_[point_index])
        new_reach = np.minimum(np.take(self.reachability_, unproc), rdists)
        self.reachability_[unproc] = new_reach

        # Define return order based on reachability distance
        return (unproc[quick_scan(np.take(self.reachability_, unproc),
                                  dists)])

    def extract_dbscan(self, eps):
        """Performs DBSCAN extraction for an arbitrary epsilon.

        Extraction runs in linear time. Note that if the `max_eps` OPTICS
        parameter was set to < inf for extracting reachability and ordering
        arrays, DBSCAN extractions will be unstable for `eps` values close to
        `max_eps`. Setting `eps` < (`max_eps` / 5.0) will guarantee
        extraction parity with DBSCAN.

        Parameters
        ----------
        eps : float or int, required
            DBSCAN `eps` parameter. Must be set to < `max_eps`. Equivalence
            with DBSCAN algorithm is achieved if `eps` is < (`max_eps` / 5)

        Returns
        -------
        core_sample_indices_ : array, shape (n_core_samples,)
            The indices of the core samples.

        labels_ : array, shape (n_samples,)
            The estimated labels.
        """
        check_is_fitted(self, 'reachability_')

        if eps > self.max_eps:
            raise ValueError('Specify an epsilon smaller than %s. Got %s.'
                             % (self.max_eps, eps))

        if eps * 5.0 > (self.max_eps * 1.05):
            warnings.warn(
                "Warning, max_eps (%s) is close to eps (%s): "
                "Output may be unstable." % (self.max_eps, eps),
                RuntimeWarning, stacklevel=2)
        # Stability warning is documented in _extract_dbscan method...

        return _extract_dbscan(self.ordering_, self.core_distances_,
                               self.reachability_, eps)

    def extract_xi(self, xi, return_clusters=False):
        """Automatically extract clusters according to the
        :func:`Xi-steep method _extract_xi`.

        Parameters
        ----------
        xi : float, between 0 and 1
            The main parameter in the Xi-steep method.

        return_clusters : bool (default=False)
            Return the clusters as well as the labels. If `False`, it will
            only return the labels.

        Returns
        -------
        labels : array, shape (n_samples)
            The labels assigned to samples. Points which are not included
            in any cluster are labeled as -1.

        clusters : list, if return_clusters is True
            The list of clusters in the form of (start, end) tuples,
            with all indices inclusive. The clusters are ordered in a way that
            larger clusters encompassing smaller clusters, come after those
            smaller clusters.
        """
        check_is_fitted(self, 'reachability_')

        return _extract_xi(self.reachability_,
                           self.ordering_,
                           self.min_samples,
                           self.min_cluster_size,
                           xi, return_clusters)


def _extract_dbscan(ordering, core_distances, reachability, eps):
    """Performs DBSCAN extraction for an arbitrary epsilon (`eps`).

    Parameters
    ----------
    ordering : array, shape (n_samples,)
        OPTICS ordered point indices (`ordering_`)
    core_distances : array, shape (n_samples,)
        Distances at which points become core (`core_distances_`)
    reachability : array, shape (n_samples,)
        Reachability distances calculated by OPTICS (`reachability_`)
    eps : float or int
        DBSCAN `eps` parameter

    Returns
    -------
    core_sample_indices_ : array, shape (n_core_samples,)
        The indices of the core samples.

    labels_ : array, shape (n_samples,)
        The estimated labels.
    """

    n_samples = len(core_distances)
    is_core = np.zeros(n_samples, dtype=bool)
    labels = np.zeros(n_samples, dtype=int)

    far_reach = reachability > eps
    near_core = core_distances <= eps
    labels[ordering] = np.cumsum(far_reach[ordering] & near_core[ordering]) - 1
    labels[far_reach & ~near_core] = -1
    is_core[near_core] = True
    return np.arange(n_samples)[is_core], labels


def _extract_optics(ordering, reachability, maxima_ratio=.75,
                    rejection_ratio=.7, similarity_threshold=0.4,
                    significant_min=.003, min_cluster_size=.005,
                    min_maxima_ratio=0.001):
    """Performs automatic cluster extraction for variable density data.

    Parameters
    ----------
    ordering : array, shape (n_samples,)
        OPTICS ordered point indices (`ordering_`)

    reachability : array, shape (n_samples,)
        Reachability distances calculated by OPTICS (`reachability_`)

    maxima_ratio : float, optional
        The maximum ratio we allow of average height of clusters on the
        right and left to the local maxima in question. The higher the
        ratio, the more generous the algorithm is to preserving local
        minima, and the more cuts the resulting tree will have.

    rejection_ratio : float, optional
        Adjusts the fitness of the clustering. When the maxima_ratio is
        exceeded, determine which of the clusters to the left and right to
        reject based on rejection_ratio. Higher values will result in points
        being more readily classified as noise; conversely, lower values will
        result in more points being clustered.

    similarity_threshold : float, optional
        Used to check if nodes can be moved up one level, that is, if the
        new cluster created is too "similar" to its parent, given the
        similarity threshold. Similarity can be determined by 1) the size
        of the new cluster relative to the size of the parent node or
        2) the average of the reachability values of the new cluster
        relative to the average of the reachability values of the parent
        node. A lower value for the similarity threshold means less levels
        in the tree.

    significant_min : float, optional
        Sets a lower threshold on how small a significant maxima can be.

    min_cluster_size : int > 1 or float between 0 and 1
        Minimum number of samples in an OPTICS cluster, expressed as an
        absolute number or a fraction of the number of samples (rounded
        to be at least 2).

    min_maxima_ratio : float, optional
        Used to determine neighborhood size for minimum cluster membership.

    Returns
    -------
    core_sample_indices_ : array, shape (n_core_samples,)
        The indices of the core samples.

    labels_ : array, shape (n_samples,)
        The estimated labels.
    """

    # Extraction wrapper
    reachability = reachability / np.max(reachability[1:])
    reachability_plot = reachability[ordering].tolist()
    root_node = _automatic_cluster(reachability_plot, ordering,
                                   maxima_ratio, rejection_ratio,
                                   similarity_threshold, significant_min,
                                   min_cluster_size, min_maxima_ratio)
    leaves = _get_leaves(root_node, [])
    # Start cluster id's at 0
    clustid = 0
    n_samples = len(reachability)
    is_core = np.zeros(n_samples, dtype=bool)
    labels = np.full(n_samples, -1, dtype=int)
    # Start all points as non-core noise
    for leaf in leaves:
        index = ordering[leaf.start:leaf.end]
        labels[index] = clustid
        is_core[index] = 1
        clustid += 1
    return np.arange(n_samples)[is_core], labels


def _automatic_cluster(reachability_plot, ordering,
                       maxima_ratio, rejection_ratio,
                       similarity_threshold, significant_min,
                       min_cluster_size, min_maxima_ratio):
    """Converts reachability plot to cluster tree and returns root node.

    Parameters
    ----------

    reachability_plot : list, required
        Reachability distances ordered by OPTICS ordering index.

    """

    min_neighborhood_size = 2
    if min_cluster_size <= 1:
        min_cluster_size = max(2, min_cluster_size * len(ordering))
    neighborhood_size = int(min_maxima_ratio * len(ordering))

    # Again, should this check < min_samples, should the parameter be public?
    if neighborhood_size < min_neighborhood_size:
        neighborhood_size = min_neighborhood_size

    local_maxima_points = _find_local_maxima(reachability_plot,
                                             neighborhood_size)
    root_node = _TreeNode(ordering, 0, len(ordering), None)
    _cluster_tree(root_node, None, local_maxima_points,
                  reachability_plot, ordering, min_cluster_size,
                  maxima_ratio, rejection_ratio,
                  similarity_threshold, significant_min)

    return root_node


class _TreeNode(object):
    # automatic cluster helper classes and functions
    def __init__(self, points, start, end, parent_node):
        self.points = points
        self.start = start
        self.end = end
        self.parent_node = parent_node
        self.children = []
        self.split_point = -1


def _is_local_maxima(index, reachability_plot, neighborhood_size):
    right_idx = slice(index + 1, index + neighborhood_size + 1)
    left_idx = slice(max(1, index - neighborhood_size - 1), index)
    return (np.all(reachability_plot[index] >= reachability_plot[left_idx]) and
            np.all(reachability_plot[index] >= reachability_plot[right_idx]))


def _find_local_maxima(reachability_plot, neighborhood_size):
    local_maxima_points = {}
    # 1st and last points on Reachability Plot are not taken
    # as local maxima points
    for i in range(1, len(reachability_plot) - 1):
        # if the point is a local maxima on the reachability plot with
        # regard to neighborhood_size, insert it into priority queue and
        # maxima list
        if (reachability_plot[i] > reachability_plot[i - 1] and
            reachability_plot[i] >= reachability_plot[i + 1] and
            _is_local_maxima(i, np.array(reachability_plot),
                             neighborhood_size) == 1):
            local_maxima_points[i] = reachability_plot[i]

    return sorted(local_maxima_points,
                  key=local_maxima_points.__getitem__, reverse=True)


def _cluster_tree(node, parent_node, local_maxima_points,
                  reachability_plot, reachability_ordering,
                  min_cluster_size, maxima_ratio, rejection_ratio,
                  similarity_threshold, significant_min):
    """Recursively builds cluster tree to hold hierarchical cluster structure

    node is a node or the root of the tree in the first call
    parent_node is parent node of N or None if node is root of the tree
    local_maxima_points is list of local maxima points sorted in
    descending order of reachability
    """

    if len(local_maxima_points) == 0:
        return  # parent_node is a leaf

    # take largest local maximum as possible separation between clusters
    s = local_maxima_points[0]
    node.split_point = s
    local_maxima_points = local_maxima_points[1:]

    # create two new nodes and add to list of nodes
    node_1 = _TreeNode(reachability_ordering[node.start:s],
                       node.start, s, node)
    node_2 = _TreeNode(reachability_ordering[s + 1:node.end],
                       s + 1, node.end, node)
    local_max_1 = []
    local_max_2 = []

    for i in local_maxima_points:
        if i < s:
            local_max_1.append(i)
        if i > s:
            local_max_2.append(i)

    node_list = []
    node_list.append((node_1, local_max_1))
    node_list.append((node_2, local_max_2))

    if reachability_plot[s] < significant_min:
        node.split_point = -1
        # if split_point is not significant, ignore this split and continue
        return

    # only check a certain ratio of points in the child
    # nodes formed to the left and right of the maxima
    # ...should check_ratio be a user settable parameter?
    check_ratio = .8
    check_value_1 = int(np.round(check_ratio * len(node_1.points)))
    check_value_2 = int(np.round(check_ratio * len(node_2.points)))
    avg_reach1 = np.mean(reachability_plot[(node_1.end -
                                            check_value_1):node_1.end])
    avg_reach2 = np.mean(reachability_plot[node_2.start:(node_2.start
                                                         + check_value_2)])

    if ((avg_reach1 / reachability_plot[s]) > maxima_ratio or
            (avg_reach2 / reachability_plot[s]) > maxima_ratio):

        if (avg_reach1 / reachability_plot[s]) < rejection_ratio:
            # reject node 2
            node_list.remove((node_2, local_max_2))
        if (avg_reach2 / reachability_plot[s]) < rejection_ratio:
            # reject node 1
            node_list.remove((node_1, local_max_1))
        if ((avg_reach1 / reachability_plot[s]) >= rejection_ratio and
                (avg_reach2 / reachability_plot[s]) >= rejection_ratio):
            # since split_point is not significant,
            # ignore this split and continue (reject both child nodes)
            node.split_point = -1
            _cluster_tree(node, parent_node, local_maxima_points,
                          reachability_plot, reachability_ordering,
                          min_cluster_size, maxima_ratio, rejection_ratio,
                          similarity_threshold, significant_min)
            return

    # remove clusters that are too small
    if (len(node_1.points) < min_cluster_size and
            node_list.count((node_1, local_max_1)) > 0):
        # cluster 1 is too small
        node_list.remove((node_1, local_max_1))
    if (len(node_2.points) < min_cluster_size and
            node_list.count((node_2, local_max_2)) > 0):
        # cluster 2 is too small
        node_list.remove((node_2, local_max_2))
    if not node_list:
        # parent_node will be a leaf
        node.split_point = -1
        return

    # Check if nodes can be moved up one level - the new cluster created
    # is too "similar" to its parent, given the similarity threshold.
    bypass_node = 0
    if parent_node is not None:
        if ((node.end - node.start) / (parent_node.end - parent_node.start) >
                similarity_threshold):

            parent_node.children.remove(node)
            bypass_node = 1

    for nl in node_list:
        if bypass_node == 1:
            parent_node.children.append(nl[0])
            _cluster_tree(nl[0], parent_node, nl[1],
                          reachability_plot, reachability_ordering,
                          min_cluster_size, maxima_ratio, rejection_ratio,
                          similarity_threshold, significant_min)
        else:
            node.children.append(nl[0])
            _cluster_tree(nl[0], node, nl[1], reachability_plot,
                          reachability_ordering, min_cluster_size,
                          maxima_ratio, rejection_ratio,
                          similarity_threshold, significant_min)


def _get_leaves(node, arr):
    if node is not None:
        if node.split_point == -1:
            arr.append(node)
        for n in node.children:
            _get_leaves(n, arr)
    return arr


def _extract_xi(reachability, ordering, min_samples, min_cluster_size, xi,
                return_clusters):
    """Automatically extract clusters according to the Xi-steep method.

    Parameters
    ----------
    reachability_plot : array, shape (n_samples)
        The reachability plot, i.e. reachability ordered according to
        the calculated ordering, all computed by OPTICS.

    xi : float, between 0 and 1
        The main parameter in the Xi-steep method.

    min_samples : integer
       The same as the min_samples given to OPTICS. Up and down steep
       regions can't have more then `min_samples` consecutive non-steep
       points.

    min_cluster_size : int > 1 or float between 0 and 1
        Minimum number of samples in an OPTICS cluster, expressed as an
        absolute number or a fraction of the number of samples (rounded
        to be at least 2).

    return_clusters : bool
        Return the clusters as well as the labels. If `False`, it will
        only return the labels.

    Returns
    -------
    labels : array, shape (n_samples)
        The labels assigned to samples. Points which are not included
        in any cluster are labeled as -1.

    clusters : list, if return_clusters is True
        The list of clusters in the form of (start, end) tuples,
        with all indices inclusive. The clusters are ordered in a way that
        larger clusters encompassing smaller clusters, come after those
        smaller clusters.
    """
    clusters = _xi_cluster(reachability[ordering], xi, min_samples,
                           min_cluster_size)
    labels = _extract_xi_labels(ordering, clusters)
    if return_clusters:
        return labels, clusters
    else:
        return labels


def _steep_upward(reachability_plot, p, ixi):
    """Check if point p is a xi steep up area (definition 9).
    ixi is the inverse xi, i.e. `1 - xi`"""
    return reachability_plot[p] <= reachability_plot[p + 1] * ixi


def _steep_downward(reachability_plot, p, ixi):
    """Check if point p is a xi steep down area (definition 9).
    ixi is the inverse xi, i.e. `1 - xi`"""
    return reachability_plot[p] * ixi >= reachability_plot[p + 1]


class _Area:
    """An (upward or downward) area.

    Attributes
    ----------
    start : integer
        The start of the region.

    end : integer
        The end of the region.

    maximum : float
        The maximum reachability in this region, which is the
        start of the region for a downward area and the end of
        the region for an upward area.

    mib : float
        Maximum in between value, i.e. the maximum value between
        the end of a steep down area and the current index. It is
        irrelevant for steep up areas.
    """
    def __init__(self, start, end, maximum, mib):
        self.start = start
        self.end = end
        # maximum is the beginning of a downward area
        # and the end of an upward area
        self.maximum = maximum
        # maximum in between (ref Sec. 4.3.2)
        self.mib = mib

    def __repr__(self):
        return ("start: %s, end: %s, max: %4g, mib: %4g" %
                (self.start, self.end, self.maximum, self.mib))


def _extend_downward(reachability_plot, start, ixi, min_samples, n_samples):
    """Extend the downward area until it's maximal.

    Parameters
    ----------
    reachability_plot : array, shape (n_samples)
        The reachability plot, i.e. reachability ordered according to
        the calculated ordering, all computed by OPTICS.

    start : integer
        The start of the downward region.

    ixi: float, between 0 and 1
        The inverse xi, i.e. `1 - xi`

    min_samples : integer
       The same as the min_samples given to OPTICS. Up and down steep
       regions can't have more then `min_samples` consecutive non-steep
       points.

    n_samples : integer
        Total number of samples.

    Returns
    -------
    index : integer
        The current index iterating over all the samples.

    end : integer
        The end of the downward region, which can be behind the index.
    """
    # print("extend down start")
    non_downward_points = 0
    index = start
    end = start
    # find a maximal downward area
    while index + 1 < n_samples:
        # print("index", index)
        # print("r", reachability_plot[index], "r + 1",
        #       reachability_plot[index + 1])
        index += 1
        if _steep_downward(reachability_plot, index, ixi):
            # print("steep")
            non_downward_points = 0
            end = index + 1
        elif reachability_plot[index] >= reachability_plot[index + 1]:
            # print("just down")
            # it's not a steep downward point, but still goes down.
            non_downward_points += 1
            # region should include no more than min_samples consecutive
            # non downward points.
            if non_downward_points == min_samples:
                # print("non downward")
                break
        else:
            break
    # print("extend end")
    return index, end


def _extend_upward(reachability_plot, start, ixi, min_samples, n_samples):
    """Extend the upward area until it's maximal.

    Parameters
    ----------
    reachability_plot : array, shape (n_samples)
        The reachability plot, i.e. reachability ordered according to
        the calculated ordering, all computed by OPTICS.

    start : integer
        The start of the upward region.

    ixi: float, between 0 and 1
        The inverse xi, i.e. `1 - xi`

    min_samples : integer
       The same as the min_samples given to OPTICS. Up and down steep
       regions can't have more then `min_samples` consecutive non-steep
       points.

    n_samples : integer
        Total number of samples.

    Returns
    -------
    index : integer
        The current index iterating over all the samples.

    end : integer
        The end of the upward region, which can be behind the index.
    """
    # print("extend up start")
    non_upward_points = 0
    index = start
    end = start
    # find a maximal upward area
    while index + 1 < n_samples:
        # print("index", index)
        # print("r", reachability_plot[index], "r + 1",
        #       reachability_plot[index + 1])
        index += 1
        if _steep_upward(reachability_plot, index, ixi):
            # print("steep")
            non_upward_points = 0
            end = index + 1
        elif reachability_plot[index] <= reachability_plot[index + 1]:
            # print("just up")
            # it's not a steep upward point, but still goes up.
            non_upward_points += 1
            # region should include no more than min_samples consecutive
            # non downward points.
            if non_upward_points == min_samples:
                # print("non upward")
                break
        else:
            break
    # print("extend end")
    return index, end


def _update_fileter_sdas(sdas, mib, ixi):
    """Update steep down areas (SDAs) using the new
    maximum in between (mib) value, and the given inverse xi, i.e. `1 - xi`
    """
    res = [sda for sda in sdas if mib <= sda.maximum * ixi]
    for sda in res:
        sda.mib = max(sda.mib, mib)
    return res


def _xi_cluster(reachability_plot, xi, min_samples, min_cluster_size):
    """Automatically extract clusters according to the Xi-steep method.

    Parameters
    ----------
    reachability_plot : array, shape (n_samples)
        The reachability plot, i.e. reachability ordered according to
        the calculated ordering, all computed by OPTICS.

    xi : float, between 0 and 1
        The main parameter in the Xi-steep method.

    min_samples : integer
       The same as the min_samples given to OPTICS. Up and down steep
       regions can't have more then `min_samples` consecutive non-steep
       points.

    min_cluster_size : int > 1 or float between 0 and 1
        Minimum number of samples in an OPTICS cluster, expressed as an
        absolute number or a fraction of the number of samples (rounded
        to be at least 2).

    Returns
    -------
    clusters : list
        The list of clusters in the form of (start, end) tuples,
        with all indices inclusive. The clusters are ordered in a way that
        larger clusters encompassing smaller clusters, come after those
        smaller clusters.
    """

    # all indices are inclusive (specially at the end)
    n_samples = len(reachability_plot)
    # add an inf to the end of reachability plot
    # but process the data only until the last actual point
    # this is fine since the last point is considered upward anyway
    reachability_plot = np.array(reachability_plot)
    reachability_plot = np.hstack((reachability_plot, np.inf))

    if min_cluster_size <= 1:
        min_cluster_size = max(2, min_cluster_size * n_samples)

    ixi = 1 - xi
    sdas = list()
    clusters = list()
    index = int(0)
    mib = 0.  # maximum in between
    while index + 1 < n_samples:
        # print("index", index)
        # print("r", reachability_plot[index])
        mib = max(mib, reachability_plot[index])
        # print("mib up there:", mib)

        # check if a steep downward area starts
        if _steep_downward(reachability_plot, index, ixi):
            # print("steep downward")
            # print("sdas", sdas)
            # print("filter mib:", mib)
            sdas = _update_fileter_sdas(sdas, mib, ixi)
            # print("sdas", sdas)
            D_start = index
            index, end = _extend_downward(reachability_plot, D_start, ixi,
                                          min_samples, n_samples)
            D = _Area(start=D_start, end=end,
                      maximum=reachability_plot[D_start], mib=0.)
            # print("D", D, "r.s %.4g" % reachability_plot[D.start],
            #       "r.e %.4g" % reachability_plot[D.end])
            sdas.append(D)
            mib = reachability_plot[index]

        elif _steep_upward(reachability_plot, index, ixi):
            # print("steep upward")
            # print("sdas", sdas)
            # print("filter mib:", mib)
            sdas = _update_fileter_sdas(sdas, mib, ixi)
            # print("sdas", sdas)
            U_start = index
            index, end = _extend_upward(reachability_plot, U_start, ixi,
                                        min_samples, n_samples)
            U = _Area(start=U_start, end=end, maximum=reachability_plot[end],
                      mib=-1)
            # if np.isinf(reachability_plot[index + 1]):
            #     U.maximum = np.inf
            #     index += 1
            # print("U", U, "r.s %.4g" % reachability_plot[U.start],
            #       "r.e %.4g" % reachability_plot[U.end])
            mib = reachability_plot[end - 1]
            # print('mib %.4g' % mib)
            # print(sdas)

            U_clusters = list()
            for D in sdas:
                c_start = D.start
                c_end = min(U.end, n_samples - 1)
                # print("D", D, "U", U)
                # print("start", c_start, "end", c_end)

                # 3.b
                if D.mib > mib * ixi:
                    continue
                # print("3b pass")

                # 4
                if D.maximum * ixi >= U.maximum:
                    while (reachability_plot[c_start + 1] >
                           U.maximum
                           and c_start < c_end):
                        c_start += 1
                elif U.maximum * ixi >= D.maximum:
                    while (reachability_plot[c_end - 1] > D.maximum
                           and c_end > c_start):
                        c_end -= 1
                # print('after 4', c_start, c_end)

                if _steep_upward(reachability_plot, index - 1, ixi):
                    c_end -= 1
                # print('check last point', c_end, 'index', index)

                # 3.a
                if c_end - c_start + 1 < min_cluster_size:
                    continue
                # print('min pts pass')

                U_clusters.append((c_start, c_end))
                # print('U clusters', U_clusters)

            # add smaller clusters first.
            U_clusters.reverse()
            clusters.extend(U_clusters)
            # print("set of clusters:", clusters)

        else:
            # print("just else", index)
            index += 1

    return clusters


def _extract_xi_labels(ordering, clusters):
    """Extracts the labels from the clusters returned by `_xi_cluster`.
    We rely on the fact that clusters are stored
    with the smaller clusters coming before the larger ones.

    Parameters
    ----------
    ordering : array, shape (n_samples)
        The ordering of points calculated by OPTICS

    clusters : list
        List of clusters i.e. (start, end) tuples,
        as returned by `_xi_cluster`.

    Returns
    -------
    labels : array, shape (n_samples)
    """

    labels = np.zeros(len(ordering), dtype=np.int)
    label = 1
    for c in clusters:
        if all(labels[c[0]:(c[1] + 1)] == 0):
            labels[c[0]:(c[1] + 1)] = label
            label += 1
    labels[ordering] = labels - 1
    return labels
