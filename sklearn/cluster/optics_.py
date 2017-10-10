# -*- coding: utf-8 -*-
"""Ordering Points To Identify the Clustering Structure (OPTICS)

These routines execute the OPTICS algorithm, and implement various
cluster extraction methods of the ordered list.

Authors: Shane Grigsby <refuge@rocktalus.com>
         Amy X. Zhang <axz@mit.edu>
License: BSD 3 clause
"""

import warnings

import numpy as np
from scipy import iterable

from ..utils import check_array
from ..utils.validation import check_is_fitted
from ..neighbors.base import NeighborsBase, UnsupervisedMixin
from ..neighbors.base import KNeighborsMixin, RadiusNeighborsMixin
from ..base import ClusterMixin
from ..metrics.pairwise import pairwise_distances
from ._optics_inner import quick_scan


def optics(X, min_samples=5, max_bound=np.inf, metric='euclidean',
           algorithm='ball_tree', leaf_size=30, p=2,
           metric_params=None, maxima_ratio=.75,
           rejection_ratio=.7, similarity_threshold=0.4,
           significant_min=.003, min_cluster_size_ratio=.005,
           min_maxima_ratio=0.001, n_jobs=1):
    """Perform OPTICS clustering from vector array

    OPTICS: Ordering Points To Identify the Clustering Structure
    Equivalent to DBSCAN, finds core sample of high density and expands
    clusters from them. Unlike DBSCAN, keeps cluster hierarchy for varying
    epsilon values. Optimized for usage on large point datasets.

    Parameters
    ----------
    X : array, shape (n_samples, n_features)
        The data.

    min_samples : int
        The number of samples in a neighborhood for a point to be considered
        as a core point.

    max_bound : float, optional
        The maximum distance between two samples for them to be considered
        as in the same neighborhood. This is also the largest object size
        expected within the dataset. Default value of "np.inf" will identify
        clusters across all scales; reducing `max_bound` will result in
        shorter run times.

    metric : string or callable, optional
        The distance metric to use for neighborhood lookups. Default is
        "minkowski". Other options include “euclidean”, “manhattan”,
        “chebyshev”, “haversine”, “seuclidean”, “hamming”, “canberra”,
        and “braycurtis”. The “wminkowski” and “mahalanobis” metrics are
        also valid with an additional argument.

    algorithm : {'auto', 'ball_tree', 'kd_tree', 'brute'}, optional
        Algorithm used to compute the nearest neighbors:
        - 'ball_tree' will use :class:`BallTree`
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

    p : integer, optional (default=2)
        Parameter for the Minkowski metric from
        :ref:`sklearn.metrics.pairwise.pairwise_distances`. When p = 1, this is
        equivalent to using manhattan_distance (l1), and euclidean_distance
        (l2) for p = 2. For arbitrary p, minkowski_distance (l_p) is used.

    metric_params : dict, optional (default=None)
        Additional keyword arguments for the metric function.

    n_jobs : int, optional (default=1)
        The number of parallel jobs to run for neighbors search.
        If ``-1``, then the number of jobs is set to the number of CPU cores.
        Affects only :meth:`kneighbors` and :meth:`kneighbors_graph` methods.

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
        result in more points being classified.

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
        Sets a lower threshold on how small a significant maxima can be

    min_cluster_size_ratio : float, optional
        Minimum percentage of dataset expected for cluster membership.

    min_maxima_ratio : float, optional
        Used to determine neighborhood size for minimum cluster membership.

    Returns
    -------
    core_sample_indices_ : array, shape (n_core_samples,)
        The indices of the core samples.

    labels_ : array, shape (n_samples,)
        The estimated labels.

    References
    ----------
    Ankerst, Mihael, Markus M. Breunig, Hans-Peter Kriegel, and Jörg Sander.
    "OPTICS: ordering points to identify the clustering structure." ACM SIGMOD
    Record 28, no. 2 (1999): 49-60.
    """

    clust = OPTICS.fit_predict(min_samples, max_bound, metric, algorithm,
                               leaf_size, p, metric_params, maxima_ratio,
                               rejection_ratio, similarity_threshold,
                               significant_min, min_cluster_size_ratio,
                               min_maxima_ratio, n_jobs)
    return clust.core_sample_indices_, clust.labels_


class OPTICS(NeighborsBase, KNeighborsMixin,
             RadiusNeighborsMixin, ClusterMixin, UnsupervisedMixin):
    """Estimate clustering structure from vector array

    OPTICS: Ordering Points To Identify the Clustering Structure
    Equivalent to DBSCAN, finds core sample of high density and expands
    clusters from them. Unlike DBSCAN, keeps cluster hierarchy for varying
    epsilon values. Optimized for usage on large point datasets.

    Parameters
    ----------
    min_samples : int
        The number of samples in a neighborhood for a point to be considered
        as a core point.

    max_bound : float, optional
        The maximum distance between two samples for them to be considered
        as in the same neighborhood. This is also the largest object size
        expected within the dataset. Default value of "np.inf" will identify
        clusters across all scales; reducing `max_bound` will result in
        shorter run times.

    metric : string or callable, optional
        The distance metric to use for neighborhood lookups. Default is
        "minkowski". Other options include “euclidean”, “manhattan”,
        “chebyshev”, “haversine”, “seuclidean”, “hamming”, “canberra”,
        and “braycurtis”. The “wminkowski” and “mahalanobis” metrics are
        also valid with an additional argument.

    algorithm : {'auto', 'ball_tree', 'kd_tree', 'brute'}, optional
        Algorithm used to compute the nearest neighbors:
        - 'ball_tree' will use :class:`BallTree`
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

    p : integer, optional (default=2)
        Parameter for the Minkowski metric from
        :ref:`sklearn.metrics.pairwise.pairwise_distances`. When p = 1, this is
        equivalent to using manhattan_distance (l1), and euclidean_distance
        (l2) for p = 2. For arbitrary p, minkowski_distance (l_p) is used.

    metric_params : dict, optional (default=None)
        Additional keyword arguments for the metric function.

    n_jobs : int, optional (default=1)
        The number of parallel jobs to run for neighbors search.
        If ``-1``, then the number of jobs is set to the number of CPU cores.
        Affects only :meth:`kneighbors` and :meth:`kneighbors_graph` methods.

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
        result in more points being classified.

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
        Sets a lower threshold on how small a significant maxima can be

    min_cluster_size_ratio : float, optional
        Minimum percentage of dataset expected for cluster membership.

    min_maxima_ratio : float, optional
        Used to determine neighborhood size for minimum cluster membership.

    Attributes
    ----------
    `core_sample_indices_` : array, shape (n_core_samples,)
        Indices of core samples.

    `labels_` : array, shape (n_samples,)
        Cluster labels for each point in the dataset given to fit().
        Noisy samples are given the label -1.

    `reachability_` : array, shape (n_samples,)
        Reachability distances per sample.

    `ordering_` : array, shape (n_samples,)
        The cluster ordered list of sample indices

    `core_dists_` : array, shape (n_samples,)
        Distance at which each sample becomes a core point.
        Points which will never be core have a distance of inf.

    References
    ----------
    Ankerst, Mihael, Markus M. Breunig, Hans-Peter Kriegel, and Jörg Sander.
    "OPTICS: ordering points to identify the clustering structure." ACM SIGMOD
    Record 28, no. 2 (1999): 49-60.
    """

    def __init__(self, min_samples=5, max_bound=np.inf, metric='euclidean',
                 algorithm='ball_tree', leaf_size=30, p=2,
                 metric_params=None, maxima_ratio=.75,
                 rejection_ratio=.7, similarity_threshold=0.4,
                 significant_min=.003, min_cluster_size_ratio=.005,
                 min_maxima_ratio=0.001, n_jobs=1):

        self._init_params(n_neighbors=min_samples,
                          algorithm=algorithm,
                          leaf_size=leaf_size, metric=metric, p=p,
                          metric_params=metric_params, n_jobs=n_jobs)

        self.max_bound = max_bound
        # min_samples is set via n_neighbors attribute in init params
        self.maxima_ratio = maxima_ratio
        self.rejection_ratio = rejection_ratio
        self.similarity_threshold = similarity_threshold
        self.significant_min = significant_min
        self.min_cluster_size_ratio = min_cluster_size_ratio
        self.min_maxima_ratio = min_maxima_ratio

    @property
    def n_neighbors(self):
        return self.min_samples

    @n_neighbors.setter
    def n_neighbors(self, value):
        self.min_samples = value

    def fit(self, X, y=None):
        """Perform OPTICS clustering

        Extracts an ordered list of points and reachability distances, and
        performs initial clustering using `max_bound` distance specified at
        OPTICS object instantiation.

        Parameters
        ----------
        X : array, shape (n_samples, n_features)
            The data.

        Returns
        -------
        self : instance of OPTICS
            The instance.
        """
        X = check_array(X, dtype=np.float)

        n_samples = len(X)
        # Start all points as 'unprocessed' ##
        self._processed = np.zeros((n_samples, 1), dtype=bool)
        self.reachability_ = np.ones(n_samples) * np.inf
        self.core_dists_ = np.ones(n_samples) * np.nan
        # Start all points as noise ##
        self.labels_ = -np.ones(n_samples, dtype=int)
        self.ordering_ = []

        # Check for valid n_samples relative to min_samples
        if self.min_samples > n_samples:
            raise ValueError("Number of training samples (%d) must be greater "
                             "than min_samples (%d) used for clustering." %
                             (n_samples, self.min_samples))

        super(OPTICS, self).fit(X)

        self.core_dists_[:] = self.kneighbors(X, self.min_samples)[0][:, -1]

        # Main OPTICS loop. Not parallelizable. The order that entries are
        # written to the 'ordering_' list is important!
        for point in range(n_samples):
            if not self._processed[point]:
                self._expand_cluster_order(point, X)

        indices_, self.labels_ = _extract_optics(self.ordering_,
                                                 self.reachability_,
                                                 self.maxima_ratio,
                                                 self.rejection_ratio,
                                                 self.similarity_threshold,
                                                 self.significant_min,
                                                 self.min_cluster_size_ratio,
                                                 self.min_maxima_ratio)
        self.core_sample_indices_ = indices_
        self.n_clusters_ = np.max(self.labels_)
        return self

    # OPTICS helper functions; these should not be public #

    def _expand_cluster_order(self, point, X):
        # As above, not parallelizable. Parallelizing would allow items in
        # the 'unprocessed' list to switch to 'processed'
        if self.core_dists_[point] <= self.max_bound:
            while not self._processed[point]:
                self._processed[point] = True
                self.ordering_.append(point)
                point = self._set_reach_dist(point, X)
        else:
            self.ordering_.append(point)  # For very noisy points
            self._processed[point] = True

    def _set_reach_dist(self, point_index, X):
        P = np.array(X[point_index]).reshape(1, -1)
        indices = self.radius_neighbors(P, radius=self.max_bound,
                                        return_distance=False)[0]

        # Checks to see if there more than one member in the neighborhood
        if iterable(indices):
            # Masking processed values; n_pr is 'not processed'
            n_pr = np.compress((np.take(self._processed,
                                        indices, axis=0) < 1).ravel(),
                               indices, axis=0)
            # n_pr = indices[(self._processed[indices] < 1).ravel()]
            # Keep n_jobs = 1 in the following lines...please
            if len(n_pr) > 0:
                dists = pairwise_distances(P, np.take(X, n_pr, axis=0),
                                           self.metric, n_jobs=1).ravel()

                rdists = np.maximum(dists, self.core_dists_[point_index])
                new_reach = np.minimum(np.take(self.reachability_,
                                               n_pr, axis=0), rdists)
                self.reachability_[n_pr] = new_reach

            # Checks to see if everything is already processed;
            # if so, return control to main loop
            if n_pr.size > 0:
                # Define return order based on reachability distance
                return(n_pr[quick_scan(np.take(self.reachability_,
                                               n_pr, axis=0), dists)])
            else:
                return point_index

    def extract_dbscan(self, eps):
        """Performs DBSCAN extraction for an arbitrary epsilon.

        Extraction runs in linear time. Note that if the `max_bound` OPTICS
        parameter was set to < inf for extracting reachability and ordering
        arrays, DBSCAN extractions will be unstable for `eps` values close to
        `max_bound`. Setting `eps` < (`max_bound` / 5.0) will guarantee
        extraction parity with DBSCAN.

        Parameters
        ----------
        eps : float or int, required
            DBSCAN `eps` parameter. Must be set to < `max_bound`. Equivalence
            with DBSCAN algorithm is achieved if `eps` is < (`max_bound` / 5)

        Returns
        -------
        core_sample_indices_ : array, shape (n_core_samples,)
            The indices of the core samples.

        labels_ : array, shape (n_samples,)
            The estimated labels.
        """
        check_is_fitted(self, 'reachability_')

        if eps > self.max_bound:
            raise ValueError('Specify an epsilon smaller than %s. Got %s.'
                             % (self.max_bound, eps))

        if eps * 5.0 > (self.max_bound * 1.05):
            warnings.warn(
                "Warning, max_bound (%s) is close to eps (%s): "
                "Output may be unstable." % (self.max_bound, eps),
                RuntimeWarning, stacklevel=2)
        # Stability warning is documented in _extract_dbscan method...

        return _extract_dbscan(self.ordering_, self.core_dists_,
                               self.reachability_, eps)


def _extract_dbscan(ordering, core_dists, reachability, eps):
    """Performs DBSCAN extraction for an arbitrary epsilon (`eps`).

    Parameters
    ----------
    ordering : array, shape (n_samples,)
        OPTICS ordered point indices (`ordering_`)
    core_dists : array, shape (n_samples,)
        Distances at which points become core (`core_dists_`)
    reachability : array, shape (n_samples,)
        Reachability distances calculated by OPTICS (`reachability_`)
    epsilon_prime : float or int
        DBSCAN `eps` parameter

    Returns
    -------
    core_sample_indices_ : array, shape (n_core_samples,)
        The indices of the core samples.

    labels_ : array, shape (n_samples,)
        The estimated labels.
    """

    n_samples = len(core_dists)
    is_core = np.zeros(n_samples, dtype=bool)
    labels = np.ones(n_samples, dtype=int) * -1
    cluster_id = 0

    for entry in ordering:
        if reachability[entry] > eps:
            if core_dists[entry] <= eps:
                cluster_id += 1
                labels[entry] = cluster_id
            else:
                # This is only needed for compatibility for repeated scans.
                labels[entry] = -1   # Noise points
                is_core[entry] = 0
        else:
            labels[entry] = cluster_id
            if core_dists[entry] <= eps:
                is_core[entry] = 1   # True
            else:
                is_core[entry] = 0   # False
    return np.arange(n_samples)[is_core], labels


def _extract_optics(ordering, reachability, maxima_ratio=.75,
                    rejection_ratio=.7, similarity_threshold=0.4,
                    significant_min=.003, min_cluster_size_ratio=.005,
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
        result in more points being classified.

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
        Sets a lower threshold on how small a significant maxima can be

    min_cluster_size_ratio : float, optional
        Minimum percentage of dataset expected for cluster membership.

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
    reachability_plot = reachability[ordering].tolist()
    root_node = _automatic_cluster(reachability_plot, ordering,
                                   maxima_ratio, rejection_ratio,
                                   similarity_threshold, significant_min,
                                   min_cluster_size_ratio, min_maxima_ratio)
    leaves = _get_leaves(root_node, [])
    # Start cluster id's at 1
    clustid = 1
    n_samples = len(reachability)
    is_core = np.zeros(n_samples, dtype=bool)
    labels = np.ones(n_samples, dtype=int) * -1
    # Start all points as non-core noise
    for leaf in leaves:
        index = ordering[leaf.start:leaf.end]
        labels[index] = clustid
        is_core[index] = 1
        clustid += 1
    return(np.arange(n_samples)[is_core], labels)


def _automatic_cluster(reachability_plot, reachability_ordering,
                       maxima_ratio, rejection_ratio,
                       similarity_threshold, significant_min,
                       min_cluster_size_ratio, min_maxima_ratio):

    min_neighborhood_size = 2
    min_cluster_size = int(min_cluster_size_ratio *
                           len(reachability_ordering))
    neighborhood_size = int(min_maxima_ratio * len(reachability_ordering))

    # Should this check for < min_samples? Should this be public?
    if min_cluster_size < 5:
        min_cluster_size = 5

    # Again, should this check < min_samples, should the parameter be public?
    if neighborhood_size < min_neighborhood_size:
        neighborhood_size = min_neighborhood_size

    local_maxima_points = _find_local_maxima(reachability_plot,
                                             reachability_ordering,
                                             neighborhood_size)
    root_node = _TreeNode(reachability_ordering, 0,
                          len(reachability_ordering), None)
    _cluster_tree(root_node, None, local_maxima_points,
                  reachability_plot, reachability_ordering, min_cluster_size,
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

    def assign_split_point(self, split_point):
        self.split_point = split_point

    def add_child(self, child):
        self.children.append(child)


def _is_local_maxima(index, reachability_plot, reachability_ordering,
                     neighborhood_size):
    # 0 = point at index is not local maxima
    # 1 = point at index is local maxima
    for i in range(1, neighborhood_size + 1):
        # process objects to the right of index
        if index + i < len(reachability_plot):
            if (reachability_plot[index] < reachability_plot[index + i]):
                return 0
        # process objects to the left of index
        if index - i >= 0:
            if (reachability_plot[index] < reachability_plot[index - i]):
                return 0
    return 1


def _find_local_maxima(reachability_plot, reachability_ordering,
                       neighborhood_size):
    local_maxima_points = {}
    # 1st and last points on Reachability Plot are not taken
    # as local maxima points
    for i in range(1, len(reachability_ordering) - 1):
        # if the point is a local maxima on the reachability plot with
        # regard to neighborhood_size, insert it into priority queue and
        # maxima list
        if (reachability_plot[i] > reachability_plot[i - 1] and
            reachability_plot[i] >= reachability_plot[i + 1] and
            _is_local_maxima(i, reachability_plot, reachability_ordering,
                             neighborhood_size) == 1):
            local_maxima_points[i] = reachability_plot[i]

    return sorted(local_maxima_points,
                  key=local_maxima_points.__getitem__, reverse=True)


def _cluster_tree(node, parent_node, local_maxima_points,
                  reachability_plot, reachability_ordering,
                  min_cluster_size, maxima_ratio, rejection_ratio,
                  similarity_threshold, significant_min):
    # node is a node or the root of the tree in the first call
    # parent_node is parent node of N or None if node is root of the tree
    # local_maxima_points is list of local maxima points sorted in
    # descending order of reachability
    if len(local_maxima_points) == 0:
        return  # parent_node is a leaf

    # take largest local maximum as possible separation between clusters
    s = local_maxima_points[0]
    node.assign_split_point(s)
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
        node.assign_split_point(-1)
        # if split_point is not significant, ignore this split and continue
        _cluster_tree(node, parent_node, local_maxima_points,
                      reachability_plot, reachability_ordering,
                      min_cluster_size, maxima_ratio, rejection_ratio,
                      similarity_threshold, significant_min)
        return

    # only check a certain ratio of points in the child
    # nodes formed to the left and right of the maxima
    # ...should check_ratio be a user settable parameter?
    check_ratio = .8
    check_value_1 = int(np.round(check_ratio * len(node_1.points)))
    check_value_2 = int(np.round(check_ratio * len(node_2.points)))
    if check_value_2 == 0:
        check_value_2 = 1
    avg_reach_value_1 = float(np.average(reachability_plot[(node_1.end -
                                         check_value_1):node_1.end]))
    avg_reach_value_2 = float(np.average(reachability_plot[node_2.start:(
                                node_2.start + check_value_2)]))

    if (float(avg_reach_value_1 / float(reachability_plot[s])) >
        maxima_ratio or float(avg_reach_value_2 /
                              float(reachability_plot[s])) > maxima_ratio):

        if float(avg_reach_value_1 /
                 float(reachability_plot[s])) < rejection_ratio:
            # reject node 2
            node_list.remove((node_2, local_max_2))
        if float(avg_reach_value_2 /
                 float(reachability_plot[s])) < rejection_ratio:
            # reject node 1
            node_list.remove((node_1, local_max_1))
        if (float(avg_reach_value_1 /
            float(reachability_plot[s])) >= rejection_ratio and float(
                avg_reach_value_2 /
                float(reachability_plot[s])) >= rejection_ratio):
            node.assign_split_point(-1)
            # since split_point is not significant,
            # ignore this split and continue (reject both child nodes)
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
            node_list.count((node_2, local_max_1)) > 0):
        # cluster 2 is too small
        node_list.remove((node_2, local_max_2))
    if len(node_list) == 0:
        # parent_node will be a leaf
        node.assign_split_point(-1)
        return

    # Check if nodes can be moved up one level - the new cluster created
    # is too "similar" to its parent, given the similarity threshold.
    bypass_node = 0
    if parent_node is not None:
        if (float(float(node.end - node.start) /
                  float(parent_node.end - parent_node.start)) >
                similarity_threshold):

            parent_node.children.remove(node)
            bypass_node = 1

    for nl in node_list:
        if bypass_node == 1:
            parent_node.add_child(nl[0])
            _cluster_tree(nl[0], parent_node, nl[1],
                          reachability_plot, reachability_ordering,
                          min_cluster_size, maxima_ratio, rejection_ratio,
                          similarity_threshold, significant_min)
        else:
            node.add_child(nl[0])
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
