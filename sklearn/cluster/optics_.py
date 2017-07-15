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
from ..neighbors import BallTree
from ..base import BaseEstimator, ClusterMixin
from ..metrics.pairwise import pairwise_distances
from ._optics_inner import quick_scan


class _SetOfObjects(BallTree):
    """Build balltree data structure with processing index.

    It's built from given data in preparation for OPTICS Algorithm

    Parameters
    ----------
    X: array, shape (n_samples, n_features)
        Input data.

    metric : str
        The metric used.

    **kwargs:
        kwargs passed to BallTree init.
    """

    def __init__(self, X, metric, **kwargs):

        super(_SetOfObjects, self).__init__(X,
                                            metric=metric, **kwargs)

        self._n = len(self.data)
        self.metric = metric
        # Start all points as 'unprocessed' ##
        self._processed = np.zeros((self._n, 1), dtype=bool)
        self.reachability_ = np.ones(self._n) * np.inf
        self.core_dists_ = np.ones(self._n) * np.nan
        # Start all points as noise ##
        self._cluster_id = -np.ones(self._n, dtype=int)
        self._is_core = np.zeros(self._n, dtype=bool)
        self.ordering_ = []


# OPTICS helper functions; these should not be public #
def _prep_optics(self, min_samples):
    self.core_dists_[:] = self.query(self.get_arrays()[0],
                                     k=min_samples)[0][:, -1]


def _build_optics(setofobjects, epsilon):
    # Main OPTICS loop. Not parallelizable. The order that entries are
    # written to the 'ordering_' list is important!
    for point in range(setofobjects._n):
        if not setofobjects._processed[point]:
            _expand_cluster_order(setofobjects, point, epsilon)


def _expand_cluster_order(setofobjects, point, epsilon):
    # As above, not parallelizable. Parallelizing would allow items in
    # the 'unprocessed' list to switch to 'processed'
    if setofobjects.core_dists_[point] <= epsilon:
        while not setofobjects._processed[point]:
            setofobjects._processed[point] = True
            setofobjects.ordering_.append(point)
            point = _set_reach_dist(setofobjects, point, epsilon)
    else:
        setofobjects.ordering_.append(point)  # For very noisy points
        setofobjects._processed[point] = True


def _set_reach_dist(setofobjects, point_index, epsilon):
    X = np.array(setofobjects.data[point_index]).reshape(1, -1)
    indices = setofobjects.query_radius(X, r=epsilon,
                                        return_distance=False,
                                        count_only=False,
                                        sort_results=False)[0]

    # Checks to see if there more than one member in the neighborhood
    if iterable(indices):
        # Masking processed values; n_pr is 'not processed'
        n_pr = np.compress((np.take(setofobjects._processed,
                                    indices, axis=0) < 1).ravel(),
                           indices, axis=0)
        # n_pr = indices[(setofobjects._processed[indices] < 1).ravel()]
        if len(n_pr) > 0:
            dists = pairwise_distances(X,
                                       np.take(setofobjects.get_arrays()[0],
                                               n_pr,
                                               axis=0),
                                       setofobjects.metric, n_jobs=1).ravel()

            rdists = np.maximum(dists, setofobjects.core_dists_[point_index])
            new_reach = np.minimum(np.take(setofobjects.reachability_,
                                           n_pr, axis=0),
                                   rdists)
            setofobjects.reachability_[n_pr] = new_reach

        # Checks to see if everything is already processed;
        # if so, return control to main loop
        if n_pr.size > 0:
            # Define return order based on reachability distance
            return(n_pr[quick_scan(np.take(setofobjects.reachability_,
                                           n_pr, axis=0), dists)])
        else:
            return point_index


class OPTICS(BaseEstimator, ClusterMixin):

    """Estimate clustering structure from vector array

    OPTICS: Ordering Points To Identify the Clustering Structure
    Equivalent to DBSCAN, finds core sample of high density and expands
    clusters from them. Unlike DBSCAN, keeps cluster hierarchy for varying
    epsilon values. Optimized for usage on large point datasets.

    Parameters
    ----------
    eps : float, optional
        The maximum distance between two samples for them to be considered
        as in the same neighborhood. This is also the largest object size
        expected within the dataset. Lower eps values can be used after
        OPTICS is run the first time, with fast returns of labels. Default
        value of "np.inf" will identify clusters across all scales; reducing
        eps will result in shorter run times.

    min_samples : int, optional
        The number of samples in a neighborhood for a point to be considered
        as a core point.

    metric : string or callable, optional
        The distance metric to use for neighborhood lookups. Default is
        "minkowski". Other options include “euclidean”, “manhattan”,
        “chebyshev”, “haversine”, “seuclidean”, “hamming”, “canberra”,
        and “braycurtis”. The “wminkowski” and “mahalanobis” metrics are
        also valid with an additional argument.

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

    def __init__(self, eps=np.inf, min_samples=5, metric='euclidean'):
        self.eps = eps
        self.min_samples = min_samples
        self.metric = metric

    def fit(self, X, y=None):
        """Perform OPTICS clustering

        Extracts an ordered list of points and reachability distances, and
        performs initial clustering using 'eps' distance specified at OPTICS
        object instantiation.

        Parameters
        ----------
        X : array, shape (n_samples, n_features)
            The data.

        Returns
        -------
        self : instance of OPTICS
            The instance.
        """
        X = check_array(X)

        n_samples = len(X)

        # Check for valid n_samples relative to min_samples
        if self.min_samples > n_samples:
            raise ValueError("Number of training samples (%d) must be greater "
                             "than min_samples (%d) used for clustering." %
                             (n_samples, self.min_samples))

        self.tree_ = _SetOfObjects(X, self.metric)
        _prep_optics(self.tree_, self.min_samples)
        _build_optics(self.tree_, self.eps * 5.0)
        self.extract_auto()
        self.core_sample_indices_ = np.arange(self._is_core)
        self.n_clusters_ = max(self._cluster_id)
        return self

    @property
    def reachability_(self):
        return self.tree_.reachability_[:]

    @property
    def core_dists_(self):
        return self.tree_.core_dists_[:]

    @property
    def _cluster_id(self):
        return self.tree_._cluster_id[:]

    @property
    def _is_core(self):
        return self.tree_._is_core[:]

    @property
    def _index(self):
        return self.tree_._index[:]

    @property
    def ordering_(self):
        return self.tree_.ordering_[:]

    @property
    def labels_(self):
        return self._cluster_id[:]

    def extract(self, epsilon_prime, clustering='auto', **kwargs):
        """Performs Automatic extraction for an arbitrary epsilon.

        It can also do DBSCAN equivalent and it can be run multiple times.

        See extract_auto() for full description of parameters

        Parameters
        ----------
        epsilon_prime: float or int, optional
            Used for 'dbscan' clustering. Must be less than or equal to what
            was used for prep and build steps

        clustering: {'auto', 'dbscan'}, optional
            Type of cluster extraction to perform; defaults to 'auto'.

        Returns
        -------
        core_sample_indices_: array, shape (n_core_samples,)
            The indices of the core samples.

        labels_ : array, shape (n_samples,)
            The estimated labels.
        """
        check_is_fitted(self, 'reachability_')

        if epsilon_prime > self.eps * 5.0:
            raise ValueError('Specify an epsilon smaller than %s. Got %s.'
                             % (self.eps * 5, epsilon_prime))

        if clustering == 'dbscan':
            self.eps_prime = epsilon_prime
            _extract_DBSCAN(self, epsilon_prime)
        elif clustering == 'auto':
            self.extract_auto(**kwargs)

        self.core_sample_indices_ = \
            np.arange(len(self.reachability_))[self._is_core == 1]
        self.n_clusters = max(self._cluster_id)

        if epsilon_prime > (self.eps * 1.05):
            warnings.warn(
                "Warning, eps (%s) is close to epsilon_prime (%s): "
                "Output may be unstable." % (self.eps, epsilon_prime),
                RuntimeWarning, stacklevel=2)
            # XXX explain 5% thing

        return self.core_sample_indices_, self.labels_

    def extract_auto(self,
                     maxima_ratio=.75,
                     rejection_ratio=.7,
                     similarity_threshold=0.4,
                     significant_min=.003,
                     min_cluster_size_ratio=.005,
                     min_maxima_ratio=0.001):
        """Performs automatic cluster extraction for variable density data.

        Can be run multiple times with adjusted parameters. Only returns
        core and noise labels.

        Parameters
        ----------
        maxima_ratio: float, optional
            The maximum ratio we allow of average height of clusters on the
            right and left to the local maxima in question. The higher the
            ratio, the more generous the algorithm is to preserving local
            minimums, and the more cuts the resulting tree will have.

        rejection_ratio: float, optional
            Adjusts the fitness of the clustering. When the maxima_ratio is
            exceeded, determine which of the clusters to the left and right to
            reject based on rejection_ratio

        similarity_threshold: float, optional
            Used to check if nodes can be moved up one level, that is, if the
            new cluster created is too "similar" to its parent, given the
            similarity threshold. Similarity can be determined by 1) the size
            of the new cluster relative to the size of the parent node or
            2) the average of the reachability values of the new cluster
            relative to the average of the reachability values of the parent
            node. A lower value for the similarity threshold means less levels
            in the tree.

        significant_min: float, optional
            Sets a lower threshold on how small a significant maxima can be
            min_cluster_size_ratio: float, optional
            Minimum percentage of dataset expected for cluster membership.

        min_maxima_ratio: float, optional
            Used to determine neighborhood size for minimum cluster membership.
        """
        check_is_fitted(self, 'reachability_')

        # Extraction wrapper
        RPlot = self.reachability_[self.ordering_].tolist()
        RPoints = self.ordering_
        root_node = _automatic_cluster(RPlot, RPoints, maxima_ratio,
                                       rejection_ratio, similarity_threshold,
                                       significant_min, min_cluster_size_ratio,
                                       min_maxima_ratio)
        leaves = _get_leaves(root_node, [])
        # Start cluster id's at 1
        clustid = 1
        # Start all points as non-core noise
        self._cluster_id[:] = -1
        self._is_core[:] = 0
        for leaf in leaves:
            index = self.ordering_[leaf.start:leaf.end]
            self._cluster_id[index] = clustid
            self._is_core[index] = 1
            clustid += 1


def _automatic_cluster(RPlot, RPoints, maxima_ratio, rejection_ratio,
                       similarity_threshold, significant_min,
                       min_cluster_size_ratio, min_maxima_ratio):

    min_neighborhood_size = 2
    min_cluster_size = int(min_cluster_size_ratio * len(RPoints))
    neighborhood_size = int(min_maxima_ratio * len(RPoints))

    # Should this check for < min_samples? Should this be public?
    if min_cluster_size < 5:
        min_cluster_size = 5

    # Again, should this check < min_samples, should the parameter be public?
    if neighborhood_size < min_neighborhood_size:
        neighborhood_size = min_neighborhood_size

    local_maxima_points = _find_local_maxima(RPlot, RPoints, neighborhood_size)
    root_node = _TreeNode(RPoints, 0, len(RPoints), None)
    _cluster_tree(root_node, None, local_maxima_points,
                  RPlot, RPoints, min_cluster_size,
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


def _is_local_maxima(index, RPlot, RPoints, neighborhood_size):
    # 0 = point at index is not local maxima
    # 1 = point at index is local maxima
    for i in range(1, neighborhood_size + 1):
        # process objects to the right of index
        if index + i < len(RPlot):
            if (RPlot[index] < RPlot[index + i]):
                return 0
        # process objects to the left of index
        if index - i >= 0:
            if (RPlot[index] < RPlot[index - i]):
                return 0
    return 1


def _find_local_maxima(RPlot, RPoints, neighborhood_size):
    local_maxima_points = {}
    # 1st and last points on Reachability Plot are not taken
    # as local maxima points
    for i in range(1, len(RPoints) - 1):
        # if the point is a local maxima on the reachability plot with
        # regard to neighborhood_size, insert it into priority queue and
        # maxima list
        if (RPlot[i] > RPlot[i - 1] and RPlot[i] >= RPlot[i + 1] and
                _is_local_maxima(i, RPlot, RPoints, neighborhood_size) == 1):
            local_maxima_points[i] = RPlot[i]

    return sorted(local_maxima_points,
                  key=local_maxima_points.__getitem__, reverse=True)


def _cluster_tree(node, parent_node, local_maxima_points, RPlot, RPoints,
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
    node_1 = _TreeNode(RPoints[node.start:s], node.start, s, node)
    node_2 = _TreeNode(RPoints[s + 1:node.end], s + 1, node.end, node)
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

    if RPlot[s] < significant_min:
        node.assign_split_point(-1)
        # if split_point is not significant, ignore this split and continue
        _cluster_tree(node, parent_node, local_maxima_points, RPlot, RPoints,
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
    avg_reach_value_1 = float(np.average(RPlot[(node_1.end -
                                                check_value_1):node_1.end]))
    avg_reach_value_2 = float(np.average(RPlot[node_2.start:(node_2.start +
                                                             check_value_2)]))

    if (float(avg_reach_value_1 / float(RPlot[s])) >
        maxima_ratio or float(avg_reach_value_2 /
                              float(RPlot[s])) > maxima_ratio):

        if float(avg_reach_value_1 / float(RPlot[s])) < rejection_ratio:
            # reject node 2
            node_list.remove((node_2, local_max_2))
        if float(avg_reach_value_2 / float(RPlot[s])) < rejection_ratio:
            # reject node 1
            node_list.remove((node_1, local_max_1))
        if (float(avg_reach_value_1 / float(RPlot[s])) >=
            rejection_ratio and float(avg_reach_value_2 /
                                      float(RPlot[s])) >= rejection_ratio):
            node.assign_split_point(-1)
            # since split_point is not significant,
            # ignore this split and continue (reject both child nodes)
            _cluster_tree(node, parent_node, local_maxima_points,
                          RPlot, RPoints, min_cluster_size,
                          maxima_ratio, rejection_ratio,
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
            _cluster_tree(nl[0], parent_node, nl[1], RPlot, RPoints,
                          min_cluster_size, maxima_ratio, rejection_ratio,
                          similarity_threshold, significant_min)
        else:
            node.add_child(nl[0])
            _cluster_tree(nl[0], node, nl[1], RPlot, RPoints, min_cluster_size,
                          maxima_ratio, rejection_ratio,
                          similarity_threshold, significant_min)


def _get_leaves(node, arr):
    if node is not None:
        if node.split_point == -1:
            arr.append(node)
        for n in node.children:
            _get_leaves(n, arr)
    return arr


def _extract_DBSCAN(setofobjects, epsilon_prime):
    # Extract DBSCAN Equivalent cluster structure
    # Important: Epsilon prime should be less than epsilon used in OPTICS
    cluster_id = 0
    for entry in setofobjects.ordering_:
        if setofobjects.reachability_[entry] > epsilon_prime:
            if setofobjects.core_dists_[entry] <= epsilon_prime:
                cluster_id += 1
                setofobjects._cluster_id[entry] = cluster_id
            else:
                # This is only needed for compatibility for repeated scans.
                setofobjects._cluster_id[entry] = -1   # Noise points
                setofobjects._is_core[entry] = 0
        else:
            setofobjects._cluster_id[entry] = cluster_id
            if setofobjects.core_dists_[entry] <= epsilon_prime:
                setofobjects._is_core[entry] = 1   # True
            else:
                setofobjects._is_core[entry] = 0   # False
