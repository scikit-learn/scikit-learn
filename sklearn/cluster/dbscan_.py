# -*- coding: utf-8 -*-
"""
DBSCAN: Density-Based Spatial Clustering of Applications with Noise
"""

# Author: Robert Layton <robertlayton@gmail.com>
#         Joel Nothman <joel.nothman@gmail.com>
#
# License: BSD 3 clause

import numpy as np

from ..base import BaseEstimator, ClusterMixin
from ..metrics import pairwise_distances
from ..utils import check_random_state, check_array, check_consistent_length
from ..neighbors import NearestNeighbors


def dbscan(X, eps=0.5, min_samples=5, metric='minkowski',
           algorithm='auto', leaf_size=30, p=2, sample_weight=None,
           random_state=None):
    """Perform DBSCAN clustering from vector array or distance matrix.

    Parameters
    ----------
    X : array or sparse (CSR) matrix of shape (n_samples, n_features), or \
            array of shape (n_samples, n_samples)
        A feature array, or array of distances between samples if
        ``metric='precomputed'``.

    eps : float, optional
        The maximum distance between two samples for them to be considered
        as in the same neighborhood.

    min_samples : int, optional
        The number of samples (or total weight) in a neighborhood for a point
        to be considered as a core point.

    metric : string, or callable
        The metric to use when calculating distance between instances in a
        feature array. If metric is a string or callable, it must be one of
        the options allowed by metrics.pairwise.pairwise_distances for its
        metric parameter.
        If metric is "precomputed", X is assumed to be a distance matrix and
        must be square.

    algorithm : {'auto', 'ball_tree', 'kd_tree', 'brute'}, optional
        The algorithm to be used by the NearestNeighbors module
        to compute pointwise distances and find nearest neighbors.
        See NearestNeighbors module documentation for details.

    leaf_size : int, optional (default = 30)
        Leaf size passed to BallTree or cKDTree. This can affect the speed
        of the construction and query, as well as the memory required
        to store the tree. The optimal value depends
        on the nature of the problem.

    p : float, optional
        The power of the Minkowski metric to be used to calculate distance
        between points.

    sample_weight : array, shape (n_samples,), optional
        Weight of each sample, such that a sample with weight greater
        than ``min_samples`` is automatically a core sample; a sample with
        negative weight may inhibit its eps-neighbor from being core.
        Note that weights are absolute, and default to 1.

    random_state: numpy.RandomState, optional
        The generator used to shuffle the samples. Defaults to numpy.random.

    Returns
    -------
    core_samples : array [n_core_samples]
        Indices of core samples.

    labels : array [n_samples]
        Cluster labels for each point.  Noisy samples are given the label -1.

    Notes
    -----
    See examples/cluster/plot_dbscan.py for an example.

    References
    ----------
    Ester, M., H. P. Kriegel, J. Sander, and X. Xu, "A Density-Based
    Algorithm for Discovering Clusters in Large Spatial Databases with Noise".
    In: Proceedings of the 2nd International Conference on Knowledge Discovery
    and Data Mining, Portland, OR, AAAI Press, pp. 226-231. 1996
    """
    if not eps > 0.0:
        raise ValueError("eps must be positive.")

    X = check_array(X, accept_sparse='csr')
    if sample_weight is not None:
        sample_weight = np.asarray(sample_weight)
        check_consistent_length(X, sample_weight)

    # If index order not given, create random order.
    random_state = check_random_state(random_state)

    # Calculate neighborhood for all samples. This leaves the original point
    # in, which needs to be considered later (i.e. point i is in the
    # neighborhood of point i. While True, its useless information)
    if metric == 'precomputed':
        D = pairwise_distances(X, metric=metric)
        neighborhoods = [np.where(x <= eps)[0] for x in D]
    else:
        neighbors_model = NearestNeighbors(radius=eps, algorithm=algorithm,
                                           leaf_size=leaf_size,
                                           metric=metric, p=p)
        neighbors_model.fit(X)
        neighborhoods = neighbors_model.radius_neighbors(X, eps,
                                                         return_distance=False)
        neighborhoods = np.array(neighborhoods)
    if sample_weight is None:
        n_neighbors = np.array([len(neighbors) for neighbors in neighborhoods])
    else:
        n_neighbors = np.array([np.sum(sample_weight[neighbors])
                                for neighbors in neighborhoods])

    # Initially, all samples are noise.
    labels = -np.ones(X.shape[0], dtype=np.int)

    # A list of all core samples found.
    core_samples = np.flatnonzero(n_neighbors > min_samples)
    index_order = core_samples[random_state.permutation(core_samples.shape[0])]

    # label_num is the label given to the new cluster
    label_num = 0

    # Look at all samples and determine if they are core.
    # If they are then build a new cluster from them.
    for index in index_order:
        # Already classified
        if labels[index] != -1:
            continue

        labels[index] = label_num

        # candidates for new core samples in the cluster.
        candidates = [index]
        while len(candidates) > 0:
            # The tolist() is needed for NumPy 1.6.
            cand_neighbors = np.concatenate(np.take(neighborhoods, candidates,
                                                    axis=0).tolist())
            cand_neighbors = np.unique(cand_neighbors)
            noise = cand_neighbors[labels.take(cand_neighbors) == -1]
            labels[noise] = label_num
            # A candidate is a core point in the current cluster that has
            # not yet been used to expand the current cluster.
            candidates = np.intersect1d(noise, core_samples,
                                        assume_unique=True)
        # Current cluster finished.
        # Next core point found will start a new cluster.
        label_num += 1
    return core_samples, labels


class DBSCAN(BaseEstimator, ClusterMixin):
    """Perform DBSCAN clustering from vector array or distance matrix.

    DBSCAN - Density-Based Spatial Clustering of Applications with Noise.
    Finds core samples of high density and expands clusters from them.
    Good for data which contains clusters of similar density.

    Parameters
    ----------
    eps : float, optional
        The maximum distance between two samples for them to be considered
        as in the same neighborhood.
    min_samples : int, optional
        The number of samples (or total weight) in a neighborhood for a point
        to be considered as a core point.
    metric : string, or callable
        The metric to use when calculating distance between instances in a
        feature array. If metric is a string or callable, it must be one of
        the options allowed by metrics.pairwise.calculate_distance for its
        metric parameter.
        If metric is "precomputed", X is assumed to be a distance matrix and
        must be square.
    random_state : numpy.RandomState, optional
        The generator used to shuffle the samples. Defaults to numpy.random.
    algorithm : {'auto', 'ball_tree', 'kd_tree', 'brute'}, optional
        The algorithm to be used by the NearestNeighbors module
        to compute pointwise distances and find nearest neighbors.
        See NearestNeighbors module documentation for details.
    leaf_size : int, optional (default = 30)
        Leaf size passed to BallTree or cKDTree. This can affect the speed
        of the construction and query, as well as the memory required
        to store the tree. The optimal value depends
        on the nature of the problem.
 
    Attributes
    ----------
    core_sample_indices_ : array, shape = [n_core_samples]
        Indices of core samples.

    components_ : array, shape = [n_core_samples, n_features]
        Copy of each core sample found by training.

    labels_ : array, shape = [n_samples]
        Cluster labels for each point in the dataset given to fit().
        Noisy samples are given the label -1.

    Notes
    -----
    See examples/cluster/plot_dbscan.py for an example.

    References
    ----------
    Ester, M., H. P. Kriegel, J. Sander, and X. Xu, "A Density-Based
    Algorithm for Discovering Clusters in Large Spatial Databases with Noise".
    In: Proceedings of the 2nd International Conference on Knowledge Discovery
    and Data Mining, Portland, OR, AAAI Press, pp. 226-231. 1996
    """

    def __init__(self, eps=0.5, min_samples=5, metric='euclidean',
                 algorithm='auto', leaf_size=30, p=None, random_state=None):
        self.eps = eps
        self.min_samples = min_samples
        self.metric = metric
        self.algorithm = algorithm
        self.leaf_size = leaf_size
        self.p = p
        self.random_state = random_state

    def fit(self, X, sample_weight=None):
        """Perform DBSCAN clustering from features or distance matrix.

        Parameters
        ----------
        X : array or sparse (CSR) matrix of shape (n_samples, n_features), or \
                array of shape (n_samples, n_samples)
            A feature array, or array of distances between samples if
            ``metric='precomputed'``.
        sample_weight : array, shape (n_samples,), optional
            Weight of each sample, such that a sample with weight greater
            than ``min_samples`` is automatically a core sample; a sample with
            negative weight may inhibit its eps-neighbor from being core.
            Note that weights are absolute, and default to 1.
        """
        X = check_array(X, accept_sparse='csr')
        clust = dbscan(X, sample_weight=sample_weight, **self.get_params())
        self.core_sample_indices_, self.labels_ = clust
        if len(self.core_sample_indices_):
            # fix for scipy sparse indexing issue
            self.components_ = X[self.core_sample_indices_].copy()
        else:
            # no core samples
            self.components_ = np.empty((0, X.shape[1]))
        return self

    def fit_predict(self, X, y=None, sample_weight=None):
        """Performs clustering on X and returns cluster labels.

        Parameters
        ----------
        X : array or sparse (CSR) matrix of shape (n_samples, n_features), or \
                array of shape (n_samples, n_samples)
            A feature array, or array of distances between samples if
            ``metric='precomputed'``.
        sample_weight : array, shape (n_samples,), optional
            Weight of each sample, such that a sample with weight greater
            than ``min_samples`` is automatically a core sample; a sample with
            negative weight may inhibit its eps-neighbor from being core.
            Note that weights are absolute, and default to 1.

        Returns
        -------
        y : ndarray, shape (n_samples,)
            cluster labels
        """
        self.fit(X, sample_weight=sample_weight)
        return self.labels_
