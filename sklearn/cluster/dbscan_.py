# -*- coding: utf-8 -*-
"""
DBSCAN: Density-Based Spatial Clustering of Applications with Noise
"""

# Author: Robert Layton <robertlayton@gmail.com>
#
# License: BSD 3 clause

import numpy as np

from ..base import BaseEstimator, ClusterMixin
from ..metrics import pairwise_distances
from ..utils import check_random_state
from ..neighbors import NearestNeighbors


def dbscan(X, eps=0.5, min_samples=5, metric='minkowski',
           algorithm='auto', leaf_size=30, p=2, random_state=None):
    """Perform DBSCAN clustering from vector array or distance matrix.

    Parameters
    ----------
    X: array [n_samples, n_samples] or [n_samples, n_features]
        Array of distances between samples, or a feature array.
        The array is treated as a feature array unless the metric is given as
        'precomputed'.

    eps: float, optional
        The maximum distance between two samples for them to be considered
        as in the same neighborhood.

    min_samples: int, optional
        The number of samples in a neighborhood for a point to be considered
        as a core point.

    metric: string, or callable
        The metric to use when calculating distance between instances in a
        feature array. If metric is a string or callable, it must be one of
        the options allowed by metrics.pairwise.pairwise_distances for its
        metric parameter.
        If metric is "precomputed", X is assumed to be a distance matrix and
        must be square.

    algorithm: {'auto', 'ball_tree', 'kd_tree', 'brute'}, optional
        The algorithm to be used by the NearestNeighbors module
        to compute pointwise distances and find nearest neighbors.
        See NearestNeighbors module documentation for details.

    leaf_size: int, optional (default = 30)
        Leaf size passed to BallTree or cKDTree. This can affect the speed
        of the construction and query, as well as the memory required
        to store the tree. The optimal value depends
        on the nature of the problem.

    p: float, optional
        The power of the Minkowski metric to be used to calculate distance
        between points.

    random_state: numpy.RandomState, optional
        The generator used to initialize the centers. Defaults to numpy.random.

    Returns
    -------
    core_samples: array [n_core_samples]
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

    X = np.asarray(X)
    n = X.shape[0]

    # If index order not given, create random order.
    random_state = check_random_state(random_state)
    index_order = random_state.permutation(n)

    # check for known metric powers
    distance_matrix = True
    if metric == 'precomputed':
        D = pairwise_distances(X, metric=metric)
    else:
        distance_matrix = False
        neighbors_model = NearestNeighbors(radius=eps, algorithm=algorithm,
                                           leaf_size=leaf_size,
                                           metric=metric, p=p)
        neighbors_model.fit(X)

    # Calculate neighborhood for all samples. This leaves the original point
    # in, which needs to be considered later (i.e. point i is the
    # neighborhood of point i. While True, its useless information)
    neighborhoods = []
    if distance_matrix:
        neighborhoods = [np.where(x <= eps)[0] for x in D]

    # Initially, all samples are noise.
    labels = -np.ones(n, dtype=np.int)

    # A list of all core samples found.
    core_samples = []

    # label_num is the label given to the new cluster
    label_num = 0

    # Look at all samples and determine if they are core.
    # If they are then build a new cluster from them.
    for index in index_order:
        # Already classified
        if labels[index] != -1:
            continue

        # get neighbors from neighborhoods or ballTree
        index_neighborhood = []
        if distance_matrix:
            index_neighborhood = neighborhoods[index]
        else:
            index_neighborhood = neighbors_model.radius_neighbors(
                X[index], eps, return_distance=False)[0]

        # Too few samples to be core
        if len(index_neighborhood) < min_samples:
            continue

        core_samples.append(index)
        labels[index] = label_num
        # candidates for new core samples in the cluster.
        candidates = [index]

        while len(candidates) > 0:
            new_candidates = []
            # A candidate is a core point in the current cluster that has
            # not yet been used to expand the current cluster.
            for c in candidates:
                c_neighborhood = []
                if distance_matrix:
                    c_neighborhood = neighborhoods[c]
                else:
                    c_neighborhood = neighbors_model.radius_neighbors(
                        X[c], eps, return_distance=False)[0]
                noise = np.where(labels[c_neighborhood] == -1)[0]
                noise = c_neighborhood[noise]
                labels[noise] = label_num
                for neighbor in noise:
                    n_neighborhood = []
                    if distance_matrix:
                        n_neighborhood = neighborhoods[neighbor]
                    else:
                        n_neighborhood = neighbors_model.radius_neighbors(
                            X[neighbor], eps, return_distance=False)[0]
                    # check if its a core point as well
                    if len(n_neighborhood) >= min_samples:
                        # is new core point
                        new_candidates.append(neighbor)
                        core_samples.append(neighbor)
            # Update candidates for next round of cluster expansion.
            candidates = new_candidates
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
        The number of samples in a neighborhood for a point to be considered
        as a core point.
    metric : string, or callable
        The metric to use when calculating distance between instances in a
        feature array. If metric is a string or callable, it must be one of
        the options allowed by metrics.pairwise.calculate_distance for its
        metric parameter.
        If metric is "precomputed", X is assumed to be a distance matrix and
        must be square.
    random_state : numpy.RandomState, optional
        The generator used to initialize the centers. Defaults to numpy.random.

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
    See examples/plot_dbscan.py for an example.

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

    def fit(self, X):
        """Perform DBSCAN clustering from features or distance matrix.

        Parameters
        ----------
        X: array [n_samples, n_samples] or [n_samples, n_features]
            Array of distances between samples, or a feature array.
            The array is treated as a feature array unless the metric is
            given as 'precomputed'.
        params: dict
            Overwrite keywords from __init__.
        """
        X = np.asarray(X)
        clust = dbscan(X, **self.get_params())
        self.core_sample_indices_, self.labels_ = clust
        self.components_ = X[self.core_sample_indices_].copy()
        return self
