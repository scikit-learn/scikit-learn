# -*- coding: utf-8 -*-
"""
DBSCAN: Density-Based Spatial Clustering of Applications with Noise
"""

# Author: Robert Layton <robertlayton@gmail.com>
#
# License: BSD

import numpy as np

from ..base import BaseEstimator
from . import calculate_similarity


def dbscan(X, eps=0.5, min_points=5, metric='euclidean',
           index_order=None, verbose=False):
    """Perform DBSCAN clustering from vector array or similarity matrix.

    Parameters
    ----------
    X: array [n_points, n_points] or [n_points, n_features]
        Array of similarities between points, or a feature array.
        The array is treated as a feature array unless the metric is given as
        'precomputed'.
    eps: float, optional
        The minimum similarity for two points to be considered
        in the same neighborhood.
    min_points: int, optional
        The number of points in a neighborhood for a point to be considered
        as a core point.
    metric: string, or callable
        The metric to use when calculating distance between instances in a
        feature array. If metric is a string, it must be one of the options
        allowed by scipy.spatial.distance.pdist for its metric parameter.
        If metric is "precomputed", X is assumed to be a similarity matrix and
        must be square.
        Alternatively, if metric is a callable function, it is called on each
        pair of instances (rows) and the resulting value recorded.
    index_order: [n_points] or None
        Order to observe points for clustering.
        If None, a random order is given.
        To look at points in order, use range(n).
    verbose: boolean, optional
        The verbosity level

    Returns
    -------
    core_points: array [n_core_points]
        Indices of core points.

    labels : array [n_points]
        Cluster labels for each point.  Noisy points are given the label -1.

    Notes
    -----
    See examples/plot_dbscan.py for an example.

    Reference:
    Ester, M., H. P. Kriegel, J. Sander, and X. Xu, “A Density-Based
    Algorithm for Discovering Clusters in Large Spatial Databases with Noise”.
    In: Proceedings of the 2nd International Conference on Knowledge Discovery
    and Data Mining, Portland, OR, AAAI Press, pp. 226–231. 1996
    """

    n = X.shape[0]
    # If index order not given, create random order.
    if index_order is None:
        index_order = np.arange(n)
        np.random.shuffle(index_order)
    assert len(index_order) == n, ("Index order must be of length n"
                                   " (%d expected, %d given)"
                                   % (n, len(index_order)))
    S = calculate_similarity(X, metric=metric)
    # Calculate neighborhood for all points. This leaves the original point
    # in, which needs to be considered later (i.e. point i is the
    # neighborhood of point i. While True, its useless information)
    neighborhoods = [np.where(x >= eps)[0] for x in S]
    # Initially, all points are noise.
    labels = np.array([-1] * n)
    # A list of all core points found.
    core_points = []
    # label_num is the label given to the new cluster
    label_num = 0
    # Look at all points and determine if they are core.
    # If they are then build a new cluster from them.
    for index in index_order:
        if labels[index] != -1 or len(neighborhoods[index]) < min_points:
            # This point is already classified, or not enough for a core point.
            continue
        core_points.append(index)
        labels[index] = label_num
        # candidates for new core points in the cluster.
        candidates = [index]
        while len(candidates) > 0:
            new_candidates = []
            # A candidate is a core point in the current cluster that has
            # not yet been used to expand the current cluster.
            for c in candidates:
                noise = np.where(labels[neighborhoods[c]] == -1)[0]
                noise = neighborhoods[c][noise]
                labels[noise] = label_num
                for neighbor in noise:
                    # check if its a core point as well
                    if len(neighborhoods[neighbor]) >= min_points:
                        # is new core point
                        new_candidates.append(neighbor)
                        core_points.append(neighbor)
            # Update candidates for next round of cluster expansion.
            candidates = new_candidates
        # Current cluster finished.
        # Next core point found will start a new cluster.
        label_num += 1
    return core_points, labels


class DBSCAN(BaseEstimator):
    """Perform DBSCAN clustering from vector array or similarity matrix.

    DBSCAN - Density-Based Spatial Clustering of Applications with Noise.
    Finds core points of high density and expands clusters from them.
    Good for data which contains clusters of similar density.

    Parameters
    ----------
    eps: float, optional
        The distance for two points to be considered in the same neighborhood
    min_points: int, optional
        The number of points in a neighborhood for a point to be considered
        as a core point.
    metric: string, or callable
        The metric to use when calculating distance between instances in a
        feature array. If metric is a string, it must be one of the options
        allowed by scipy.spatial.distance.pdist for its metric parameter.
        If metric is "precomputed", X is assumed to be a similarity matrix and
        must be square.
        Alternatively, if metric is a callable function, it is called on each
        pair of instances (rows) and the resulting value recorded.
    index_order: [n_points] or None
        Order to observe points for clustering.
        If None, a random order is given.
        To look at points in order, use range(n).
    verbose: boolean, optional
        The verbosity level

    Methods
    -------
    fit:
        Compute the clustering

    Attributes
    ----------
    core_points: array, shape = [n_core_points]
        Indices of core points.

    labels_ : array, shape = [n_points]
        Cluster labels for each point. Noisy points are given the label -1.

    Notes
    -----
    See examples/plot_dbscan.py for an example.

    Reference:
    Ester, M., H. P. Kriegel, J. Sander, and X. Xu, “A Density-Based
    Algorithm for Discovering Clusters in Large Spatial Databases with Noise”.
    In: Proceedings of the 2nd International Conference on Knowledge Discovery
    and Data Mining, Portland, OR, AAAI Press, pp. 226–231. 1996
    """

    def __init__(self, eps=0.5, min_points=5, metric='euclidean',
                 verbose=False, index_order=None):
        self.eps = eps
        self.min_points = min_points
        self.metric = metric
        self.verbose = verbose
        self.index_order = index_order
        self.verbose = verbose

    def fit(self, X, **params):
        """Compute DBSCAN labels for points, using similarity array S.

        Parameters
        ----------
        X: array [n_points, n_points] or [n_points, n_features]
            Array of similarities between points, or a feature array.
            The array is treated as a feature array unless the metric is
            given as 'precomputed'.
        params: dict
            Overwrite keywords from __init__.
        """

        self._set_params(**params)
        self.core_points_, self.labels_ = dbscan(X, **self._get_params())
        return self
