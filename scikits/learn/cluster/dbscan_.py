# -*- coding: utf-8 -*-
""" Algorithms for clustering : DBSCAN

    DBSCAN: (Density-Based Spatial Clustering of Applications with Noise)
"""
# Author: Robert Layton robertlayton@gmail.com
#
# License: BSD

import numpy as np
from scipy.spatial import distance

from ..base import BaseEstimator


def dbscan(S, eps=0.5, min_points=5, metric='euclidean',
           index_order=None, verbose=False, is_similarity=None,):
    """Perform DBSCAN Clustering of data

    Parameters
    ----------

    S: array [n_points, n_points] or [n_points, n_features]
        Matrix of similarities between points, or a feature matrix.
        If the matrix is square, it is treated as a similarity matrix,
        otherwise it is treated as a feature matrix. Use is_similarity to
        override this pattern.
    eps: float, optional
        The minimum similarity for two points to be considered
        in the same neighbourhood.
    min_points: int, optional
        The number of points in a neighbourhood for a point to be considered
        as a core point.
    metric: string, or callable
        The metric to use when calculating distance between instances in a
        feature matrix. If metric is a string, it must be one of the options
        allowed by scipy.spatial.distance.pdist for its metric parameter.
        Alternatively, if metric is a callable function, it is called on each
        pair of instances (rows) and the resulting value recorded.
    index_order: [n_points] or None
        Order to observe points for clustering.
        If None, a random order is given.
        To look at points in order, use range(n).
    verbose: boolean, optional
        The verbosity level
    is_similarity: boolean, optional (default=None)
        Overrides the behaviour of the matrix handling of S.
        If is_similarity is None, any square matrix is handled as a similarity
        matrix and any non-square matrix is a feature matrix.
        If is_similarity is True, any matrix is handled as a similarity matrix,
        and the procedure will raise a ValueError if the matrix is not square.
        If is_similarity is False, any matrix will be handled as a feature
        matrix, including square matrices.

    Returns
    -------

    core_points: array [n_core_points]
        index of core points

    labels : array [n_points]
        cluster labels for each point
        Noisey points are given the label -1

    Notes
    -----
    See examples/plot_dbscan.py for an example.

    Reference:
    Ester, M., H. P. Kriegel, J. Sander, and X. Xu, “A Density-Based
    Algorithm for Discovering Clusters in Large Spatial Databases with Noise”.
    In: Proceedings of the 2nd International Conference on Knowledge Discovery
    and Data Mining, Portland, OR, AAAI Press, pp. 226–231. 1996


    """
    n = S.shape[0]
    # If index order not given, create random order.
    if index_order is None:
        index_order = np.arange(n)
        np.random.shuffle(index_order)
    assert len(index_order) == n, ("Index order must be of length n"
                                   " (%d expected, %d given)"
                                   % (n, len(index_order)))
    S = calculateSimilarity(S, metric=metric, is_similarity=is_similarity)
    # Calculate neighbourhood for all points. This leaves the original point
    # in, which needs to be considered later (i.e. point i is the
    # neighbourhood of point i. While True, its useless information)
    neighbourhoods = [np.where(x >= eps)[0] for x in S]
    # Initially, all points are noise.
    labels = np.zeros((n,), dtype='int') - 1
    # A list of all core points found.
    core_points = []
    # Look at all points and determine if they are core.
    # If they are then build a new cluster from them.
    for index in index_order:
        if labels[index] != -1 or len(neighbourhoods[index]) < min_points:
            # This point is already classified, or not enough for a core point.
            continue
        core_points.append(index)
        # label_num is the label given to the new cluster
        label_num = np.max(labels) + 1
        labels[index] = label_num
        # candidates for new core points in the cluster.
        candidates = [index]
        while len(candidates) > 0:
            new_candidates = []
            # A candidate is a core point in the current cluster that has
            # not yet been used to expand the current cluster.
            for c in candidates:
                for neighbour in neighbourhoods[c]:
                    if labels[neighbour] == -1:
                        # neighbour is part of the current cluster iff
                        # it is not part of another cluster already.
                        labels[neighbour] = label_num
                        # check if its a core point as well
                        if len(neighbourhoods[neighbour]) >= min_points:
                            # is new core point
                            new_candidates.append(neighbour)
                            core_points.append(neighbour)
            # Update candidates for next round of cluster expansion.
            candidates = new_candidates
    return core_points, labels


def calculateSimilarity(S, metric=None, is_similarity=None):
    n, d = S.shape
    # If the matrix looks square, it may be a similarity matrix.
    if n == d:
        if is_similarity is None or is_similarity:
            return S
    else:
        # Matrix is not square, so it cannot be a similarity matrix.
        if is_similarity:
            raise ValueError("Matrix not square, "
                             "cannot be a similarity matrix."
                             " Size: %d x %d." % (n, d))
    # In all other cases, the matrix is to be considered as a feature matrix.
    D = distance.squareform(distance.pdist(S, metric=metric))
    S = 1. - (D / np.max(D))
    return S


###############################################################################

class DBSCAN(BaseEstimator):
    """Perform DBSCAN Clustering of data

    DBSCAN - Density-Based Spatial Clustering of Applications with Noise.
    Finds core points of high density and expands clusters from them.
    Good for data which contains clusters of similar density.

    Parameters
    ----------

    eps: float, optional
        The distance for two points to be considered in the same neighbourhood
    min_points: int, optional
        The number of points in a neighbourhood for a point to be considered
        as a core point.
    metric: string, or callable
        The metric to use when calculating distance between instances in a
        feature matrix. If metric is a string, it must be one of the options
        allowed by scipy.spatial.distance.pdist for its metric parameter.
        Alternatively, if metric is a callable function, it is called on each
        pair of instances (rows) and the resulting value recorded.
    index_order: [n_points] or None
        Order to observe points for clustering.
        If None, a random order is given.
        To look at points in order, use range(n).
    verbose: boolean, optional
        The verbosity level
    is_similarity: boolean, optional (default=None)
        Overrides the behaviour of the matrix handling of S.
        If is_similarity is None, any square matrix is handled as a similarity
        matrix and any non-square matrix is a feature matrix.
        If is_similarity is True, any matrix is handled as a similarity matrix,
        and the procedure will raise a ValueError if the matrix is not square.
        If is_similarity is False, any matrix will be handled as a feature
        matrix, including square matrices.

    Methods
    -------

    fit:
        Compute the clustering

    Attributes
    ----------

    core_points: array [n_core_points]
        index of core points

    labels : array [n_points]
        cluster labels for each point
        Noisey points are given the label -1


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
                 verbose=False, index_order=None,
                 is_similarity=None):
        self.eps = eps
        self.min_points = min_points
        self.metric = metric
        self.verbose = verbose
        self.index_order = index_order
        self.verbose = verbose
        self.is_similarity = is_similarity

    def fit(self, S, **params):
        """Compute DBSCAN labels for points, using similarity matrix S.

        Parameters
        ----------

        S: array [n_points, n_points] or [n_points, n_features]
            Matrix of similarities between points, or a feature matrix.
            If the matrix is square, it is treated as a similarity matrix,
            otherwise it is treated as a feature matrix. Use is_similarity to
            override this pattern.
        params: Overwrite keywords from __init__

        """
        self._set_params(**params)
        self.core_points_, self.labels_ = dbscan(S, eps=self.eps,
                                                 min_points=self.min_points,
                                                 verbose=self.verbose,
                                                 metric=self.metric,
                                                 index_order=self.index_order,
                                                 is_similarity=\
                                                 self.is_similarity)
        return self
