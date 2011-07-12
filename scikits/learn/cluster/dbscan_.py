# -*- coding: utf-8 -*-
""" Algorithms for clustering : Meanshift,  Affinity propagation and spectral
clustering.

"""
# Author: Robert Layton robertlayton@gmail.com
#
# License: BSD

import numpy as np

from ..base import BaseEstimator


def dbscan(S, eps=0.5, min_points=5, index_order=None,
           verbose=False):
    """Perform Affinity Propagation Clustering of data

    Parameters
    ----------

    S: array [n_points, n_points]
        Matrix of similarities between points
    eps: float, optional
        The minimum similarity for two points to be considered
        in the same neighbourhood.
    min_points: int, optional
        The number of points in a neighbourhood for a point to be considered
        as a core point.
    index_order: [n_points] or None
        Order to observe points for clustering. If None, a random order is given.
        To look at points in order, use range(n).
    verbose: boolean, optional
        The verbosity level

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
    and Data Mining, Portland, OR, AAAI Press, pp. 226–231. 2006


    """
    n = S.shape[0]
    # if index order not given, create random order
    if index_order is None:
        index_order = np.arange(n)
        np.random.shuffle(index_order)
    assert len(index_order) == n, ("Index order must be of length n"
                                   " (%d expected, %d given)"
                                   % (n, len(index_order)))                                   
    # Calculate neighbourhood for all points. This leaves the original point
    # in which needs to be considered later (i.e. point i is the the neighbourhood
    # of point i. While True, its useless information)
    neighbourhoods = [np.where(x >= eps)[0] for x in S]
    # Initially, all points are noise
    labels = np.zeros((n,), dtype='int') - 1
    # list of all core points found
    core_points = []
    # main DBSCAN loop: look at all points, consider if they are core, then
    # build a cluster from them if they are
    for index in index_order:
        if labels[index] != -1 or len(neighbourhoods[index]) < min_points:
            # point already classified, or not enough for a core point
            continue
        core_points.append(index)
        label_num = np.max(labels) + 1 # new cluster number
        labels[index] = label_num # start of new cluster
        candidates = [index,] # candidate points
        while len(candidates) > 0:
            new_candidates = []
            # a candidate is a core point in the current cluster that has
            # not yet been used to expand the current cluster
            for c in candidates:
                for neighbour in neighbourhoods[c]:
                    if labels[neighbour] == -1:
                        # neighbour is part of cluster iff its not part of another
                        # cluster already
                        labels[neighbour] = label_num
                        # check if its a core point as well
                        if len(neighbourhoods[neighbour]) >= min_points:
                            # is new core point
                            new_candidates.append(neighbour)
                            core_points.append(neighbour)
            # update candidates for next round of cluster expansion
            candidates = new_candidates
    return core_points, labels


###############################################################################

class DBSCAN(BaseEstimator):
    """Perform DBSCAN Clustering of data

    Parameters
    ----------

    eps: float, optional
        The distance for two points to be considered in the same neighbourhood
    min_points: int, optional
        The number of points in a neighbourhood for a point to be considered
        as a core point.

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
    and Data Mining, Portland, OR, AAAI Press, pp. 226–231. 2006
    
    """

    def __init__(self, eps=0.5, min_points=5, verbose=False):
        self.eps = eps
        self.min_points = min_points
        self.verbose = verbose


    def fit(self, S, **params):
        """compute DBSCAN

        Parameters
        ----------

        S: array [n_points, n_points]
            Matrix of similarities between points
        eps: float, optional
            The distance for two points to be considered in the same neighbourhood
        min_points: int, optional
            The number of points in a neighbourhood for a point to be considered
            as a core point.
        verbose: boolean, optional
            The verbosity level

        """
        self._set_params(**params)
        self.core_points_, self.labels_ = dbscan(S, eps=self.eps,
                                                 min_points=self.min_points,
                                                 verbose=self.verbose)
        return self
