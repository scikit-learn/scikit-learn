# -*- coding: utf-8 -*-
"""
DBSCAN: Density-Based Spatial Clustering of Applications with Noise
"""

# Author: Robert Layton <robertlayton@gmail.com>
#
# License: BSD

import numpy as np

from ..base import BaseEstimator
from ..metrics import pairwise_distances
from ..utils import check_random_state


def dbscan(X, eps=0.5, min_samples=5, metric='euclidean',
           random_state=None, verbose=False):
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
        the options allowed by metrics.pairwise.calculate_distance for its
        metric parameter.
        If metric is "precomputed", X is assumed to be a distance matrix and
        must be square.
    random_state: numpy.RandomState, optional
        The generator used to initialize the centers. Defaults to numpy.random.
    verbose: boolean, optional
        The verbosity level

    Returns
    -------
    core_samples: array [n_core_samples]
        Indices of core samples.

    labels : array [n_samples]
        Cluster labels for each point.  Noisy samples are given the label -1.

    Notes
    -----
    See examples/plot_dbscan.py for an example.

    **References**:
    Ester, M., H. P. Kriegel, J. Sander, and X. Xu, “A Density-Based
    Algorithm for Discovering Clusters in Large Spatial Databases with Noise”.
    In: Proceedings of the 2nd International Conference on Knowledge Discovery
    and Data Mining, Portland, OR, AAAI Press, pp. 226–231. 1996
    """
    n = X.shape[0]
    # If index order not given, create random order.
    random_state = check_random_state(random_state)
    index_order = np.arange(n)
    random_state.shuffle(index_order)
    assert len(index_order) == n, ("Index order must be of length n"
                                   " (%d expected, %d given)"
                                   % (n, len(index_order)))
    D = pairwise_distances(X, metric=metric)
    # Calculate neighborhood for all samples. This leaves the original point
    # in, which needs to be considered later (i.e. point i is the
    # neighborhood of point i. While True, its useless information)
    neighborhoods = [np.where(x <= eps)[0] for x in D]
    # Initially, all samples are noise.
    labels = -np.ones(n)
    # A list of all core samples found.
    core_samples = []
    # label_num is the label given to the new cluster
    label_num = 0
    # Look at all samples and determine if they are core.
    # If they are then build a new cluster from them.
    for index in index_order:
        if labels[index] != -1 or len(neighborhoods[index]) < min_samples:
            # This point is already classified, or not enough for a core point.
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
                noise = np.where(labels[neighborhoods[c]] == -1)[0]
                noise = neighborhoods[c][noise]
                labels[noise] = label_num
                for neighbor in noise:
                    # check if its a core point as well
                    if len(neighborhoods[neighbor]) >= min_samples:
                        # is new core point
                        new_candidates.append(neighbor)
                        core_samples.append(neighbor)
            # Update candidates for next round of cluster expansion.
            candidates = new_candidates
        # Current cluster finished.
        # Next core point found will start a new cluster.
        label_num += 1
    return core_samples, labels


class DBSCAN(BaseEstimator):
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
    verbose : boolean, optional
        The verbosity level

    Attributes
    ----------
    `core_sample_indices_` : array, shape = [n_core_samples]
        Indices of core samples.

    `components_` : array, shape = [n_core_samples, n_features]
        Copy of each core sample found by training.

    `labels_` : array, shape = [n_samples]
        Cluster labels for each point in the dataset given to fit().
        Noisy samples are given the label -1.

    Notes
    -----
    See examples/plot_dbscan.py for an example.

    **References**:
    Ester, M., H. P. Kriegel, J. Sander, and X. Xu, “A Density-Based
    Algorithm for Discovering Clusters in Large Spatial Databases with Noise”.
    In: Proceedings of the 2nd International Conference on Knowledge Discovery
    and Data Mining, Portland, OR, AAAI Press, pp. 226–231. 1996
    """

    def __init__(self, eps=0.5, min_samples=5, metric='euclidean',
                 verbose=False, random_state=None):
        self.eps = eps
        self.min_samples = min_samples
        self.metric = metric
        self.verbose = verbose
        self.random_state = check_random_state(random_state)
        self.verbose = verbose

    def fit(self, X, **params):
        """Perform DBSCAN clustering from vector array or distance matrix.

        Parameters
        ----------
        X: array [n_samples, n_samples] or [n_samples, n_features]
            Array of distances between samples, or a feature array.
            The array is treated as a feature array unless the metric is
            given as 'precomputed'.
        params: dict
            Overwrite keywords from __init__.
        """

        self.set_params(**params)
        self.core_sample_indices_, self.labels_ = dbscan(X,
                                                         **self._get_params())
        self.components_ = X[self.core_sample_indices_].copy()
        return self
