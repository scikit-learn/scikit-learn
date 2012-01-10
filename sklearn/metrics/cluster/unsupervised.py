""" Unsupervised evaluation metrics. """

# Authors: Robert Layton <robertlayton@gmail.com>
#
# License: BSD Style.

import numpy as np

from ...utils import check_random_state
from ..pairwise import pairwise_distances


def silhouette_score(X, labels, metric='euclidean',
                      sample_size=None, random_state=None, **kwds):
    """Compute the mean Silhouette Coefficient of all samples.

    The Silhouette Coefficient is calculated using the mean intra-cluster
    distance (a) and the mean nearest-cluster distance (b) for each sample.
    The Silhouette Coefficient for a sample is ``(b - a) / max(a, b)``.
    To clarrify, b is the distance between a sample and the nearest cluster
    that b is not a part of.

    This function returns the mean Silhoeutte Coefficient over all samples.
    To obtain the values for each sample, use silhouette_samples

    The best value is 1 and the worst value is -1. Values near 0 indicate
    overlapping clusters. Negative values generally indicate that a sample has
    been assigned to the wrong cluster, as a different cluster is more similar.

    Parameters
    ----------
    X : array [n_samples_a, n_samples_a] if metric == "precomputed", or, \
             [n_samples_a, n_features] otherwise
        Array of pairwise distances between samples, or a feature array.

    labels : array, shape = [n_samples]
             label values for each sample

    metric : string, or callable
        The metric to use when calculating distance between instances in a
        feature array. If metric is a string, it must be one of the options
        allowed by metrics.pairwise.pairwise_distances. If X is the distance
        array itself, use "precomputed" as the metric.

    sample_size : int or None
        The size of the sample to use when computing the Silhouette
        Coefficient. If sample_size is None, no sampling is used.

    random_state : integer or numpy.RandomState, optional
        The generator used to initialize the centers. If an integer is
        given, it fixes the seed. Defaults to the global numpy random
        number generator.

    `**kwds` : optional keyword parameters
        Any further parameters are passed directly to the distance function.
        If using a scipy.spatial.distance metric, the parameters are still
        metric dependent. See the scipy docs for usage examples.

    Returns
    -------
    silhouette : float
        Mean Silhouette Coefficient for all samples.

    Notes
    -----
    **References**:

    Peter J. Rousseeuw (1987). "Silhouettes: a Graphical Aid to the
        Interpretation and Validation of Cluster Analysis". Computational
        and Applied Mathematics 20: 53-65. doi:10.1016/0377-0427(87)90125-7.

    http://en.wikipedia.org/wiki/Silhouette_(clustering)

    """
    if sample_size is not None:
        random_state = check_random_state(random_state)
        indices = random_state.permutation(X.shape[0])[:sample_size]
        if metric == "precomputed":
            X, labels = X[indices].T[indices].T, labels[indices]
        else:
            X, labels = X[indices], labels[indices]
    return np.mean(silhouette_samples(X, labels, metric=metric, **kwds))


def silhouette_samples(X, labels, metric='euclidean', **kwds):
    """Compute the Silhouette Coefficient for each sample.

    The Silhoeutte Coefficient is a measure of how well samples are clustered
    with samples that are similar to themselves. Clustering models with a high
    Silhouette Coefficient are said to be dense, where samples in the same
    cluster are similar to each other, and well separated, where samples in
    different clusters are not very similar to each other.

    The Silhouette Coefficient is calculated using the mean intra-cluster
    distance (a) and the mean nearest-cluster distance (b) for each sample.
    The Silhouette Coefficient for a sample is ``(b - a) / max(a, b)``.

    This function returns the Silhoeutte Coefficient for each sample.

    The best value is 1 and the worst value is -1. Values near 0 indicate
    overlapping clusters.

    Parameters
    ----------
    X : array [n_samples_a, n_samples_a] if metric == "precomputed", or, \
             [n_samples_a, n_features] otherwise
        Array of pairwise distances between samples, or a feature array.

    labels : array, shape = [n_samples]
             label values for each sample

    metric : string, or callable
        The metric to use when calculating distance between instances in a
        feature array. If metric is a string, it must be one of the options
        allowed by metrics.pairwise.pairwise_distances. If X is the distance
        array itself, use "precomputed" as the metric.

    `**kwds` : optional keyword parameters
        Any further parameters are passed directly to the distance function.
        If using a scipy.spatial.distance metric, the parameters are still
        metric dependent. See the scipy docs for usage examples.

    Returns
    -------
    silhouette : array, shape = [n_samples]
        Silhouette Coefficient for each samples.

    Notes
    -----
    **References**:

    Peter J. Rousseeuw (1987). "Silhouettes: a Graphical Aid to the
        Interpretation and Validation of Cluster Analysis". Computational
        and Applied Mathematics 20: 53-65. doi:10.1016/0377-0427(87)90125-7.

    http://en.wikipedia.org/wiki/Silhouette_(clustering)

    """
    distances = pairwise_distances(X, metric=metric, **kwds)
    n = labels.shape[0]
    A = np.array([_intra_cluster_distance(distances[i], labels, i)
                  for i in range(n)])
    B = np.array([_nearest_cluster_distance(distances[i], labels, i)
                  for i in range(n)])
    return (B - A) / np.maximum(A, B)


def _intra_cluster_distance(distances_row, labels, i):
    """Calculate the mean intra-cluster distance for sample i.

    Parameters
    ----------
    distances_row : array, shape = [n_samples]
        Pairwise distance matrix between sample i and each sample.

    labels : array, shape = [n_samples]
        label values for each sample

    i : int
        Sample index being calculated. It is excluded from calculation and
        used to determine the current label

    Returns
    -------
    a : float
        Mean intra-cluster distance for sample i
    """
    mask = labels == labels[i]
    mask[i] = False
    a = np.mean(distances_row[mask])
    return a


def _nearest_cluster_distance(distances_row, labels, i):
    """Calculate the mean nearest-cluster distance for sample i.

    Parameters
    ----------
    distances_row : array, shape = [n_samples]
        Pairwise distance matrix between sample i and each sample.

    labels : array, shape = [n_samples]
        label values for each sample

    i : int
        Sample index being calculated. It is used to determine the current
        label.

    Returns
    -------
    b : float
        Mean nearest-cluster distance for sample i
    """
    label = labels[i]
    b = np.min([np.mean(distances_row[labels == cur_label])
               for cur_label in set(labels) if not cur_label == label])
    return b
