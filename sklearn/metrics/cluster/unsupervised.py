""" Unsupervised evaluation metrics


"""

# Authors: Robert Layton <robertlayton@gmail.com>
#
# License: BSD Style.

import numpy as np


def silhouette_score(distances, labels):
    """Compute the mean Silhouette Coefficient of all samples.

    The Silhouette Coefficient is calculated using the mean intra-cluster
    distance (a) and the mean nearest-cluster distance (b) for each sample.
    The Silhouette Coefficient for a sample is (b - a) / max(a, b).
    To clarrify, b is the distance between a sample and the nearest cluster
    that b is not a part of.

    This function returns the mean Silhoeutte Coefficient over all samples.
    To obtain the values for each sample, use silhouette_samples

    The best value is 1 and the worst value is -1. Values near 0 indicate
    overlapping clusters. Negative values generally indicate that a sample has
    been assigned to the wrong cluster, as a different cluster is more similar.

    Parameters
    ----------
    distances : array, shape = [n_samples, n_samples]
                Pairwise distance matrix between each sample.

    labels : array, shape = [n_samples]
             label values for each sample

    Returns
    -------
    silhouette : float
        Mean Silhouette Coefficient for all samples.

    References
    ----------
    Peter J. Rousseeuw (1987). "Silhouettes: a Graphical Aid to the
        Interpretation and Validation of Cluster Analysis". Computational
        and Applied Mathematics 20: 53–65. doi:10.1016/0377-0427(87)90125-7.

    http://en.wikipedia.org/wiki/Silhouette_(clustering)

    """
    return np.mean(silhouette_samples(distances, labels))


def silhouette_samples(distances, labels):
    """Compute the Silhouette Coefficient for each sample.

    The Silhoeutte Coefficient is a measure of how well samples are clustered
    with samples that are similar to themselves. Clustering models with a high
    Silhouette Coefficient are said to be dense, where samples in the same
    cluster are similar to each other, and well separated, where samples in
    different clusters are not very similar to each other.

    The Silhouette Coefficient is calculated using the mean intra-cluster
    distance (a) and the mean nearest-cluster distance (b) for each sample.
    The Silhouette Coefficient for a sample is (b - a) / max(a, b).

    This function returns the Silhoeutte Coefficient for each sample.

    The best value is 1 and the worst value is -1. Values near 0 indicate
    overlapping clusters.

    Parameters
    ----------
    distances : array, shape = [n_samples, n_samples]
                Pairwise distance matrix between each sample.

    labels : array, shape = [n_samples]
             label values for each sample

    Returns
    -------
    silhouette : array, shape = [n_samples]
        Silhouette Coefficient for each samples.

    References
    ----------
    Peter J. Rousseeuw (1987). "Silhouettes: a Graphical Aid to the
        Interpretation and Validation of Cluster Analysis". Computational
        and Applied Mathematics 20: 53–65. doi:10.1016/0377-0427(87)90125-7.

    http://en.wikipedia.org/wiki/Silhouette_(clustering)

    """
    n = labels.shape[0]
    A = np.array([_intra_cluster_distance(distances[i], labels, i)
                  for i in range(n)])
    B = np.array([_nearest_cluster_distance(distances[i], labels, i)
                  for i in range(n)])
    return (B - A) / np.maximum(A, B)


def _intra_cluster_distance(distances_row, labels, i):
    """ Calculate the mean intra-cluster distance for sample i.

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
    label = labels[i]
    a = np.mean([distances_row[j] for j in range(len(distances_row))
                 if labels[j] == label and not i == j])
    return a


def _nearest_cluster_distance(distances_row, labels, i):
    """ Calculate the mean nearest-cluster distance for sample i.

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
    b = np.min([np.mean([distances_row[j] for j in range(len(distances_row))
                         if labels[j] == cur_label])
               for cur_label in set(labels) if not cur_label == label])
    return b
