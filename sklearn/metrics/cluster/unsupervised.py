""" Unsupervised evaluation metrics. """

# Authors: Robert Layton <robertlayton@gmail.com>
#          Joel Nothman <joel.nothman@gmail.com>
#
# License: BSD 3 clause

import numpy as np

from ...utils import check_random_state, check_array, check_X_y
from ..pairwise import pairwise_distances


def silhouette_score(X, labels, metric='euclidean', sample_size=None,
                     sample_weight=None, random_state=None, **kwds):
    """Compute the mean Silhouette Coefficient of all samples.

    The Silhouette Coefficient is calculated using the mean intra-cluster
    distance (``a``) and the mean nearest-cluster distance (``b``) for each
    sample.  The Silhouette Coefficient for a sample is ``(b - a) / max(a,
    b)``.  To clarify, ``b`` is the distance between a sample and the nearest
    cluster that the sample is not a part of.
    Note that Silhouette Coefficent is only defined if number of labels
    is 2 <= n_labels <= n_samples - 1.

    This function returns the mean Silhouette Coefficient over all samples.
    To obtain the values for each sample, use :func:`silhouette_samples`.

    The best value is 1 and the worst value is -1. Values near 0 indicate
    overlapping clusters. Negative values generally indicate that a sample has
    been assigned to the wrong cluster, as a different cluster is more similar.

    Read more in the :ref:`User Guide <silhouette_coefficient>`.

    Parameters
    ----------
    X : array [n_samples_a, n_samples_a] if metric == "precomputed", or, \
             [n_samples_a, n_features] otherwise
        Array of pairwise distances between samples, or a feature array.

    labels : array, shape = [n_samples]
         Predicted labels for each sample.

    metric : string, or callable
        The metric to use when calculating distance between instances in a
        feature array. If metric is a string, it must be one of the options
        allowed by :func:`metrics.pairwise.pairwise_distances
        <sklearn.metrics.pairwise.pairwise_distances>`. If X is the distance
        array itself, use ``metric="precomputed"``.

    sample_size : int or None
        The size of the sample to use when computing the Silhouette
        Coefficient. If ``sample_size is None``, no sampling is used.
        Note this does not take ``sample_weight`` into account.

    sample_weight : array, shape = [n_samples], optional
        The weight for each sample, where samples would usually be weighted 1.
        Thus intra-cluster distance for point i will incorporate a
        zero-distance with weight ``sample_weight[i] - 1``.

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

    References
    ----------

    .. [1] `Peter J. Rousseeuw (1987). "Silhouettes: a Graphical Aid to the
       Interpretation and Validation of Cluster Analysis". Computational
       and Applied Mathematics 20: 53-65.
       <http://www.sciencedirect.com/science/article/pii/0377042787901257>`_

    .. [2] `Wikipedia entry on the Silhouette Coefficient
           <http://en.wikipedia.org/wiki/Silhouette_(clustering)>`_

    """
    n_labels = len(np.unique(labels))
    n_samples = X.shape[0]
    if not 1 < n_labels < n_samples:
        raise ValueError("Number of labels is %d. Valid values are 2 "
                         "to n_samples - 1 (inclusive)" % n_labels)

    if sample_size is not None:
        random_state = check_random_state(random_state)
        indices = random_state.permutation(X.shape[0])[:sample_size]
        if metric == "precomputed":
            X, labels = X[indices].T[indices].T, labels[indices]
        else:
            X, labels = X[indices], labels[indices]
        if sample_weight is not None:
            sample_weight = sample_weight[indices]
    return np.average(silhouette_samples(X, labels, metric=metric,
                                         sample_weight=sample_weight, **kwds),
                      weights=sample_weight)


def silhouette_samples(X, labels, metric='euclidean', sample_weight=None,
                       **kwds):
    """Compute the Silhouette Coefficient for each sample.

    The Silhouette Coefficient is a measure of how well samples are clustered
    with samples that are similar to themselves. Clustering models with a high
    Silhouette Coefficient are said to be dense, where samples in the same
    cluster are similar to each other, and well separated, where samples in
    different clusters are not very similar to each other.

    The Silhouette Coefficient is calculated using the mean intra-cluster
    distance (``a``) and the mean nearest-cluster distance (``b``) for each
    sample.  The Silhouette Coefficient for a sample is ``(b - a) / max(a,
    b)``.
    Note that Silhouette Coefficent is only defined if number of labels
    is 2 <= n_labels <= n_samples - 1.

    This function returns the Silhouette Coefficient for each sample.

    The best value is 1 and the worst value is -1. Values near 0 indicate
    overlapping clusters.

    Read more in the :ref:`User Guide <silhouette_coefficient>`.

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
        allowed by :func:`sklearn.metrics.pairwise.pairwise_distances`. If X is
        the distance array itself, use "precomputed" as the metric.

    sample_weight : array, shape = [n_samples], optional
        The weight for each sample, where samples would usually be weighted 1.
        Thus intra-cluster distance for point i will incorporate a
        zero-distance with weight ``sample_weight[i] - 1``.

    `**kwds` : optional keyword parameters
        Any further parameters are passed directly to the distance function.
        If using a ``scipy.spatial.distance`` metric, the parameters are still
        metric dependent. See the scipy docs for usage examples.

    Returns
    -------
    silhouette : array, shape = [n_samples]
        Silhouette Coefficient for each samples.

    References
    ----------

    .. [1] `Peter J. Rousseeuw (1987). "Silhouettes: a Graphical Aid to the
       Interpretation and Validation of Cluster Analysis". Computational
       and Applied Mathematics 20: 53-65.
       <http://www.sciencedirect.com/science/article/pii/0377042787901257>`_

    .. [2] `Wikipedia entry on the Silhouette Coefficient
       <http://en.wikipedia.org/wiki/Silhouette_(clustering)>`_

    """
    check_array(labels)
    distances = pairwise_distances(X, metric=metric, **kwds)
    if sample_weight is not None:
        # weight each column
        distances = distances * sample_weight
    X, labels = check_X_y(X, labels)
    n = labels.shape[0]

    # relabel to range [0, #labels - 1]:
    unique_labels, labels = np.unique(labels, return_inverse=True)
    n_labels = len(unique_labels)

    # calculate sum of distances between each point and each cluster:
    sum_distances = np.zeros((n, n_labels))
    # Shorthand available from numpy 1.8.0:
    # np.add.at(sum_distances.T, labels, distances)
    for label in range(n_labels):
        sum_distances[:, label] = np.sum(distances[:, labels == label],
                                         axis=-1)

    # Calculate distance between each point and co-clustered points
    cluster_sizes = np.bincount(labels, weights=sample_weight)
    intra_cluster = (sum_distances[np.arange(n), labels] /
                     (cluster_sizes[labels] - 1))
    intra_cluster[cluster_sizes[labels] == 1] = 0

    # Calculate minimum distance between each point and another cluster
    sum_distances[np.arange(n), labels] = np.inf
    mean_distances = sum_distances / cluster_sizes
    nearest_cluster = np.min(mean_distances, axis=1)

    sil_samples = ((nearest_cluster - intra_cluster) /
                   np.maximum(intra_cluster, nearest_cluster))
    # nan values are for clusters of size 1, and should be 0
    return np.nan_to_num(sil_samples)
