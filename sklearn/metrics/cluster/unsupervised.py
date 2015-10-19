""" Unsupervised evaluation metrics. """

# Authors: Robert Layton <robertlayton@gmail.com>
#          Alexandre Abraham <abraham.alexandre@gmail.com>
#
# License: BSD 3 clause

import warnings
from itertools import combinations

import numpy as np

from ...utils import check_random_state
from ..pairwise import pairwise_distances
from ...externals.joblib import Parallel, delayed, cpu_count


def silhouette_score(X, labels, metric='euclidean', blockwise='auto',
                     sample_size=None, n_jobs=1, random_state=None, **kwds):
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

    blockwise: {'auto', True, False}
        Enables blockwise computation of the distance matrix.
        If false, the full distance matrix is computed yielding in fast
        computation but high memory consumption. If True, it computes
        clusterwise distance matrices, dividing memory consumption by
        approximately the squared number of clusters. The latter allows
        parallelization through ``n_jobs`` parameter.
        Default is 'auto' that chooses the fastest option without
        consideration for memory consumption.

    sample_size : int or None
        The size of the sample to use when computing the Silhouette Coefficient
        on a random subset of the data.
        If ``sample_size is None``, no sampling is used.

    n_jobs : integer, optional
        The number of CPUs to use to do the computation. -1 means
        'all CPUs'. This option must be used with ``blockwise={True, 'auto'}``.
        Memory consumption is proportional to the number of CPUs.

    random_state : integer or numpy.RandomState, optional
        The generator used to randomly select a subset of samples if
        ``sample_size is not None``. If an integer is given, it fixes the seed.
        Defaults to the global numpy random number generator.

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

    if blockwise not in [True, False, 'auto']:
        raise ValueError("Blockwise parameter must be True, False or 'auto'. "
                         "You have set it to %s." % str(blockwise))

    if sample_size is not None:
        random_state = check_random_state(random_state)
        indices = random_state.permutation(X.shape[0])[:sample_size]
        if metric == "precomputed":
            #if blockwise is True:
            #    raise ValueError('Precomputed matrix is not compatible with'
            #            ' blockwise computation')
            #blockwise = False
            #n_jobs = 1
            X, labels = X[indices].T[indices].T, labels[indices]
        else:
            X, labels = X[indices], labels[indices]
    return np.mean(silhouette_samples(X, labels, metric=metric,
                                      blockwise=blockwise,
                                      n_jobs=n_jobs, **kwds))


def silhouette_samples(X, labels, metric='euclidean', blockwise='auto',
                       n_jobs=1, **kwds):
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

    blockwise: {'auto', True, False}
        Enables blockwise computation of the distance matrix.
        If false, the full distance matrix is computed yielding in fast
        computation but high memory consumption. If True, it computes
        clusterwise distance matrices, dividing memory consumption by
        approximately the squared number of clusters. The latter allows
        parallelization through ``n_jobs`` parameter.
        Default is 'auto' that chooses the fastest option without
        consideration for memory consumption.

    n_jobs : integer, optional
        The number of CPUs to use to do the computation. -1 means
        'all CPUs'. This option must be used with ``blockwise={True, 'auto'}``.
        Memory consumption is proportional to the number of CPUs.

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
    if blockwise not in [True, False, 'auto']:
        raise ValueError("Blockwise parameter must be True, False or 'auto'. "
                         "You have set it to %s." % str(blockwise))
    if blockwise is False and n_jobs != 1:
        warnings.warn('Parallelization is only available for blockwise method')
        n_jobs = 1

    if blockwise == 'auto':
        #n_labels = len(np.unique(labels))
        #if n_jobs is None:
        #    n_jobs = cpu_count()
        #blockwise = not (n_labels / n_jobs > 50 and X.shape[0] < 5 * n_labels)
        blackwise = False
    if not blockwise:
        distances = pairwise_distances(X, metric=metric, **kwds)
        n = labels.shape[0]
        A = np.array([_intra_cluster_distance(distances[i], labels, i)
                      for i in range(n)])
        B = np.array([_nearest_cluster_distance(distances[i], labels, i)
                      for i in range(n)])
    else:
        # Intra distance
        A = np.zeros(labels.size, dtype=float)
        intra_dist = Parallel(n_jobs=n_jobs)(
            delayed(_intra_cluster_distances_block)(
                X[np.where(labels == label)[0]], metric, **kwds)
            for label in np.unique(labels))
        for label, dist in zip(np.unique(labels), intra_dist):
            A[np.where(labels == label)[0]] = dist

        # Nearest cluster distance
        B = np.empty(labels.size, dtype=float)
        B.fill(np.inf)

        # Compute cluster distance between pairs of clusters
        unique_labels = np.unique(labels)
        values = Parallel(n_jobs=n_jobs)(
            delayed(_inter_cluster_distance_block)(
                    X[np.where(labels == label_a)[0]],
                    X[np.where(labels == label_b)[0]],
                    metric, **kwds)
            for label_a, label_b in combinations(unique_labels, 2)
        )

        # Take the distance to the closest cluster
        for (label_a, label_b), (values_a, values_b) in \
                    zip(combinations(unique_labels, 2), values):
                indices_a = np.where(labels == label_a)[0]
                B[indices_a] = np.minimum(values_a, B[indices_a])
                del indices_a
                indices_b = np.where(labels == label_b)[0]
                B[indices_b] = np.minimum(values_b, B[indices_b])
                del indices_b
    sil_samples = (B - A) / np.maximum(A, B)
    return sil_samples


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
    if not np.any(mask):
        # cluster of size 1
        return 0
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


def _intra_cluster_distances_block(X_cluster, metric, **kwds):
    ''' Calculate the mean intra-cluster distance for a given cluster

    Parameters
    ----------
    X_cluster : array, shape = [n_samples, n_features]
        Feature array of given cluster samples

    metric : string, or callable
        The metric to use when calculating distance between instances in a
        feature array. If metric is a string, it must be one of the options
        allowed by :func:`sklearn.metrics.pairwise.pairwise_distances`. If X is
        the distance array itself, use "precomputed" as the metric.

    `**kwds` : optional keyword parameters
        Any further parameters are passed directly to the distance function.
        If using a ``scipy.spatial.distance`` metric, the parameters are still
        metric dependent. See the scipy docs for usage examples.
    '''
    distances = pairwise_distances(X_cluster, metric=metric, **kwds)
    return distances.sum(axis=1) / (distances.shape[0] - 1)


def _inter_cluster_distance_block(X_cluster_a, X_cluster_b, metric, **kwds):
    """Calculate the mean inter-cluster distance between two clusters.

    Parameters
    ----------
    X_cluster_a : array, shape = [n_samples_a, n_features]
        Feature array of first cluster

    X_cluster_b : array, shape = [n_samples_b, n_features]
        Feature array of second cluster

    metric : string, or callable
        The metric to use when calculating distance between instances in a
        feature array. If metric is a string, it must be one of the options
        allowed by :func:`sklearn.metrics.pairwise.pairwise_distances`. If X is
        the distance array itself, use "precomputed" as the metric.

    `**kwds` : optional keyword parameters
        Any further parameters are passed directly to the distance function.
        If using a ``scipy.spatial.distance`` metric, the parameters are still
        metric dependent. See the scipy docs for usage examples.

    Returns
    -------
    (dist_a, dist_b) : tuple of array, shape = [n_samples_a] and [n_samples_b]
        Mean inter-cluster distance betweens samples of two given clusters
    """
    dist = pairwise_distances(X_cluster_a, X_cluster_b, metric=metric, **kwds)
    dist_a = dist.mean(axis=1)
    dist_b = dist.mean(axis=0)
    return dist_a, dist_b
