"""Unsupervised evaluation metrics."""

# Authors: Robert Layton <robertlayton@gmail.com>
#          Arnaud Fouchet <foucheta@gmail.com>
#          Thierry Guillemot <thierry.guillemot.work@gmail.com>
# License: BSD 3 clause


import functools

import numpy as np

from ...utils import check_random_state
from ...utils import check_X_y
from ...utils import safe_indexing
from ..pairwise import pairwise_distances_chunked
from ..pairwise import pairwise_distances
from ...preprocessing import LabelEncoder
from sklearn.utils import deprecated


def check_number_of_labels(n_labels, n_samples):
    """Check that number of labels are valid.

    Parameters
    ----------
    n_labels : int
        Number of labels

    n_samples : int
        Number of samples
    """
    if not 1 < n_labels < n_samples:
        raise ValueError("Number of labels is %d. Valid values are 2 "
                         "to n_samples - 1 (inclusive)" % n_labels)


def silhouette_score(X, labels, metric='euclidean', sample_size=None,
                     random_state=None, **kwds):
    """Compute the mean Silhouette Coefficient of all samples.

    The Silhouette Coefficient is calculated using the mean intra-cluster
    distance (``a``) and the mean nearest-cluster distance (``b``) for each
    sample.  The Silhouette Coefficient for a sample is ``(b - a) / max(a,
    b)``.  To clarify, ``b`` is the distance between a sample and the nearest
    cluster that the sample is not a part of.
    Note that Silhouette Coefficient is only defined if number of labels
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
        The size of the sample to use when computing the Silhouette Coefficient
        on a random subset of the data.
        If ``sample_size is None``, no sampling is used.

    random_state : int, RandomState instance or None, optional (default=None)
        The generator used to randomly select a subset of samples.  If int,
        random_state is the seed used by the random number generator; If
        RandomState instance, random_state is the random number generator; If
        None, the random number generator is the RandomState instance used by
        `np.random`. Used when ``sample_size is not None``.

    **kwds : optional keyword parameters
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
       <https://www.sciencedirect.com/science/article/pii/0377042787901257>`_

    .. [2] `Wikipedia entry on the Silhouette Coefficient
           <https://en.wikipedia.org/wiki/Silhouette_(clustering)>`_

    """
    if sample_size is not None:
        X, labels = check_X_y(X, labels, accept_sparse=['csc', 'csr'])
        random_state = check_random_state(random_state)
        indices = random_state.permutation(X.shape[0])[:sample_size]
        if metric == "precomputed":
            X, labels = X[indices].T[indices].T, labels[indices]
        else:
            X, labels = X[indices], labels[indices]
    return np.mean(silhouette_samples(X, labels, metric=metric, **kwds))


def _silhouette_reduce(D_chunk, start, labels, label_freqs):
    """Accumulate silhouette statistics for vertical chunk of X

    Parameters
    ----------
    D_chunk : shape (n_chunk_samples, n_samples)
        precomputed distances for a chunk
    start : int
        first index in chunk
    labels : array, shape (n_samples,)
        corresponding cluster labels, encoded as {0, ..., n_clusters-1}
    label_freqs : array
        distribution of cluster labels in ``labels``
    """
    # accumulate distances from each sample to each cluster
    clust_dists = np.zeros((len(D_chunk), len(label_freqs)),
                           dtype=D_chunk.dtype)
    for i in range(len(D_chunk)):
        clust_dists[i] += np.bincount(labels, weights=D_chunk[i],
                                      minlength=len(label_freqs))

    # intra_index selects intra-cluster distances within clust_dists
    intra_index = (np.arange(len(D_chunk)), labels[start:start + len(D_chunk)])
    # intra_clust_dists are averaged over cluster size outside this function
    intra_clust_dists = clust_dists[intra_index]
    # of the remaining distances we normalise and extract the minimum
    clust_dists[intra_index] = np.inf
    clust_dists /= label_freqs
    inter_clust_dists = clust_dists.min(axis=1)
    return intra_clust_dists, inter_clust_dists


def silhouette_samples(X, labels, metric='euclidean', **kwds):
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
    Note that Silhouette Coefficient is only defined if number of labels
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
       <https://www.sciencedirect.com/science/article/pii/0377042787901257>`_

    .. [2] `Wikipedia entry on the Silhouette Coefficient
       <https://en.wikipedia.org/wiki/Silhouette_(clustering)>`_

    """
    X, labels = check_X_y(X, labels, accept_sparse=['csc', 'csr'])
    le = LabelEncoder()
    labels = le.fit_transform(labels)
    n_samples = len(labels)
    label_freqs = np.bincount(labels)
    check_number_of_labels(len(le.classes_), n_samples)

    kwds['metric'] = metric
    reduce_func = functools.partial(_silhouette_reduce,
                                    labels=labels, label_freqs=label_freqs)
    results = zip(*pairwise_distances_chunked(X, reduce_func=reduce_func,
                                              **kwds))
    intra_clust_dists, inter_clust_dists = results
    intra_clust_dists = np.concatenate(intra_clust_dists)
    inter_clust_dists = np.concatenate(inter_clust_dists)

    denom = (label_freqs - 1).take(labels, mode='clip')
    with np.errstate(divide="ignore", invalid="ignore"):
        intra_clust_dists /= denom

    sil_samples = inter_clust_dists - intra_clust_dists
    with np.errstate(divide="ignore", invalid="ignore"):
        sil_samples /= np.maximum(intra_clust_dists, inter_clust_dists)
    # nan values are for clusters of size 1, and should be 0
    return np.nan_to_num(sil_samples)


def calinski_harabasz_score(X, labels):
    """Compute the Calinski and Harabasz score.

    It is also known as the Variance Ratio Criterion.

    The score is defined as ratio between the within-cluster dispersion and
    the between-cluster dispersion.

    Read more in the :ref:`User Guide <calinski_harabasz_index>`.

    Parameters
    ----------
    X : array-like, shape (``n_samples``, ``n_features``)
        List of ``n_features``-dimensional data points. Each row corresponds
        to a single data point.

    labels : array-like, shape (``n_samples``,)
        Predicted labels for each sample.

    Returns
    -------
    score : float
        The resulting Calinski-Harabasz score.

    References
    ----------
    .. [1] `T. Calinski and J. Harabasz, 1974. "A dendrite method for cluster
       analysis". Communications in Statistics
       <https://www.tandfonline.com/doi/abs/10.1080/03610927408827101>`_
    """
    X, labels = check_X_y(X, labels)
    le = LabelEncoder()
    labels = le.fit_transform(labels)

    n_samples, _ = X.shape
    n_labels = len(le.classes_)

    check_number_of_labels(n_labels, n_samples)

    extra_disp, intra_disp = 0., 0.
    mean = np.mean(X, axis=0)
    for k in range(n_labels):
        cluster_k = X[labels == k]
        mean_k = np.mean(cluster_k, axis=0)
        extra_disp += len(cluster_k) * np.sum((mean_k - mean) ** 2)
        intra_disp += np.sum((cluster_k - mean_k) ** 2)

    return (1. if intra_disp == 0. else
            extra_disp * (n_samples - n_labels) /
            (intra_disp * (n_labels - 1.)))


@deprecated("Function 'calinski_harabaz_score' has been renamed to "
            "'calinski_harabasz_score' "
            "and will be removed in version 0.23.")
def calinski_harabaz_score(X, labels):
    return calinski_harabasz_score(X, labels)


def davies_bouldin_score(X, labels):
    """Computes the Davies-Bouldin score.

    The score is defined as the ratio of within-cluster distances to
    between-cluster distances.

    Read more in the :ref:`User Guide <davies-bouldin_index>`.

    Parameters
    ----------
    X : array-like, shape (``n_samples``, ``n_features``)
        List of ``n_features``-dimensional data points. Each row corresponds
        to a single data point.

    labels : array-like, shape (``n_samples``,)
        Predicted labels for each sample.

    Returns
    -------
    score: float
        The resulting Davies-Bouldin score.

    References
    ----------
    .. [1] Davies, David L.; Bouldin, Donald W. (1979).
       `"A Cluster Separation Measure"
       <https://ieeexplore.ieee.org/document/4766909>`__.
       IEEE Transactions on Pattern Analysis and Machine Intelligence.
       PAMI-1 (2): 224-227
    """
    X, labels = check_X_y(X, labels)
    le = LabelEncoder()
    labels = le.fit_transform(labels)
    n_samples, _ = X.shape
    n_labels = len(le.classes_)
    check_number_of_labels(n_labels, n_samples)

    intra_dists = np.zeros(n_labels)
    centroids = np.zeros((n_labels, len(X[0])), dtype=np.float)
    for k in range(n_labels):
        cluster_k = safe_indexing(X, labels == k)
        centroid = cluster_k.mean(axis=0)
        centroids[k] = centroid
        intra_dists[k] = np.average(pairwise_distances(
            cluster_k, [centroid]))

    centroid_distances = pairwise_distances(centroids)

    if np.allclose(intra_dists, 0) or np.allclose(centroid_distances, 0):
        return 0.0

    score = (intra_dists[:, None] + intra_dists) / centroid_distances
    score[score == np.inf] = np.nan
    return np.mean(np.nanmax(score, axis=1))


def score_function(X, labels):
    """Computes Score Function.

    Score Function is a bounded Cluster Validity Index. It is based on
    standard cluster properties and is proved to be better than or atleast as
    good as existing validation indices for hyperspheroidal clusters. The
    score function is a combination of two terms: the distance between
    clusters (intercluster distance) and the distance inside a cluster (intra-
    cluster distance). Also, the score function can used for two purposes:
    1) to estimate the number of clusters; and
    2) to evaluate the quality of the clustering results.

    DRAWBACKS OF EXISTING VALIDITY INDICES :
    Most current validity indices only cover a subset of important aspects of
    clusters. Moreover, these indices are relevant only for data sets
    containing at least two clusters. As a result, these indices are useful in
    certain situations but are not of general-purpose.

    The score function is tested against four existing validity indices, that
    are Dunn index, Davies-Bouldin index, Silhouette index, Maulik-
    Bandyopadhyay index on several artificial and real-life data sets to
    evaluate the performance of the score function. It is found to be always
    as good or better than these indices in the case of hyperspheroidal
    clusters. It also works well on multidimensional data sets and is able to
    accommodate unique and sub-cluster cases.
    * Score Function has linear computational complexity.

    # DOCTEST FOR SCORE FUNCTION

    >>> from sklearn.cluster import KMeans
    >>> from sklearn.datasets import load_iris, load_wine
    >>>
    >>> # TEST ON IRIS DATA
    >>>
    >>> loader = load_iris()
    >>> data = loader.data
    >>>
    >>> n_samples, n_features = data.shape
    >>> unique_labels = len(np.unique(loader.target))
    >>> labels = loader.target
    >>> est = KMeans(init='random', n_clusters=unique_labels, n_init=10)
    >>> fit_data = est.fit(data)
    >>> score_function(data, est.labels_)
    0.521
    >>>
    >>> # TEST ON WINE DATA
    >>>
    >>> loader = load_wine()
    >>> data = loader.data
    >>>
    >>> n_samples, n_features = data.shape
    >>> unique_labels = len(np.unique(loader.target))
    >>> labels = loader.target
    >>> est = KMeans(init='random', n_clusters=unique_labels, n_init=10)
    >>> fit_data = est.fit(data)
    >>> score_function(data, est.labels_)
    0.161

    Parameters
    ----------
    X : array-like, shape (``n_samples``, ``n_features``)
        List of ``n_features``-dimensional data points. Each row corresponds
        to a single data point.

    labels : array-like, shape (``n_samples``,)
        Predicted labels for each sample.

    Returns
    -------
    sf: float
        The resulting Score Function.

    References
    ----------
    .. [1] Saitta S., Raphael B., Smith I.F.C. (2007).
        `"A Bounded Index for Cluster Validity"
        <https://dl.acm.org/citation.cfm?id=1420344>`__.
        Perner P.(eds) Machine Learning and Data Mining in Pattern Recognition.
        Lecture Notes in Computer Science, vol 4571.
    """
    # Pre-processing and validation of input data and labels
    X, labels = check_X_y(X, labels)
    le = LabelEncoder()
    labels = le.fit_transform(labels)

    n_samples, _ = X.shape
    n_labels = len(le.classes_)

    check_number_of_labels(n_labels, n_samples)

    # Mean of the entire data. Same as z_tot in the original paper.
    data_mean = np.mean(X, axis=0)

    # cluster_size : is a 1D-array of length "n_labels"
    # (i.e. number of unique labels or number of unique clusters)
    # It stores the number of data points in each cluster.
    cluster_size = np.zeros(n_labels)
    # intra_dists: It stores the mean euclidean distance between
    # all points belonging to a cluster and their cluster centroid.
    # Since there are "n_labels" unique clusters, the variable is
    # a 1D-array of length "n_labels".
    intra_dists = np.zeros(n_labels)
    # intra_centroid_dists: It stores the product of euclidean distance
    # between a cluster centroid and the centroid of
    # the entire data to the number of elements of that cluster.
    # Acc. to paper : intra_centroid_dists[i] = ||z_i - z_tot|| * n_i,
    # where "i" is a particular cluster.
    # Since there are "n_labels" unique clusters, the variable is a
    # 1D-array of length "n_labels".
    intra_centroid_dists = np.zeros(n_labels)
    # centroids: stores the centroid of each cluster.
    centroids = np.zeros((n_labels, len(X[0])), dtype=np.float)

    for k in range(n_labels):
        # Retrieve all data points belonging to cluster "k".
        cluster_k = safe_indexing(X, labels == k)
        # Find number of data points in cluster "k" and save them.
        cluster_size[k] = len(cluster_k)
        # Finding the centroid of cluster_k
        centroids[k] = cluster_k.mean(axis=0)

        # Compute intra_centroid_dists[k] acc. to the formula in the paper :
        # intra_centroid_dists[k] = ||z_k - z_tot|| * n_k,
        # where "k" is a particular cluster.
        intra_centroid_dists[k] = \
            pairwise_distances([centroids[k]], [data_mean],
                               metric='euclidean') * cluster_size[k]
        # Compute mean euclidean distance between all points belonging to
        # cluster "k" and the centroid of cluster "k".
        intra_dists[k] = np.mean(pairwise_distances(cluster_k, [centroids[k]],
                                                    metric='euclidean'))

    # Compute "between class distance"(bcd).
    bcd = np.mean(intra_centroid_dists) / n_samples
    # Compute "within class distance"(wcd).
    wcd = np.sum(intra_dists)

    # return score function
    return 1 - 1 / np.exp(np.exp(bcd - wcd))
