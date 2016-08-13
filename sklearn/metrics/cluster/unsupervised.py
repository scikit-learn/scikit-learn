"""Unsupervised evaluation metrics."""

# Authors: Robert Layton <robertlayton@gmail.com>
#          Arnaud Fouchet <foucheta@gmail.com>
#          Thierry Guillemot <thierry.guillemot.work@gmail.com>
# License: BSD 3 clause

from __future__ import division

import numpy as np

from ...utils import check_random_state
from ...utils import check_X_y
from ...utils import _get_n_jobs
from ...externals.joblib import Parallel, delayed
from ..pairwise import pairwise_distances
from ...preprocessing import LabelEncoder


def check_number_of_labels(n_labels, n_samples):
    if not 1 < n_labels < n_samples:
        raise ValueError("Number of labels is %d. Valid values are 2 "
                         "to n_samples - 1 (inclusive)" % n_labels)


DEFAULT_BLOCK_SIZE = 2 ** 26
BYTES_PER_FLOAT = 8


def silhouette_score(X, labels, metric='euclidean', sample_size=None,
                     block_size=DEFAULT_BLOCK_SIZE, n_jobs=1,
                     random_state=None, **kwds):
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

    block_size : int, optional
        The maximum number of bytes of memory per job (see ``n_jobs``) to use
        at a time for calculating pairwise distances. Default is 64MiB.

        .. versionadded:: 0.18

    n_jobs : int, optional (default = 1)
        The number of parallel jobs to run.
        If ``-1``, then the number of jobs is set to the number of CPU cores.

        .. versionadded:: 0.18

    sample_size : int or None
        The size of the sample to use when computing the Silhouette Coefficient
        on a random subset of the data.
        If ``sample_size is None``, no sampling is used.

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
           <https://en.wikipedia.org/wiki/Silhouette_(clustering)>`_

    """
    X, labels = check_X_y(X, labels, accept_sparse=['csc', 'csr'])
    le = LabelEncoder()
    labels = le.fit_transform(labels)
    n_labels = len(le.classes_)
    n_samples = X.shape[0]

    check_number_of_labels(n_labels, n_samples)

    if sample_size is not None:
        random_state = check_random_state(random_state)
        indices = random_state.permutation(X.shape[0])[:sample_size]
        if metric == "precomputed":
            X, labels = X[indices].T[indices].T, labels[indices]
        else:
            X, labels = X[indices], labels[indices]
    return np.mean(silhouette_samples(X, labels, metric=metric,
                                      block_size=block_size, **kwds))


def _process_block(X, labels, start, block_n_rows, block_range, add_at,
                   label_freqs, metric, kwds):
    # get distances from block to every other sample
    stop = min(start + block_n_rows, X.shape[0])
    if stop - start == X.shape[0]:
        # allow pairwise_distances to use fast paths
        block_dists = pairwise_distances(X, metric=metric, **kwds)
    else:
        block_dists = pairwise_distances(X[start:stop], X,
                                         metric=metric, **kwds)

    # accumulate distances from each sample to each cluster
    clust_dists = np.bincount(add_at[:block_dists.size],
                              block_dists.ravel())
    clust_dists = clust_dists.reshape((stop - start, len(label_freqs)))

    # intra_index selects intra-cluster distances within clust_dists
    intra_index = (block_range[:len(clust_dists)], labels[start:stop])
    # intra_clust_dists are averaged over cluster size outside this function
    intra_clust_dists = clust_dists[intra_index]
    # of the remaining distances we normalise and extract the minimum
    clust_dists[intra_index] = np.inf
    clust_dists /= label_freqs
    inter_clust_dists = clust_dists.min(axis=1)
    return intra_clust_dists, inter_clust_dists


def silhouette_samples(X, labels, metric='euclidean',
                       block_size=DEFAULT_BLOCK_SIZE, n_jobs=1, **kwds):
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

    block_size : int, optional
        The maximum number of bytes of memory per job (see ``n_jobs``) to use
        at a time for calculating pairwise distances. Default is 64MiB.

        .. versionadded:: 0.18

    n_jobs : int, optional (default = 1)
        The number of parallel jobs to run.
        If ``-1``, then the number of jobs is set to the number of CPU cores.

        .. versionadded:: 0.18

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
       <https://en.wikipedia.org/wiki/Silhouette_(clustering)>`_

    """
    X, labels = check_X_y(X, labels, accept_sparse=['csc', 'csr'])
    le = LabelEncoder()
    labels = le.fit_transform(labels)
    n_samples = len(labels)
    label_freqs = np.bincount(labels)

    n_jobs = _get_n_jobs(n_jobs)
    block_n_rows = block_size // (BYTES_PER_FLOAT * n_samples)
    if block_n_rows > n_samples:
        block_n_rows = min(block_n_rows, n_samples)
    if block_n_rows < 1:
        raise ValueError('block_size should be at least n_samples * %d = %d '
                         'bytes, got %r' % (BYTES_PER_FLOAT,
                                            n_samples * BYTES_PER_FLOAT,
                                            block_size))

    intra_clust_dists = []
    inter_clust_dists = []

    # We use these indices as bins to accumulate distances from each sample in
    # a block to each cluster.
    # NB: we currently use np.bincount but could use np.add.at when Numpy >=1.8
    # is minimum dependency, which would avoid materialising this index.
    block_range = np.arange(block_n_rows)
    add_at = np.ravel_multi_index((np.repeat(block_range, n_samples),
                                   np.tile(labels, block_n_rows)),
                                  dims=(block_n_rows, len(label_freqs)))
    parallel = Parallel(n_jobs=n_jobs)

    results = parallel(delayed(_process_block)(X, labels, start, block_n_rows,
                                               block_range, add_at,
                                               label_freqs, metric, kwds)
                       for start in range(0, n_samples, block_n_rows))

    intra_clust_dists, inter_clust_dists = zip(*results)
    if len(intra_clust_dists) == 1:
        intra_clust_dists = intra_clust_dists[0]
        inter_clust_dists = inter_clust_dists[0]
    else:
        intra_clust_dists = np.hstack(intra_clust_dists)
        inter_clust_dists = np.hstack(inter_clust_dists)

    denom = (label_freqs - 1).take(labels, mode='clip')
    with np.errstate(divide="ignore", invalid="ignore"):
        intra_clust_dists /= denom

    sil_samples = inter_clust_dists - intra_clust_dists
    with np.errstate(divide="ignore", invalid="ignore"):
        sil_samples /= np.maximum(intra_clust_dists, inter_clust_dists)
    # nan values are for clusters of size 1, and should be 0
    return np.nan_to_num(sil_samples)


def calinski_harabaz_score(X, labels):
    """Compute the Calinski and Harabaz score.

    The score is defined as ratio between the within-cluster dispersion and
    the between-cluster dispersion.

    Read more in the :ref:`User Guide <calinski_harabaz_index>`.

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
        The resulting Calinski-Harabaz score.

    References
    ----------
    .. [1] `T. Calinski and J. Harabasz, 1974. "A dendrite method for cluster
       analysis". Communications in Statistics
       <http://www.tandfonline.com/doi/abs/10.1080/03610927408827101>`_
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
