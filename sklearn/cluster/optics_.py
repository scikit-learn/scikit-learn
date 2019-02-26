# -*- coding: utf-8 -*-
"""Ordering Points To Identify the Clustering Structure (OPTICS)

These routines execute the OPTICS algorithm, and implement various
cluster extraction methods of the ordered list.

Authors: Shane Grigsby <refuge@rocktalus.com>
         Amy X. Zhang <axz@mit.edu>
         Erich Schubert <erich@debian.org>
License: BSD 3 clause
"""

import warnings
import numpy as np

from ..utils import check_array
from ..utils import gen_batches, get_chunk_n_rows
from ..neighbors import NearestNeighbors
from ..base import BaseEstimator, ClusterMixin
from ..metrics import pairwise_distances


def optics(X, min_samples=5, max_eps=np.inf, metric='minkowski',
           p=2, metric_params=None, extract_method='sqlnk',
           eps=0.5, maxima_ratio=.75,
           rejection_ratio=.7, similarity_threshold=0.4,
           significant_min=.003, min_cluster_size=.005,
           min_maxima_ratio=0.001, algorithm='auto',
           leaf_size=30, n_jobs=None):
    """Estimate clustering structure from vector array

    OPTICS: Ordering Points To Identify the Clustering Structure
    Closely related to DBSCAN, finds core sample of high density and expands
    clusters from them [1]_. Unlike DBSCAN, keeps cluster hierarchy for a
    variable neighborhood radius. Better suited for usage on large point
    datasets than the current sklearn implementation of DBSCAN.

    Clusters are then extracted using a DBSCAN like method [1]_, or SQLNK [2]_.

    This implementation deviates from the original OPTICS by first performing
    k-nearest-neighborhood searches on all points to identify core sizes, then
    computing only the distances to unprocessed points when constructing the
    cluster order. Note that we do not employ a heap to manage the expansion
    candidates, so the time complexity will be O(n^2).

    Read more in the :ref:`User Guide <optics>`.

    Parameters
    ----------
    X : array, shape (n_samples, n_features)
        The data.

    min_samples : int (default=5)
        The number of samples in a neighborhood for a point to be considered
        as a core point.

    max_eps : float, optional (default=np.inf)
        The maximum distance between two samples for them to be considered
        as in the same neighborhood. Default value of "np.inf" will identify
        clusters across all scales; reducing `max_eps` will result in
        shorter run times.

    metric : string or callable, optional (default='minkowski')
        metric to use for distance computation. Any metric from scikit-learn
        or scipy.spatial.distance can be used.

        If metric is a callable function, it is called on each
        pair of instances (rows) and the resulting value recorded. The callable
        should take two arrays as input and return one value indicating the
        distance between them. This works for Scipy's metrics, but is less
        efficient than passing the metric name as a string.

        Distance matrices are not supported.

        Valid values for metric are:

        - from scikit-learn: ['cityblock', 'cosine', 'euclidean', 'l1', 'l2',
          'manhattan']

        - from scipy.spatial.distance: ['braycurtis', 'canberra', 'chebyshev',
          'correlation', 'dice', 'hamming', 'jaccard', 'kulsinski',
          'mahalanobis', 'minkowski', 'rogerstanimoto', 'russellrao',
          'seuclidean', 'sokalmichener', 'sokalsneath', 'sqeuclidean',
          'yule']

        See the documentation for scipy.spatial.distance for details on these
        metrics.

    p : integer, optional (default=2)
        Parameter for the Minkowski metric from
        :class:`sklearn.metrics.pairwise_distances`. When p = 1, this is
        equivalent to using manhattan_distance (l1), and euclidean_distance
        (l2) for p = 2. For arbitrary p, minkowski_distance (l_p) is used.

    metric_params : dict, optional (default=None)
        Additional keyword arguments for the metric function.

    extract_method : string, optional (default='sqlnk')
        The extraction method used to extract clusters using the calculated
        reachability and ordering. Possible values are "dbscan"
        and "sqlnk".

    eps : float, optional (default=0.5)
        The maximum distance between two samples for them to be considered
        as in the same neighborhood. Used ony when `extract_method='dbscan'`.

    maxima_ratio : float, optional (default=.75)
        The maximum ratio we allow of average height of clusters on the
        right and left to the local maxima in question. The higher the
        ratio, the more generous the algorithm is to preserving local
        minima, and the more cuts the resulting tree will have.
        Used only when `extract_method='sqlnk'`.

    rejection_ratio : float, optional (default=.7)
        Adjusts the fitness of the clustering. When the maxima_ratio is
        exceeded, determine which of the clusters to the left and right to
        reject based on rejection_ratio. Higher values will result in points
        being more readily classified as noise; conversely, lower values will
        result in more points being clustered.
        Used only when `extract_method='sqlnk'`.

    similarity_threshold : float, optional (default=.4)
        Used to check if nodes can be moved up one level, that is, if the
        new cluster created is too "similar" to its parent, given the
        similarity threshold. Similarity can be determined by 1) the size
        of the new cluster relative to the size of the parent node or
        2) the average of the reachability values of the new cluster
        relative to the average of the reachability values of the parent
        node. A lower value for the similarity threshold means less levels
        in the tree.
        Used only when `extract_method='sqlnk'`.

    significant_min : float, optional (default=.003)
        Sets a lower threshold on how small a significant maxima can be.
        Used only when `extract_method='sqlnk'`.

    min_cluster_size : int > 1 or float between 0 and 1 (default=0.005)
        Minimum number of samples in an OPTICS cluster, expressed as an
        absolute number or a fraction of the number of samples (rounded
        to be at least 2).
        Used only when `extract_method='sqlnk'`.

    min_maxima_ratio : float, optional (default=.001)
        Used to determine neighborhood size for minimum cluster membership.
        Each local maxima should be a largest value in a neighborhood
        of the `size min_maxima_ratio * len(X)` from left and right.
        Used only when `extract_method='sqlnk'`.

    algorithm : {'auto', 'ball_tree', 'kd_tree', 'brute'}, optional
        Algorithm used to compute the nearest neighbors:

        - 'ball_tree' will use :class:`BallTree`
        - 'kd_tree' will use :class:`KDTree`
        - 'brute' will use a brute-force search.
        - 'auto' will attempt to decide the most appropriate algorithm
          based on the values passed to :meth:`fit` method. (default)

        Note: fitting on sparse input will override the setting of
        this parameter, using brute force.

    leaf_size : int, optional (default=30)
        Leaf size passed to :class:`BallTree` or :class:`KDTree`. This can
        affect the speed of the construction and query, as well as the memory
        required to store the tree. The optimal value depends on the
        nature of the problem.

    n_jobs : int or None, optional (default=None)
        The number of parallel jobs to run for neighbors search.
        ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
        ``-1`` means using all processors. See :term:`Glossary <n_jobs>`
        for more details.

    Returns
    -------
    core_sample_indices_ : array, shape (n_core_samples,)
        Indices of core samples.

    labels_ : array, shape (n_samples,)
        Cluster labels for each point in the dataset given to fit().
        Noisy samples are given the label -1.

    See also
    --------
    OPTICS
        An estimator interface for this clustering algorithm.
    dbscan
        A similar clustering for a specified neighborhood radius (eps).
        Our implementation is optimized for runtime.

    References
    ----------
    .. [1] Ankerst, Mihael, Markus M. Breunig, Hans-Peter Kriegel,
       and Jörg Sander. "OPTICS: ordering points to identify the clustering
       structure." ACM SIGMOD Record 28, no. 2 (1999): 49-60.

    .. [2] Sander, Jörg, Xuejie Qin, Zhiyong Lu, Nan Niu, and Alex Kovarsky.
       "Automatic extraction of clusters from hierarchical clustering
       representations." In Advances in Knowledge Discovery and Data
       Mining, pp. 75-87. Springer Berlin Heidelberg, 2003.

    .. [3] Schubert, Erich, Michael Gertz.
       "Improving the Cluster Structure Extracted from OPTICS Plots." Proc. of
       the Conference "Lernen, Wissen, Daten, Analysen" (LWDA) (2018): 318-329.
    """

    clust = OPTICS(min_samples, max_eps, metric, p, metric_params,
                   extract_method, eps, maxima_ratio, rejection_ratio,
                   similarity_threshold, significant_min,
                   min_cluster_size, min_maxima_ratio,
                   algorithm, leaf_size, n_jobs)
    clust.fit(X)
    return clust.core_sample_indices_, clust.labels_


class OPTICS(BaseEstimator, ClusterMixin):
    """Estimate clustering structure from vector array

    OPTICS: Ordering Points To Identify the Clustering Structure
    Closely related to DBSCAN, finds core sample of high density and expands
    clusters from them [1]_. Unlike DBSCAN, keeps cluster hierarchy for a
    variable neighborhood radius. Better suited for usage on large point
    datasets than the current sklearn implementation of DBSCAN.

    Clusters are then extracted using a DBSCAN like method [1]_, or SQLNK [2]_.

    This implementation deviates from the original OPTICS by first performing
    k-nearest-neighborhood searches on all points to identify core sizes, then
    computing only the distances to unprocessed points when constructing the
    cluster order. Note that we do not employ a heap to manage the expansion
    candidates, so the time complexity will be O(n^2).

    Read more in the :ref:`User Guide <optics>`.

    Parameters
    ----------
    min_samples : int (default=5)
        The number of samples in a neighborhood for a point to be considered
        as a core point.

    max_eps : float, optional (default=np.inf)
        The maximum distance between two samples for them to be considered
        as in the same neighborhood. Default value of "np.inf" will identify
        clusters across all scales; reducing `max_eps` will result in
        shorter run times.

    metric : string or callable, optional (default='minkowski')
        metric to use for distance computation. Any metric from scikit-learn
        or scipy.spatial.distance can be used.

        If metric is a callable function, it is called on each
        pair of instances (rows) and the resulting value recorded. The callable
        should take two arrays as input and return one value indicating the
        distance between them. This works for Scipy's metrics, but is less
        efficient than passing the metric name as a string.

        Distance matrices are not supported.

        Valid values for metric are:

        - from scikit-learn: ['cityblock', 'cosine', 'euclidean', 'l1', 'l2',
          'manhattan']

        - from scipy.spatial.distance: ['braycurtis', 'canberra', 'chebyshev',
          'correlation', 'dice', 'hamming', 'jaccard', 'kulsinski',
          'mahalanobis', 'minkowski', 'rogerstanimoto', 'russellrao',
          'seuclidean', 'sokalmichener', 'sokalsneath', 'sqeuclidean',
          'yule']

        See the documentation for scipy.spatial.distance for details on these
        metrics.

    p : integer, optional (default=2)
        Parameter for the Minkowski metric from
        :class:`sklearn.metrics.pairwise_distances`. When p = 1, this is
        equivalent to using manhattan_distance (l1), and euclidean_distance
        (l2) for p = 2. For arbitrary p, minkowski_distance (l_p) is used.

    metric_params : dict, optional (default=None)
        Additional keyword arguments for the metric function.

    extract_method : string, optional (default='sqlnk')
        The extraction method used to extract clusters using the calculated
        reachability and ordering. Possible values are "dbscan"
        and "sqlnk".

    eps : float, optional (default=0.5)
        The maximum distance between two samples for them to be considered
        as in the same neighborhood. Used ony when `extract_method='dbscan'`.

    maxima_ratio : float, optional (default=.75)
        The maximum ratio we allow of average height of clusters on the
        right and left to the local maxima in question. The higher the
        ratio, the more generous the algorithm is to preserving local
        minima, and the more cuts the resulting tree will have.
        Used only when `extract_method='sqlnk'`.

    rejection_ratio : float, optional (default=.7)
        Adjusts the fitness of the clustering. When the maxima_ratio is
        exceeded, determine which of the clusters to the left and right to
        reject based on rejection_ratio. Higher values will result in points
        being more readily classified as noise; conversely, lower values will
        result in more points being clustered.
        Used only when `extract_method='sqlnk'`.

    similarity_threshold : float, optional (default=.4)
        Used to check if nodes can be moved up one level, that is, if the
        new cluster created is too "similar" to its parent, given the
        similarity threshold. Similarity can be determined by 1) the size
        of the new cluster relative to the size of the parent node or
        2) the average of the reachability values of the new cluster
        relative to the average of the reachability values of the parent
        node. A lower value for the similarity threshold means less levels
        in the tree.
        Used only when `extract_method='sqlnk'`.

    significant_min : float, optional (default=.003)
        Sets a lower threshold on how small a significant maxima can be.
        Used only when `extract_method='sqlnk'`.

    min_cluster_size : int > 1 or float between 0 and 1 (default=0.005)
        Minimum number of samples in an OPTICS cluster, expressed as an
        absolute number or a fraction of the number of samples (rounded
        to be at least 2).
        Used only when `extract_method='sqlnk'`.

    min_maxima_ratio : float, optional (default=.001)
        Used to determine neighborhood size for minimum cluster membership.
        Each local maxima should be a largest value in a neighborhood
        of the `size min_maxima_ratio * len(X)` from left and right.
        Used only when `extract_method='sqlnk'`.

    algorithm : {'auto', 'ball_tree', 'kd_tree', 'brute'}, optional
        Algorithm used to compute the nearest neighbors:

        - 'ball_tree' will use :class:`BallTree`
        - 'kd_tree' will use :class:`KDTree`
        - 'brute' will use a brute-force search.
        - 'auto' will attempt to decide the most appropriate algorithm
          based on the values passed to :meth:`fit` method. (default)

        Note: fitting on sparse input will override the setting of
        this parameter, using brute force.

    leaf_size : int, optional (default=30)
        Leaf size passed to :class:`BallTree` or :class:`KDTree`. This can
        affect the speed of the construction and query, as well as the memory
        required to store the tree. The optimal value depends on the
        nature of the problem.

    n_jobs : int or None, optional (default=None)
        The number of parallel jobs to run for neighbors search.
        ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
        ``-1`` means using all processors. See :term:`Glossary <n_jobs>`
        for more details.

    Attributes
    ----------
    core_sample_indices_ : array, shape (n_core_samples,)
        Indices of core samples.

    labels_ : array, shape (n_samples,)
        Cluster labels for each point in the dataset given to fit().
        Noisy samples are given the label -1.

    reachability_ : array, shape (n_samples,)
        Reachability distances per sample, indexed by object order. Use
        ``clust.reachability_[clust.ordering_]`` to access in cluster order.

    ordering_ : array, shape (n_samples,)
        The cluster ordered list of sample indices.

    core_distances_ : array, shape (n_samples,)
        Distance at which each sample becomes a core point, indexed by object
        order. Points which will never be core have a distance of inf. Use
        ``clust.core_distances_[clust.ordering_]`` to access in cluster order.

    predecessor_ : array, shape (n_samples,)
        Point that a sample was reached from, indexed by object order.
        Seed points have a predecessor of -1.

    See also
    --------

    DBSCAN
        A similar clustering for a specified neighborhood radius (eps).
        Our implementation is optimized for runtime.

    References
    ----------
    .. [1] Ankerst, Mihael, Markus M. Breunig, Hans-Peter Kriegel,
       and Jörg Sander. "OPTICS: ordering points to identify the clustering
       structure." ACM SIGMOD Record 28, no. 2 (1999): 49-60.

    .. [2] Sander, Jörg, Xuejie Qin, Zhiyong Lu, Nan Niu, and Alex Kovarsky.
       "Automatic extraction of clusters from hierarchical clustering
       representations." In Advances in Knowledge Discovery and Data
       Mining, pp. 75-87. Springer Berlin Heidelberg, 2003.

    .. [3] Schubert, Erich, Michael Gertz.
       "Improving the Cluster Structure Extracted from OPTICS Plots." Proc. of
       the Conference "Lernen, Wissen, Daten, Analysen" (LWDA) (2018): 318-329.
    """

    def __init__(self, min_samples=5, max_eps=np.inf, metric='minkowski',
                 p=2, metric_params=None, extract_method='sqlnk',
                 eps=0.5, maxima_ratio=.75,
                 rejection_ratio=.7, similarity_threshold=0.4,
                 significant_min=.003, min_cluster_size=.005,
                 min_maxima_ratio=0.001, algorithm='auto',
                 leaf_size=30, n_jobs=None):

        self.max_eps = max_eps
        self.min_samples = min_samples
        self.maxima_ratio = maxima_ratio
        self.rejection_ratio = rejection_ratio
        self.similarity_threshold = similarity_threshold
        self.significant_min = significant_min
        self.min_cluster_size = min_cluster_size
        self.min_maxima_ratio = min_maxima_ratio
        self.algorithm = algorithm
        self.metric = metric
        self.metric_params = metric_params
        self.p = p
        self.leaf_size = leaf_size
        self.extract_method = extract_method
        self.eps = eps
        self.n_jobs = n_jobs

    def fit(self, X, y=None):
        """Perform OPTICS clustering

        Extracts an ordered list of points and reachability distances, and
        performs initial clustering using `max_eps` distance specified at
        OPTICS object instantiation.

        Parameters
        ----------
        X : array, shape (n_samples, n_features)
            The data.

        y : ignored

        Returns
        -------
        self : instance of OPTICS
            The instance.
        """
        X = check_array(X, dtype=np.float)

        n_samples = len(X)

        if self.min_samples > n_samples:
            raise ValueError("Number of training samples (n_samples=%d) must "
                             "be greater than min_samples (min_samples=%d) "
                             "used for clustering." %
                             (n_samples, self.min_samples))

        if self.min_cluster_size <= 0 or (self.min_cluster_size !=
                                          int(self.min_cluster_size)
                                          and self.min_cluster_size > 1):
            raise ValueError('min_cluster_size must be a positive integer or '
                             'a float between 0 and 1. Got %r' %
                             self.min_cluster_size)
        elif self.min_cluster_size > n_samples:
            raise ValueError('min_cluster_size must be no greater than the '
                             'number of samples (%d). Got %d' %
                             (n_samples, self.min_cluster_size))

        if self.extract_method not in ['dbscan', 'sqlnk']:
            raise ValueError("extract_method should be one of"
                             " 'dbscan' or 'sqlnk', but is %s" %
                             self.extract_method)

        if self.extract_method == 'dbscan':
            if self.eps > self.max_eps:
                raise ValueError('Specify an epsilon smaller than %s. Got %s.'
                                 % (self.max_eps, self.eps))

            if self.eps * 5.0 > (self.max_eps * 1.05):
                warnings.warn(
                    "Warning, max_eps (%s) is close to eps (%s): "
                    "Output may be unstable." % (self.max_eps, self.eps),
                    RuntimeWarning, stacklevel=2)
                # Stability warning is documented in _extract_dbscan method...

        (self.ordering_, self.core_distances_, self.reachability_,
         self.predecessor_) = _calculate_optics_order(
             X=X, min_samples=self.min_samples, algorithm=self.algorithm,
             leaf_size=self.leaf_size, metric=self.metric,
             metric_params=self.metric_params, p=self.p, n_jobs=self.n_jobs,
             max_eps=self.max_eps)
        # Start all points as noise ##
        self.labels_ = np.full(n_samples, -1, dtype=int)

        # Extract clusters from the calculated orders and reachability
        if self.extract_method == 'sqlnk':
            extract_params = {
                'reachability': self.reachability_,
                'ordering': self.ordering_,
                'maxima_ratio': self.maxima_ratio,
                'rejection_ratio': self.rejection_ratio,
                'similarity_threshold': self.similarity_threshold,
                'significant_min': self.significant_min,
                'min_cluster_size': self.min_cluster_size,
                'min_maxima_ratio': self.min_maxima_ratio
            }
            indices_, labels_ = _extract_optics_sqlnk(**extract_params)
        elif self.extract_method == 'dbscan':
            indices_, labels_ = _extract_optics_dbscan(self.reachability_,
                                                       self.core_distances_,
                                                       self.ordering_,
                                                       self.eps)

        self.core_sample_indices_ = indices_
        self.labels_ = labels_
        return self


# OPTICS helper functions
def _compute_core_distances_(X, neighbors, min_samples, working_memory):
    """Compute the k-th nearest neighbor of each sample

    Equivalent to neighbors.kneighbors(X, self.min_samples)[0][:, -1]
    but with more memory efficiency.

    Parameters
    ----------
    X : array, shape (n_samples, n_features)
        The data.
    neighbors : NearestNeighbors instance
        The fitted nearest neighbors estimator.
    working_memory : int, optional
        The sought maximum memory for temporary distance matrix chunks.
        When None (default), the value of
        ``sklearn.get_config()['working_memory']`` is used.

    Returns
    -------
    core_distances : array, shape (n_samples,)
        Distance at which each sample becomes a core point.
        Points which will never be core have a distance of inf.
    """
    n_samples = len(X)
    core_distances = np.empty(n_samples)
    core_distances.fill(np.nan)

    chunk_n_rows = get_chunk_n_rows(row_bytes=16 * min_samples,
                                    max_n_rows=n_samples,
                                    working_memory=working_memory)
    slices = gen_batches(n_samples, chunk_n_rows)
    for sl in slices:
        core_distances[sl] = neighbors.kneighbors(
            X[sl], min_samples)[0][:, -1]
    return core_distances


def _calculate_optics_order(X, min_samples, algorithm, leaf_size,
                            metric, metric_params, p, n_jobs, max_eps):
    n_samples = len(X)
    # Start all points as 'unprocessed' ##
    reachability_ = np.empty(n_samples)
    reachability_.fill(np.inf)
    predecessor_ = np.empty(n_samples, dtype=int)
    predecessor_.fill(-1)

    nbrs = NearestNeighbors(n_neighbors=min_samples,
                            algorithm=algorithm,
                            leaf_size=leaf_size,
                            metric=metric,
                            metric_params=metric_params,
                            p=p,
                            n_jobs=n_jobs)

    nbrs.fit(X)
    # Here we first do a kNN query for each point, this differs from
    # the original OPTICS that only used epsilon range queries.
    # TODO: handle working_memory somehow?
    core_distances_ = _compute_core_distances_(X=X, neighbors=nbrs,
                                               min_samples=min_samples,
                                               working_memory=None)
    # OPTICS puts an upper limit on these, use inf for undefined.
    core_distances_[core_distances_ > max_eps] = np.inf

    # Main OPTICS loop. Not parallelizable. The order that entries are
    # written to the 'ordering_' list is important!
    # Note that this implementation is O(n^2) theoretically, but
    # supposedly with very low constant factors.
    processed = np.zeros(X.shape[0], dtype=bool)
    ordering = np.zeros(X.shape[0], dtype=int)
    for ordering_idx in range(X.shape[0]):
        # Choose next based on smallest reachability distance
        # (And prefer smaller ids on ties, possibly np.inf!)
        index = np.where(processed == 0)[0]
        point = index[np.argmin(reachability_[index])]

        processed[point] = True
        ordering[ordering_idx] = point
        if core_distances_[point] != np.inf:
            _set_reach_dist(core_distances_=core_distances_,
                            reachability_=reachability_,
                            predecessor_=predecessor_,
                            point_index=point,
                            processed=processed, X=X, nbrs=nbrs,
                            metric=metric, metric_params=metric_params,
                            p=p, max_eps=max_eps)
    return ordering, core_distances_, reachability_, predecessor_


def _set_reach_dist(core_distances_, reachability_, predecessor_,
                    point_index, processed, X, nbrs, metric, metric_params,
                    p, max_eps):
    P = X[point_index:point_index + 1]
    # Assume that radius_neighbors is faster without distances
    # and we don't need all distances, nevertheless, this means
    # we may be doing some work twice.
    indices = nbrs.radius_neighbors(P, radius=max_eps,
                                    return_distance=False)[0]

    # Getting indices of neighbors that have not been processed
    unproc = np.compress((~np.take(processed, indices)).ravel(),
                         indices, axis=0)
    # Neighbors of current point are already processed.
    if not unproc.size:
        return

    # Only compute distances to unprocessed neighbors:
    if metric == 'precomputed':
        dists = X[point_index, unproc]
    else:
        _params = dict() if metric_params is None else metric_params.copy()
        if metric == 'minkowski' and 'p' not in _params:
            # the same logic as neighbors, p is ignored if explicitly set
            # in the dict params
            _params['p'] = p
        dists = pairwise_distances(P, np.take(X, unproc, axis=0),
                                   metric, n_jobs=None,
                                   **_params).ravel()

    rdists = np.maximum(dists, core_distances_[point_index])
    improved = np.where(rdists < np.take(reachability_, unproc))
    reachability_[unproc[improved]] = rdists[improved]
    predecessor_[unproc[improved]] = point_index


def _extract_optics_dbscan(reachability, core_distances, ordering, eps=0.5):
    """Performs DBSCAN extraction for an arbitrary epsilon.

    Extraction runs in linear time. Note that if the `max_eps` OPTICS
    parameter was set to < inf for extracting reachability and ordering
    arrays, DBSCAN extractions will be unstable for `eps` values close to
    `max_eps`. Setting `eps` < (`max_eps` / 5.0) will guarantee
    extraction parity with DBSCAN.

    Parameters
    ----------
    reachability : array, shape (n_samples,)
        Reachability distances calculated by OPTICS (`reachability_`)

    core_distances : array, shape (n_samples,)
        Distances at which points become core (`core_distances_`)

    ordering : array, shape (n_samples,)
        OPTICS ordered point indices (`ordering_`)

    eps : float, optional (default=0.5)
        DBSCAN `eps` parameter. Must be set to < `max_eps`. Results
        will be close to DBSCAN algorithm if `eps` is < (`max_eps` / 5)

    Returns
    -------
    core_sample_indices_ : array, shape (n_core_samples,)
        The indices of the core samples.

    labels_ : array, shape (n_samples,)
        The estimated labels.
    """
    n_samples = len(core_distances)
    is_core = np.zeros(n_samples, dtype=bool)
    labels = np.zeros(n_samples, dtype=int)

    far_reach = reachability > eps
    near_core = core_distances <= eps
    labels[ordering] = np.cumsum(far_reach[ordering] & near_core[ordering]) - 1
    labels[far_reach & ~near_core] = -1
    is_core[near_core] = True
    return np.arange(n_samples)[is_core], labels


def _extract_optics_sqlnk(reachability, ordering,
                          maxima_ratio=.75,
                          rejection_ratio=.7, similarity_threshold=.4,
                          significant_min=.003, min_cluster_size=.005,
                          min_maxima_ratio=.001):
    """Performs automatic cluster extraction for variable density data.
    All parameters will use the value present in the class instance if
    not provided.

    Parameters
    ----------
    reachability : array, shape (n_samples,)
        Reachability distances calculated by OPTICS (`reachability_`)

    ordering : array, shape (n_samples,)
        OPTICS ordered point indices (`ordering_`)

    maxima_ratio : float, optional (default=.75)
        The maximum ratio we allow of average height of clusters on the
        right and left to the local maxima in question. The higher the
        ratio, the more generous the algorithm is to preserving local
        minima, and the more cuts the resulting tree will have.

    rejection_ratio : float, optional (default=.7)
        Adjusts the fitness of the clustering. When the maxima_ratio is
        exceeded, determine which of the clusters to the left and right to
        reject based on rejection_ratio. Higher values will result in points
        being more readily classified as noise; conversely, lower values will
        result in more points being clustered.

    similarity_threshold : float, optional (default=.4)
        Used to check if nodes can be moved up one level, that is, if the
        new cluster created is too "similar" to its parent, given the
        similarity threshold. Similarity can be determined by 1) the size
        of the new cluster relative to the size of the parent node or
        2) the average of the reachability values of the new cluster
        relative to the average of the reachability values of the parent
        node. A lower value for the similarity threshold means less levels
        in the tree.

    significant_min : float, optional (default=.003)
        Sets a lower threshold on how small a significant maxima can be.

    min_cluster_size : int > 1 or float between 0 and 1 (default=0.005)
        Minimum number of samples in an OPTICS cluster, expressed as an
        absolute number or a fraction of the number of samples (rounded
        to be at least 2).

    min_maxima_ratio : float, optional (default=.001)
        Used to determine neighborhood size for minimum cluster membership.
        Each local maxima should be a largest value in a neighborhood
        of the `size min_maxima_ratio * len(X)` from left and right.

    Returns
    -------
    core_sample_indices_ : array, shape (n_core_samples,)
        The indices of the core samples.

    labels_ : array, shape (n_samples,)
        The estimated labels.
    """

    # Extraction wrapper
    # according to Ankerst M. et.al. 1999 (p. 5), for a small enough
    # generative distance epsilong, there should be more than one INF.
    if np.all(np.isinf(reachability)):
        raise ValueError("All reachability values are inf. Set a larger"
                         " max_eps.")
    normalization_factor = np.max(reachability[reachability < np.inf])
    reachability = reachability / normalization_factor
    reachability_plot = reachability[ordering].tolist()
    root_node = _sqlnk_cluster(reachability_plot, ordering,
                               maxima_ratio, rejection_ratio,
                               similarity_threshold, significant_min,
                               min_cluster_size, min_maxima_ratio)
    leaves = _get_leaves(root_node, [])
    # Start cluster id's at 0
    clustid = 0
    n_samples = len(reachability)
    is_core = np.zeros(n_samples, dtype=bool)
    labels = np.full(n_samples, -1, dtype=int)
    # Start all points as non-core noise
    for leaf in leaves:
        index = ordering[leaf.start:leaf.end]
        labels[index] = clustid
        is_core[index] = 1
        clustid += 1
    return np.arange(n_samples)[is_core], labels


def _sqlnk_cluster(reachability_plot, ordering,
                   maxima_ratio, rejection_ratio,
                   similarity_threshold, significant_min,
                   min_cluster_size, min_maxima_ratio):
    """Converts reachability plot to cluster tree and returns root node.

    Parameters
    ----------

    reachability_plot : list, required
        Reachability distances ordered by OPTICS ordering index.

    """

    min_neighborhood_size = 2
    if min_cluster_size <= 1:
        min_cluster_size = max(2, min_cluster_size * len(ordering))
    neighborhood_size = int(min_maxima_ratio * len(ordering))

    # Again, should this check < min_samples, should the parameter be public?
    if neighborhood_size < min_neighborhood_size:
        neighborhood_size = min_neighborhood_size

    local_maxima_points = _find_local_maxima(reachability_plot,
                                             neighborhood_size)
    root_node = _TreeNode(ordering, 0, len(ordering), None)
    _cluster_tree(root_node, None, local_maxima_points,
                  reachability_plot, ordering, min_cluster_size,
                  maxima_ratio, rejection_ratio,
                  similarity_threshold, significant_min)

    return root_node


class _TreeNode:
    # automatic cluster helper classes and functions
    def __init__(self, points, start, end, parent_node):
        self.points = points
        self.start = start
        self.end = end
        self.parent_node = parent_node
        self.children = []
        self.split_point = -1


def _is_local_maxima(index, reachability_plot, neighborhood_size):
    right_idx = slice(index + 1, index + neighborhood_size + 1)
    left_idx = slice(max(1, index - neighborhood_size - 1), index)
    return (np.all(reachability_plot[index] >= reachability_plot[left_idx]) and
            np.all(reachability_plot[index] >= reachability_plot[right_idx]))


def _find_local_maxima(reachability_plot, neighborhood_size):
    local_maxima_points = {}
    # 1st and last points on Reachability Plot are not taken
    # as local maxima points
    for i in range(1, len(reachability_plot) - 1):
        # if the point is a local maxima on the reachability plot with
        # regard to neighborhood_size, insert it into priority queue and
        # maxima list
        if (reachability_plot[i] > reachability_plot[i - 1] and
            reachability_plot[i] >= reachability_plot[i + 1] and
            _is_local_maxima(i, np.array(reachability_plot),
                             neighborhood_size) == 1):
            local_maxima_points[i] = reachability_plot[i]

    return sorted(local_maxima_points,
                  key=local_maxima_points.__getitem__, reverse=True)


def _cluster_tree(node, parent_node, local_maxima_points,
                  reachability_plot, reachability_ordering,
                  min_cluster_size, maxima_ratio, rejection_ratio,
                  similarity_threshold, significant_min):
    """Recursively builds cluster tree to hold hierarchical cluster structure

    node is a node or the root of the tree in the first call
    parent_node is parent node of N or None if node is root of the tree
    local_maxima_points is list of local maxima points sorted in
    descending order of reachability
    """

    if len(local_maxima_points) == 0:
        return  # parent_node is a leaf

    # take largest local maximum as possible separation between clusters
    s = local_maxima_points[0]
    node.split_point = s
    local_maxima_points = local_maxima_points[1:]

    # create two new nodes and add to list of nodes
    node_1 = _TreeNode(reachability_ordering[node.start:s],
                       node.start, s, node)
    node_2 = _TreeNode(reachability_ordering[s + 1:node.end],
                       s + 1, node.end, node)
    local_max_1 = []
    local_max_2 = []

    for i in local_maxima_points:
        if i < s:
            local_max_1.append(i)
        if i > s:
            local_max_2.append(i)

    node_list = []
    node_list.append((node_1, local_max_1))
    node_list.append((node_2, local_max_2))

    if reachability_plot[s] < significant_min:
        node.split_point = -1
        # if split_point is not significant, ignore this split and continue
        return

    # only check a certain ratio of points in the child
    # nodes formed to the left and right of the maxima
    # ...should check_ratio be a user settable parameter?
    check_ratio = .8
    check_value_1 = int(np.round(check_ratio * len(node_1.points)))
    check_value_2 = int(np.round(check_ratio * len(node_2.points)))
    avg_reach1 = np.mean(reachability_plot[(node_1.end -
                                            check_value_1):node_1.end])
    avg_reach2 = np.mean(reachability_plot[node_2.start:(node_2.start
                                                         + check_value_2)])

    if ((avg_reach1 / maxima_ratio) > reachability_plot[s] or
            (avg_reach2 / maxima_ratio) > reachability_plot[s]):

        if (avg_reach1 / rejection_ratio) < reachability_plot[s]:
            # reject node 2
            node_list.remove((node_2, local_max_2))
        if (avg_reach2 / rejection_ratio) < reachability_plot[s]:
            # reject node 1
            node_list.remove((node_1, local_max_1))
        if ((avg_reach1 / rejection_ratio) >= reachability_plot[s] and
                (avg_reach2 / rejection_ratio) >= reachability_plot[s]):
            # since split_point is not significant,
            # ignore this split and continue (reject both child nodes)
            node.split_point = -1
            _cluster_tree(node, parent_node, local_maxima_points,
                          reachability_plot, reachability_ordering,
                          min_cluster_size, maxima_ratio, rejection_ratio,
                          similarity_threshold, significant_min)
            return

    # remove clusters that are too small
    if (len(node_1.points) < min_cluster_size and
            node_list.count((node_1, local_max_1)) > 0):
        # cluster 1 is too small
        node_list.remove((node_1, local_max_1))
    if (len(node_2.points) < min_cluster_size and
            node_list.count((node_2, local_max_2)) > 0):
        # cluster 2 is too small
        node_list.remove((node_2, local_max_2))
    if not node_list:
        # parent_node will be a leaf
        node.split_point = -1
        return

    # Check if nodes can be moved up one level - the new cluster created
    # is too "similar" to its parent, given the similarity threshold.
    bypass_node = 0
    if parent_node is not None:
        if ((node.end - node.start) / (parent_node.end - parent_node.start) >
                similarity_threshold):

            parent_node.children.remove(node)
            bypass_node = 1

    for nl in node_list:
        if bypass_node == 1:
            parent_node.children.append(nl[0])
            _cluster_tree(nl[0], parent_node, nl[1],
                          reachability_plot, reachability_ordering,
                          min_cluster_size, maxima_ratio, rejection_ratio,
                          similarity_threshold, significant_min)
        else:
            node.children.append(nl[0])
            _cluster_tree(nl[0], node, nl[1], reachability_plot,
                          reachability_ordering, min_cluster_size,
                          maxima_ratio, rejection_ratio,
                          similarity_threshold, significant_min)


def _get_leaves(node, arr):
    if node is not None:
        if node.split_point == -1:
            arr.append(node)
        for n in node.children:
            _get_leaves(n, arr)
    return arr
