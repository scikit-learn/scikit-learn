# -*- coding: utf-8 -*-
"""Ordering Points To Identify the Clustering Structure (OPTICS)

These routines execute the OPTICS algorithm, and implement various
cluster extraction methods of the ordered list.

Authors: Shane Grigsby <refuge@rocktalus.com>
         Adrin Jalali <adrinjalali@gmail.com>
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


class OPTICS(BaseEstimator, ClusterMixin):
    """Estimate clustering structure from vector array

    OPTICS: Ordering Points To Identify the Clustering Structure Closely
    related to DBSCAN, finds core sample of high density and expands clusters
    from them [1]_. Unlike DBSCAN, keeps cluster hierarchy for a variable
    neighborhood radius. Better suited for usage on large datasets than the
    current sklearn implementation of DBSCAN.

    Clusters are then extracted using a DBSCAN like method [1]_.

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
        as in the same neighborhood. Default value of ``np.inf`` will identify
        clusters across all scales; reducing ``max_eps`` will result in
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

    cluster_method : string, optional (default='dbscan')
        The extraction method used to extract clusters using the calculated
        reachability and ordering. Possible values are "xi" and "dbscan".

    eps : float, optional (default=0.5)
        The maximum distance between two samples for them to be considered
        as in the same neighborhood. Used ony when `cluster_method='dbscan'`.

    xi : float, between 0 and 1
        For cluster_method='xi'. Determines the minimum steepness on the
        reachability plot that constitutes a cluster boundary. For example, an
        upwards point in the reachability plot is defined by the ratio from one
        point to its successor being at most 1-xi.

    min_cluster_size : int > 1 or float between 0 and 1 (default=0.005)
        Minimum number of samples in an OPTICS cluster, expressed as an
        absolute number or a fraction of the number of samples (rounded
        to be at least 2).
        Used only when ``cluster_method='xi'``.

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

    .. [2] Schubert, Erich, Michael Gertz.
       "Improving the Cluster Structure Extracted from OPTICS Plots." Proc. of
       the Conference "Lernen, Wissen, Daten, Analysen" (LWDA) (2018): 318-329.

    """

    def __init__(self, min_samples=5, max_eps=np.inf, metric='minkowski', p=2,
                 metric_params=None, cluster_method='xi', eps=0.5, xi=0.05,
                 min_cluster_size=.005, algorithm='auto', leaf_size=30,
                 n_jobs=None):

        self.max_eps = max_eps
        self.min_samples = min_samples
        self.min_cluster_size = min_cluster_size
        self.algorithm = algorithm
        self.metric = metric
        self.metric_params = metric_params
        self.p = p
        self.leaf_size = leaf_size
        self.cluster_method = cluster_method
        self.eps = eps
        self.xi = xi
        self.n_jobs = n_jobs

    def fit(self, X, y=None):
        """Perform OPTICS clustering

        Extracts an ordered list of points and reachability distances, and
        performs initial clustering using ``max_eps`` distance specified at
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

        if self.min_samples > n_samples:
            raise ValueError("Number of training samples (n_samples=%d) must "
                             "be greater than min_samples (min_samples=%d) "
                             "used for clustering." %
                             (n_samples, self.min_samples))

        if self.cluster_method not in ['dbscan', 'xi']:
            raise ValueError("cluster_method should be one of"
                             " 'dbscan' or 'xi' but is %s" %
                             self.cluster_method)

        if self.cluster_method == 'dbscan':
            if self.eps > self.max_eps:
                raise ValueError('Specify an epsilon smaller than %s. Got %s.'
                                 % (self.max_eps, self.eps))

            if self.eps * 5.0 > self.max_eps * 1.05:
                warnings.warn(
                    "Warning, max_eps (%s) is close to eps (%s): "
                    "Output may be unstable." % (self.max_eps, self.eps),
                    RuntimeWarning, stacklevel=2)
                # Stability warning is documented in cluster_optics_dbscan
                # method...

        (self.ordering_, self.core_distances_, self.reachability_,
         self.predecessor_) = compute_optics_graph(
             X=X, min_samples=self.min_samples, algorithm=self.algorithm,
             leaf_size=self.leaf_size, metric=self.metric,
             metric_params=self.metric_params, p=self.p, n_jobs=self.n_jobs,
             max_eps=self.max_eps)
        # Start all points as noise ##
        self.labels_ = np.full(n_samples, -1, dtype=int)

        # Extract clusters from the calculated orders and reachability
        if self.cluster_method == 'xi':
            # TODO: _ is the to-be-cluster_
            labels_, _ = cluster_optics_xi(self.reachability_, self.ordering_,
                                           self.min_samples,
                                           self.min_cluster_size, self.xi)
        elif self.cluster_method == 'dbscan':
            labels_ = cluster_optics_dbscan(self.reachability_,
                                            self.core_distances_,
                                            self.ordering_,
                                            self.eps)

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


def compute_optics_graph(X, min_samples, max_eps, metric, p, metric_params,
                         algorithm, leaf_size, n_jobs):
    """Computes the OPTICS reachability graph.

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
    ordering_ : array, shape (n_samples,)
        The cluster ordered list of sample indices.

    core_distances_ : array, shape (n_samples,)
        Distance at which each sample becomes a core point, indexed by object
        order. Points which will never be core have a distance of inf. Use
        ``clust.core_distances_[clust.ordering_]`` to access in cluster order.

    reachability_ : array, shape (n_samples,)
        Reachability distances per sample, indexed by object order. Use
        ``clust.reachability_[clust.ordering_]`` to access in cluster order.

    predecessor_ : array, shape (n_samples,)
        Point that a sample was reached from, indexed by object order.
        Seed points have a predecessor of -1.

    References
    ----------
    .. [1] Ankerst, Mihael, Markus M. Breunig, Hans-Peter Kriegel,
       and Jörg Sander. "OPTICS: ordering points to identify the clustering
       structure." ACM SIGMOD Record 28, no. 2 (1999): 49-60.
    """
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
    if np.all(np.isinf(reachability_)):
        warnings.warn("All reachability values are inf. Set a larger"
                      " max_eps or all data will be considered outliers.",
                      UserWarning)
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


def cluster_optics_dbscan(reachability, core_distances, ordering, eps=0.5):
    """Performs DBSCAN extraction for an arbitrary epsilon.

    Extracting the clusters runs in linear time. Note that if the `max_eps`
    OPTICS parameter was set to < inf for extracting reachability and ordering
    arrays, DBSCAN extractions will be unstable for `eps` values close to
    `max_eps`. Setting `eps` < (`max_eps` / 5.0) will guarantee extraction
    parity with DBSCAN.

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
    labels_ : array, shape (n_samples,)
        The estimated labels.

    """
    n_samples = len(core_distances)
    labels = np.zeros(n_samples, dtype=int)

    far_reach = reachability > eps
    near_core = core_distances <= eps
    labels[ordering] = np.cumsum(far_reach[ordering] & near_core[ordering]) - 1
    labels[far_reach & ~near_core] = -1
    return labels


def cluster_optics_xi(reachability, ordering, min_samples, min_cluster_size,
                      xi):
    """Automatically extract clusters according to the Xi-steep method.

    Parameters
    ----------
    reachability : array, shape (n_samples,)
        Reachability distances calculated by OPTICS (`reachability_`)

    ordering : array, shape (n_samples,)
        OPTICS ordered point indices (`ordering_`)

    min_samples : integer
       The same as the min_samples given to OPTICS. Up and down steep
       regions can't have more then `min_samples` consecutive non-steep
       points.

    min_cluster_size : int > 1 or float between 0 and 1
        Minimum number of samples in an OPTICS cluster, expressed as an
        absolute number or a fraction of the number of samples (rounded
        to be at least 2).

    xi : float, between 0 and 1
        Determines the minimum steepness on the reachability plot that
        constitutes a cluster boundary. For example, an upwards point in the
        reachability plot is defined by the ratio from one point to its
        successor being at most 1-xi.

    Returns
    -------
    labels : array, shape (n_samples)
        The labels assigned to samples. Points which are not included
        in any cluster are labeled as -1.

    clusters : list
        The list of clusters in the form of (start, end) tuples,
        with all indices inclusive. The clusters are ordered in a way that
        larger clusters encompassing smaller clusters, come after those
        smaller clusters.
    """
    clusters = _xi_cluster(reachability[ordering], xi, min_samples,
                           min_cluster_size)
    labels = _extract_xi_labels(ordering, clusters)
    print("labels: ", labels[ordering])
    return labels, clusters


def _steep_upward(reachability_plot, p, xi_complement):
    """Check if point p is a xi steep up area (definition 9).
    xi_complement is the inverse xi, i.e. `1 - xi`"""

    if np.isinf(reachability_plot[p]) and np.isinf(reachability_plot[p + 1]):
        return False

    return reachability_plot[p] <= reachability_plot[p + 1] * xi_complement


def _steep_downward(reachability_plot, p, xi_complement):
    """Check if point p is a xi steep down area (definition 9).
    xi_complement is the inverse xi, i.e. `1 - xi`"""

    if np.isinf(reachability_plot[p]) and np.isinf(reachability_plot[p + 1]):
        return False

    return reachability_plot[p] * xi_complement >= reachability_plot[p + 1]


class _Area:
    """An (upward or downward) area.

    Attributes
    ----------
    start : integer
        The start of the region.

    end : integer
        The end of the region.

    maximum : float
        The maximum reachability in this region, which is the
        start of the region for a downward area and the end of
        the region for an upward area.

    mib : float
        Maximum in between value, i.e. the maximum value between
        the end of a steep down area and the current index. It is
        irrelevant for steep up areas.
    """
    def __init__(self, start, end, maximum, mib):
        self.start = start
        self.end = end
        # maximum is the beginning of a downward area
        # and the end of an upward area
        self.maximum = maximum
        # maximum in between (ref Sec. 4.3.2)
        self.mib = mib

    def __repr__(self):
        return ("start: %s, end: %s, max: %4g, mib: %4g" %
                (self.start, self.end, self.maximum, self.mib))


def _extend_downward(reachability_plot, start, xi_complement, min_samples,
                     n_samples):
    """Extend the downward area until it's maximal.

    Parameters
    ----------
    reachability_plot : array, shape (n_samples)
        The reachability plot, i.e. reachability ordered according to
        the calculated ordering, all computed by OPTICS.

    start : integer
        The start of the downward region.
    
    xi_complement : float, between 0 and 1
        The inverse xi, i.e. `1 - xi`

    min_samples : integer
       The same as the min_samples given to OPTICS. Up and down steep
       regions can't have more then `min_samples` consecutive non-steep
       points.

    n_samples : integer
        Total number of samples.

    Returns
    -------
    index : integer
        The current index iterating over all the samples.

    end : integer
        The end of the downward region, which can be behind the index.
    """
    # print("extend down start")
    non_downward_points = 0
    index = start
    end = start
    # find a maximal downward area assuming the first point is xi-downward
    # and there's an inf attached to the end of the reachability plot.
    while index < n_samples:
        # print("index", index)
        # print("r", reachability_plot[index], "r + 1",
        #       reachability_plot[index + 1])
        if _steep_downward(reachability_plot, index, xi_complement):
            non_downward_points = 0
            end = index
            # print("steep, end:", end)
        elif reachability_plot[index] >= reachability_plot[index + 1]:
            # print("just down")
            # it's not a steep downward point, but still goes down.
            non_downward_points += 1
            # region should include no more than min_samples consecutive
            # non steep downward points.
            if non_downward_points > min_samples:
                # print("too many downward")
                break
        else:
            # this is not a steep down area anymore
            # print("not going down, index:", index)
            return index, end

        index += 1
    # print("extend end")
    return index + 1, end


def _extend_upward(reachability_plot, start, xi_complement, min_samples,
                   n_samples):
    """Extend the upward area until it's maximal.

    Parameters
    ----------
    reachability_plot : array, shape (n_samples)
        The reachability plot, i.e. reachability ordered according to
        the calculated ordering, all computed by OPTICS.

    start : integer
        The start of the upward region.

    xi_complement : float, between 0 and 1
        The inverse xi, i.e. `1 - xi`

    min_samples : integer
       The same as the min_samples given to OPTICS. Up and down steep
       regions can't have more then `min_samples` consecutive non-steep
       points.

    n_samples : integer
        Total number of samples.

    Returns
    -------
    index : integer
        The current index iterating over all the samples.

    end : integer
        The end of the upward region, which can be behind the index.
    """
    # print("extend up start")
    non_upward_points = 0
    index = start
    end = start
    # find a maximal upward area
    while index < n_samples:
        # print("index", index)
        # print("r", reachability_plot[index], "r + 1",
        #       reachability_plot[index + 1])
        if _steep_upward(reachability_plot, index, xi_complement):
            non_upward_points = 0
            end = index
            # print("steep, end:", end)
        elif reachability_plot[index] <= reachability_plot[index + 1]:
            # print("just up")
            # it's not a steep upward point, but still goes up.
            non_upward_points += 1
            # region should include no more than min_samples consecutive
            # non steep upward points.
            if non_upward_points > min_samples:
                # print("non upward")
                break
        else:
            # print("not going up")
            return min(index, n_samples - 1), end

        index += 1
    # print("extend end")
    return min(index + 1, n_samples - 1), end


def _update_filter_sdas(sdas, mib, xi_complement):
    """Update steep down areas (SDAs) using the new
    maximum in between (mib) value, and the given inverse xi, i.e. `1 - xi`
    """
    res = [sda for sda in sdas if mib <= sda.maximum * xi_complement]
    for sda in res:
        sda.mib = max(sda.mib, mib)
    return res


def _xi_cluster(reachability_plot, xi, min_samples, min_cluster_size):
    """Automatically extract clusters according to the Xi-steep method.

    This is rouphly an implementation of Figure 19 of the OPTICS paper.

    Parameters
    ----------
    reachability_plot : array, shape (n_samples)
        The reachability plot, i.e. reachability ordered according to
        the calculated ordering, all computed by OPTICS.

    xi : float, between 0 and 1
        Determines the minimum steepness on the reachability plot that
        constitutes a cluster boundary. For example, an upwards point in the
        reachability plot is defined by the ratio from one point to its
        successor being at most 1-xi.

    min_samples : integer
       The same as the min_samples given to OPTICS. Up and down steep
       regions can't have more then `min_samples` consecutive non-steep
       points.

    min_cluster_size : int > 1 or float between 0 and 1
        Minimum number of samples in an OPTICS cluster, expressed as an
        absolute number or a fraction of the number of samples (rounded
        to be at least 2).

    Returns
    -------
    clusters : list
        The list of clusters in the form of (start, end) tuples,
        with all indices inclusive. The clusters are ordered in a way that
        larger clusters encompassing smaller clusters, come after those
        smaller clusters.
    """
    
    for i in range(len(reachability_plot) - 1):
        print("%d\t%d %d %.6f\n" % (
            i, _steep_downward(reachability_plot, i, 1 - xi),
            _steep_upward(reachability_plot, i, 1 - xi),
            reachability_plot[i]))
    print(reachability_plot[-1])

    # all indices are inclusive (specially at the end)
    n_samples = len(reachability_plot)
    # add an inf to the end of reachability plot
    # but process the data only until the last actual point
    # this is fine since the last point is considered upward anyway
    reachability_plot = np.hstack((reachability_plot, np.inf))

    if min_cluster_size <= 1:
        min_cluster_size = max(2, min_cluster_size * n_samples)

    xi_complement = 1 - xi
    sdas = [] # steep down areas, introduced in section 4.3.2 of the paper
    clusters = []
    index = 0
    mib = 0.  # maximum in between
    while index < n_samples - 1:
        # print("index", index)
        # print("r", reachability_plot[index])
        mib = max(mib, reachability_plot[index])
        # print("mib up there:", mib)

        # check if a steep downward area starts
        if _steep_downward(reachability_plot, index, xi_complement):
            print("steep downward")
            # print("sdas", sdas)
            # print("filter mib:", mib)
            sdas = _update_filter_sdas(sdas, mib, xi_complement)
            # print("sdas", sdas)
            D_start = index
            index, end = _extend_downward(reachability_plot, D_start,
                                          xi_complement,
                                          min_samples, n_samples)
            D = _Area(start=D_start, end=end,
                      maximum=reachability_plot[D_start], mib=0.)
            print("D", D, "r.s %.4g" % reachability_plot[D.start],
                  "r.e %.4g" % reachability_plot[D.end])
            sdas.append(D)
            mib = reachability_plot[index]

        elif _steep_upward(reachability_plot, index, xi_complement):
            print("steep upward")
            # print("sdas", sdas)
            # print("filter mib:", mib)
            sdas = _update_filter_sdas(sdas, mib, xi_complement)
            # print("sdas", sdas)
            U_start = index
            index, end = _extend_upward(reachability_plot, U_start,
                                        xi_complement,
                                        min_samples, n_samples)
            U = _Area(start=U_start, end=end, maximum=reachability_plot[end],
                      mib=-1)
            # if np.isinf(reachability_plot[index + 1]):
            #     U.maximum = np.inf
            #     index += 1
            # print("U", U, "r.s %.4g" % reachability_plot[U.start],
            #       "r.e %.4g" % reachability_plot[U.end])
            mib = reachability_plot[index]
            # print('mib %.4g' % mib)
            # print(sdas)

            U_clusters = []
            for D in sdas:
                c_start = D.start
                c_end = min(U.end, n_samples - 1)
                print("D", D, "U", U)
                print("start", c_start, "end", c_end)

                # line (**)
                if reachability_plot[c_end + 1] * xi_complement < D.mib:
                    continue

                # 3.b
                if D.mib > mib * xi_complement:
                    continue
                # print("3b pass")

                # 4
                if D.maximum * xi_complement >= reachability_plot[c_end + 1]:
                    while (reachability_plot[c_start + 1] >
                           reachability_plot[c_end + 1]
                           and c_start < c_end):
                        c_start += 1
                elif reachability_plot[c_end] * xi_complement >= D.maximum:
                    print("d max: ", D.maximum)
                    while (reachability_plot[c_end] > D.maximum
                           and c_end > c_start):
                        c_end -= 1
                print('after 4', c_start, c_end)

                #if _steep_upward(reachability_plot, index - 1, xi_complement):
                #    c_end -= 1
                # print('check last point', c_end, 'index', index)

                # 3.a
                if c_end - c_start + 1 < min_cluster_size:
                    continue
                # print('min pts pass')

                # 1
                if c_start > D.end:
                    continue

                # 2
                if c_end < U.start:
                    continue

                U_clusters.append((c_start, c_end))
                print('U clusters', U_clusters)

            # add smaller clusters first.
            U_clusters.reverse()
            clusters.extend(U_clusters)
            print("set of clusters:", clusters)

        else:
            # print("just else", index)
            index += 1

    return clusters


def _extract_xi_labels(ordering, clusters):
    """Extracts the labels from the clusters returned by `_xi_cluster`.
    We rely on the fact that clusters are stored
    with the smaller clusters coming before the larger ones.

    Parameters
    ----------
    ordering : array, shape (n_samples)
        The ordering of points calculated by OPTICS

    clusters : list
        List of clusters i.e. (start, end) tuples,
        as returned by `_xi_cluster`.

    Returns
    -------
    labels : array, shape (n_samples)
    """

    labels = np.zeros(len(ordering), dtype=np.int)
    label = 1
    for c in clusters:
        if all(labels[c[0]:(c[1] + 1)] == 0):
            labels[c[0]:(c[1] + 1)] = label
            label += 1
    labels[ordering] = labels - 1
    return labels
