"""Nearest Neighbors graph functions"""

# Author: Jake Vanderplas <vanderplas@astro.washington.edu>
#
# License: BSD 3 clause (C) INRIA, University of Amsterdam

from .base import KNeighborsMixin, RadiusNeighborsMixin
from .unsupervised import NearestNeighbors


def _check_params(X, metric, p, metric_params):
    """Check the validity of the input parameters"""
    params = zip(['metric', 'p', 'metric_params'],
                 [metric, p, metric_params])
    est_params = X.get_params()
    for param_name, func_param in params:
        if func_param != est_params[param_name]:
            raise ValueError(
                "Got %s for %s, while the estimator has %s for "
                "the same parameter." % (
                    func_param, param_name, est_params[param_name]))


def _query_include_self(X, include_self):
    """Return the query based on include_self param"""
    if include_self:
        query = X._fit_X
    else:
        query = None

    return query


def kneighbors_graph(X, n_neighbors, mode='connectivity', metric='minkowski',
                     p=2, metric_params=None, include_self=False, n_jobs=1):
    """Computes the (weighted) graph of k-Neighbors for points in X

    Read more in the :ref:`User Guide <unsupervised_neighbors>`.

    Parameters
    ----------
    X : array-like or BallTree, shape = [n_samples, n_features]
        Sample data, in the form of a numpy array or a precomputed
        :class:`BallTree`.

    n_neighbors : int
        Number of neighbors for each sample.

    mode : {'connectivity', 'distance'}, optional
        Type of returned matrix: 'connectivity' will return the connectivity
        matrix with ones and zeros, and 'distance' will return the distances
        between neighbors according to the given metric.

    metric : string, default 'minkowski'
        The distance metric used to calculate the k-Neighbors for each sample
        point. The DistanceMetric class gives a list of available metrics.
        The default distance is 'euclidean' ('minkowski' metric with the p
        param equal to 2.)

    include_self: bool, default=False.
        Whether or not to mark each sample as the first nearest neighbor to
        itself. If `None`, then True is used for mode='connectivity' and False
        for mode='distance' as this will preserve backwards compatibilty.

    p : int, default 2
        Power parameter for the Minkowski metric. When p = 1, this is
        equivalent to using manhattan_distance (l1), and euclidean_distance
        (l2) for p = 2. For arbitrary p, minkowski_distance (l_p) is used.

    metric_params: dict, optional
        additional keyword arguments for the metric function.

    n_jobs : int, optional (default = 1)
        The number of parallel jobs to run for neighbors search.
        If ``-1``, then the number of jobs is set to the number of CPU cores.

    Returns
    -------
    A : sparse matrix in CSR format, shape = [n_samples, n_samples]
        A[i, j] is assigned the weight of edge that connects i to j.

    Examples
    --------
    >>> X = [[0], [3], [1]]
    >>> from sklearn.neighbors import kneighbors_graph
    >>> A = kneighbors_graph(X, 2, mode='connectivity', include_self=True)
    >>> A.toarray()
    array([[ 1.,  0.,  1.],
           [ 0.,  1.,  1.],
           [ 1.,  0.,  1.]])

    See also
    --------
    radius_neighbors_graph
    """
    if not isinstance(X, KNeighborsMixin):
        X = NearestNeighbors(n_neighbors, metric=metric, p=p,
                             metric_params=metric_params, n_jobs=n_jobs).fit(X)
    else:
        _check_params(X, metric, p, metric_params)

    query = _query_include_self(X, include_self)
    return X.kneighbors_graph(X=query, n_neighbors=n_neighbors, mode=mode)


def radius_neighbors_graph(X, radius, mode='connectivity', metric='minkowski',
                           p=2, metric_params=None, include_self=False, n_jobs=1):
    """Computes the (weighted) graph of Neighbors for points in X

    Neighborhoods are restricted the points at a distance lower than
    radius.

    Read more in the :ref:`User Guide <unsupervised_neighbors>`.

    Parameters
    ----------
    X : array-like or BallTree, shape = [n_samples, n_features]
        Sample data, in the form of a numpy array or a precomputed
        :class:`BallTree`.

    radius : float
        Radius of neighborhoods.

    mode : {'connectivity', 'distance'}, optional
        Type of returned matrix: 'connectivity' will return the connectivity
        matrix with ones and zeros, and 'distance' will return the distances
        between neighbors according to the given metric.

    metric : string, default 'minkowski'
        The distance metric used to calculate the neighbors within a
        given radius for each sample point. The DistanceMetric class
        gives a list of available metrics. The default distance is
        'euclidean' ('minkowski' metric with the param equal to 2.)

    include_self: bool, default=False
        Whether or not to mark each sample as the first nearest neighbor to
        itself. If `None`, then True is used for mode='connectivity' and False
        for mode='distance' as this will preserve backwards compatibilty.

    p : int, default 2
        Power parameter for the Minkowski metric. When p = 1, this is
        equivalent to using manhattan_distance (l1), and euclidean_distance
        (l2) for p = 2. For arbitrary p, minkowski_distance (l_p) is used.

    metric_params: dict, optional
        additional keyword arguments for the metric function.

    n_jobs : int, optional (default = 1)
        The number of parallel jobs to run for neighbors search.
        If ``-1``, then the number of jobs is set to the number of CPU cores.

    Returns
    -------
    A : sparse matrix in CSR format, shape = [n_samples, n_samples]
        A[i, j] is assigned the weight of edge that connects i to j.

    Examples
    --------
    >>> X = [[0], [3], [1]]
    >>> from sklearn.neighbors import radius_neighbors_graph
    >>> A = radius_neighbors_graph(X, 1.5, mode='connectivity', include_self=True)
    >>> A.toarray()
    array([[ 1.,  0.,  1.],
           [ 0.,  1.,  0.],
           [ 1.,  0.,  1.]])

    See also
    --------
    kneighbors_graph
    """
    if not isinstance(X, RadiusNeighborsMixin):
        X = NearestNeighbors(radius=radius, metric=metric, p=p,
                             metric_params=metric_params, n_jobs=n_jobs).fit(X)
    else:
        _check_params(X, metric, p, metric_params)

    query = _query_include_self(X, include_self)
    return X.radius_neighbors_graph(query, radius, mode)
