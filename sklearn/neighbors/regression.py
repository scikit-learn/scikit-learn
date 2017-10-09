"""Nearest Neighbor Regression"""

# Authors: Jake Vanderplas <vanderplas@astro.washington.edu>
#          Fabian Pedregosa <fabian.pedregosa@inria.fr>
#          Alexandre Gramfort <alexandre.gramfort@inria.fr>
#          Sparseness support by Lars Buitinck
#          Multi-output support by Arnaud Joly <a.joly@ulg.ac.be>
#
# License: BSD 3 clause (C) INRIA, University of Amsterdam

import numpy as np

from .base import _get_weights, _check_weights, NeighborsBase, KNeighborsMixin
from .base import RadiusNeighborsMixin, SupervisedFloatMixin
from ..base import RegressorMixin
from ..utils import check_array


class KNeighborsRegressor(NeighborsBase, KNeighborsMixin,
                          SupervisedFloatMixin,
                          RegressorMixin):
    """Regression based on k-nearest neighbors.

    The target is predicted by local interpolation of the targets
    associated of the nearest neighbors in the training set.

    Read more in the :ref:`User Guide <regression>`.

    Parameters
    ----------
    n_neighbors : int, optional (default = 5)
        Number of neighbors to use by default for :meth:`kneighbors` queries.

    weights : str or callable
        weight function used in prediction.  Possible values:

        - 'uniform' : uniform weights.  All points in each neighborhood
          are weighted equally.
        - 'distance' : weight points by the inverse of their distance.
          in this case, closer neighbors of a query point will have a
          greater influence than neighbors which are further away.
        - [callable] : a user-defined function which accepts an
          array of distances, and returns an array of the same shape
          containing the weights.

        Uniform weights are used by default.

    algorithm : {'auto', 'ball_tree', 'kd_tree', 'brute'}, optional
        Algorithm used to compute the nearest neighbors:

        - 'ball_tree' will use :class:`BallTree`
        - 'kd_tree' will use :class:`KDTree`
        - 'brute' will use a brute-force search.
        - 'auto' will attempt to decide the most appropriate algorithm
          based on the values passed to :meth:`fit` method.

        Note: fitting on sparse input will override the setting of
        this parameter, using brute force.

    leaf_size : int, optional (default = 30)
        Leaf size passed to BallTree or KDTree.  This can affect the
        speed of the construction and query, as well as the memory
        required to store the tree.  The optimal value depends on the
        nature of the problem.

    p : integer, optional (default = 2)
        Power parameter for the Minkowski metric. When p = 1, this is
        equivalent to using manhattan_distance (l1), and euclidean_distance
        (l2) for p = 2. For arbitrary p, minkowski_distance (l_p) is used.

    metric : string or callable, default 'minkowski'
        the distance metric to use for the tree.  The default metric is
        minkowski, and with p=2 is equivalent to the standard Euclidean
        metric. See the documentation of the DistanceMetric class for a
        list of available metrics.

    metric_params : dict, optional (default = None)
        Additional keyword arguments for the metric function.

    n_jobs : int, optional (default = 1)
        The number of parallel jobs to run for neighbors search.
        If ``-1``, then the number of jobs is set to the number of CPU cores.
        Doesn't affect :meth:`fit` method.

    Examples
    --------
    >>> X = [[0], [1], [2], [3]]
    >>> y = [0, 0, 1, 1]
    >>> from sklearn.neighbors import KNeighborsRegressor
    >>> neigh = KNeighborsRegressor(n_neighbors=2)
    >>> neigh.fit(X, y) # doctest: +ELLIPSIS
    KNeighborsRegressor(...)
    >>> print(neigh.predict([[1.5]]))
    [ 0.5]

    See also
    --------
    NearestNeighbors
    RadiusNeighborsRegressor
    KNeighborsClassifier
    RadiusNeighborsClassifier

    Notes
    -----
    See :ref:`Nearest Neighbors <neighbors>` in the online documentation
    for a discussion of the choice of ``algorithm`` and ``leaf_size``.

    .. warning::

       Regarding the Nearest Neighbors algorithms, if it is found that two
       neighbors, neighbor `k+1` and `k`, have identical distances but
       different labels, the results will depend on the ordering of the
       training data.

    https://en.wikipedia.org/wiki/K-nearest_neighbor_algorithm
    """

    def __init__(self, n_neighbors=5, weights='uniform',
                 algorithm='auto', leaf_size=30,
                 p=2, metric='minkowski', metric_params=None, n_jobs=1,
                 **kwargs):
        self._init_params(n_neighbors=n_neighbors,
                          algorithm=algorithm,
                          leaf_size=leaf_size, metric=metric, p=p,
                          metric_params=metric_params, n_jobs=n_jobs, **kwargs)
        self.weights = _check_weights(weights)

    def predict(self, X):
        """Predict the target for the provided data

        Parameters
        ----------
        X : array-like, shape (n_query, n_features), \
                or (n_query, n_indexed) if metric == 'precomputed'
            Test samples.

        Returns
        -------
        y : array of int, shape = [n_samples] or [n_samples, n_outputs]
            Target values
        """
        X = check_array(X, accept_sparse='csr')

        neigh_dist, neigh_ind = self.kneighbors(X)

        weights = _get_weights(neigh_dist, self.weights)

        _y = self._y
        if _y.ndim == 1:
            _y = _y.reshape((-1, 1))

        if weights is None:
            y_pred = np.mean(_y[neigh_ind], axis=1)
        else:
            y_pred = np.empty((X.shape[0], _y.shape[1]), dtype=np.float64)
            denom = np.sum(weights, axis=1)

            for j in range(_y.shape[1]):
                num = np.sum(_y[neigh_ind, j] * weights, axis=1)
                y_pred[:, j] = num / denom

        if self._y.ndim == 1:
            y_pred = y_pred.ravel()

        return y_pred


class RadiusNeighborsRegressor(NeighborsBase, RadiusNeighborsMixin,
                               SupervisedFloatMixin,
                               RegressorMixin):
    """Regression based on neighbors within a fixed radius.

    The target is predicted by local interpolation of the targets
    associated of the nearest neighbors in the training set.

    Read more in the :ref:`User Guide <regression>`.

    Parameters
    ----------
    radius : float, optional (default = 1.0)
        Range of parameter space to use by default for :meth:`radius_neighbors`
        queries.

    weights : str or callable
        weight function used in prediction.  Possible values:

        - 'uniform' : uniform weights.  All points in each neighborhood
          are weighted equally.
        - 'distance' : weight points by the inverse of their distance.
          in this case, closer neighbors of a query point will have a
          greater influence than neighbors which are further away.
        - [callable] : a user-defined function which accepts an
          array of distances, and returns an array of the same shape
          containing the weights.

        Uniform weights are used by default.

    algorithm : {'auto', 'ball_tree', 'kd_tree', 'brute'}, optional
        Algorithm used to compute the nearest neighbors:

        - 'ball_tree' will use :class:`BallTree`
        - 'kd_tree' will use :class:`KDTree`
        - 'brute' will use a brute-force search.
        - 'auto' will attempt to decide the most appropriate algorithm
          based on the values passed to :meth:`fit` method.

        Note: fitting on sparse input will override the setting of
        this parameter, using brute force.

    leaf_size : int, optional (default = 30)
        Leaf size passed to BallTree or KDTree.  This can affect the
        speed of the construction and query, as well as the memory
        required to store the tree.  The optimal value depends on the
        nature of the problem.

    p : integer, optional (default = 2)
        Power parameter for the Minkowski metric. When p = 1, this is
        equivalent to using manhattan_distance (l1), and euclidean_distance
        (l2) for p = 2. For arbitrary p, minkowski_distance (l_p) is used.

    metric : string or callable, default 'minkowski'
        the distance metric to use for the tree.  The default metric is
        minkowski, and with p=2 is equivalent to the standard Euclidean
        metric. See the documentation of the DistanceMetric class for a
        list of available metrics.

    metric_params : dict, optional (default = None)
        Additional keyword arguments for the metric function.

    Examples
    --------
    >>> X = [[0], [1], [2], [3]]
    >>> y = [0, 0, 1, 1]
    >>> from sklearn.neighbors import RadiusNeighborsRegressor
    >>> neigh = RadiusNeighborsRegressor(radius=1.0)
    >>> neigh.fit(X, y) # doctest: +ELLIPSIS
    RadiusNeighborsRegressor(...)
    >>> print(neigh.predict([[1.5]]))
    [ 0.5]

    See also
    --------
    NearestNeighbors
    KNeighborsRegressor
    KNeighborsClassifier
    RadiusNeighborsClassifier

    Notes
    -----
    See :ref:`Nearest Neighbors <neighbors>` in the online documentation
    for a discussion of the choice of ``algorithm`` and ``leaf_size``.

    https://en.wikipedia.org/wiki/K-nearest_neighbor_algorithm
    """

    def __init__(self, radius=1.0, weights='uniform',
                 algorithm='auto', leaf_size=30,
                 p=2, metric='minkowski', metric_params=None, **kwargs):
        self._init_params(radius=radius,
                          algorithm=algorithm,
                          leaf_size=leaf_size,
                          p=p, metric=metric, metric_params=metric_params,
                          **kwargs)
        self.weights = _check_weights(weights)

    def predict(self, X):
        """Predict the target for the provided data

        Parameters
        ----------
        X : array-like, shape (n_query, n_features), \
                or (n_query, n_indexed) if metric == 'precomputed'
            Test samples.

        Returns
        -------
        y : array of int, shape = [n_samples] or [n_samples, n_outputs]
            Target values
        """
        X = check_array(X, accept_sparse='csr')

        neigh_dist, neigh_ind = self.radius_neighbors(X)

        weights = _get_weights(neigh_dist, self.weights)

        _y = self._y
        if _y.ndim == 1:
            _y = _y.reshape((-1, 1))

        if weights is None:
            y_pred = np.array([np.mean(_y[ind, :], axis=0)
                               for ind in neigh_ind])
        else:
            y_pred = np.array([(np.average(_y[ind, :], axis=0,
                                           weights=weights[i]))
                               for (i, ind) in enumerate(neigh_ind)])

        if self._y.ndim == 1:
            y_pred = y_pred.ravel()

        return y_pred
