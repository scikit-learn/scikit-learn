"""Base and mixin classes for nearest neighbors"""
# Authors: Jake Vanderplas <vanderplas@astro.washington.edu>
#          Fabian Pedregosa <fabian.pedregosa@inria.fr>
#          Alexandre Gramfort <alexandre.gramfort@inria.fr>
#          Sparseness support by Lars Buitinck <L.J.Buitinck@uva.nl>
#          Multi-output support by Arnaud Joly <a.joly@ulg.ac.be>
#
# License: BSD 3 clause (C) INRIA, University of Amsterdam
import warnings
from abc import ABCMeta, abstractmethod

import numpy as np
from scipy.sparse import csr_matrix, issparse

from .ball_tree import BallTree
from .kd_tree import KDTree
from ..base import BaseEstimator
from ..metrics import pairwise_distances
from ..metrics.pairwise import PAIRWISE_DISTANCE_FUNCTIONS
from ..utils import check_X_y, check_array
from ..utils.fixes import argpartition
from ..utils.validation import DataConversionWarning
from ..externals import six


VALID_METRICS = dict(ball_tree=BallTree.valid_metrics,
                     kd_tree=KDTree.valid_metrics,
                     # The following list comes from the
                     # sklearn.metrics.pairwise doc string
                     brute=(list(PAIRWISE_DISTANCE_FUNCTIONS.keys()) +
                            ['braycurtis', 'canberra', 'chebyshev',
                             'correlation', 'cosine', 'dice', 'hamming',
                             'jaccard', 'kulsinski', 'mahalanobis',
                             'matching', 'minkowski', 'rogerstanimoto',
                             'russellrao', 'seuclidean', 'sokalmichener',
                             'sokalsneath', 'sqeuclidean',
                             'yule', 'wminkowski']))


VALID_METRICS_SPARSE = dict(ball_tree=[],
                            kd_tree=[],
                            brute=PAIRWISE_DISTANCE_FUNCTIONS.keys())


class NeighborsWarning(UserWarning):
    pass


# Make sure that NeighborsWarning are displayed more than once
warnings.simplefilter("always", NeighborsWarning)


def _check_weights(weights):
    """Check to make sure weights are valid"""
    if weights in (None, 'uniform', 'distance'):
        return weights
    elif callable(weights):
        return weights
    else:
        raise ValueError("weights not recognized: should be 'uniform', "
                         "'distance', or a callable function")


def _get_weights(dist, weights):
    """Get the weights from an array of distances and a parameter ``weights``

    Parameters
    ===========
    dist: ndarray
        The input distances
    weights: {'uniform', 'distance' or a callable}
        The kind of weighting used

    Returns
    ========
    weights_arr: array of the same shape as ``dist``
        if ``weights == 'uniform'``, then returns None
    """
    if weights in (None, 'uniform'):
        return None
    elif weights == 'distance':
        with np.errstate(divide='ignore'):
            dist = 1. / dist
        return dist
    elif callable(weights):
        return weights(dist)
    else:
        raise ValueError("weights not recognized: should be 'uniform', "
                         "'distance', or a callable function")


class NeighborsBase(six.with_metaclass(ABCMeta, BaseEstimator)):
    """Base class for nearest neighbors estimators."""

    @abstractmethod
    def __init__(self):
        pass

    def _init_params(self, n_neighbors=None, radius=None,
                     algorithm='auto', leaf_size=30, metric='minkowski',
                     p=2, metric_params=None, **kwargs):
        if kwargs:
            warnings.warn("Passing additional arguments to the metric "
                          "function as **kwargs is deprecated "
                          "and will no longer be supported in 0.18. "
                          "Use metric_params instead.",
                          DeprecationWarning, stacklevel=3)
            if metric_params is None:
                metric_params = {}
            metric_params.update(kwargs)

        self.n_neighbors = n_neighbors
        self.radius = radius
        self.algorithm = algorithm
        self.leaf_size = leaf_size
        self.metric = metric
        self.metric_params = metric_params
        self.p = p

        if algorithm not in ['auto', 'brute',
                             'kd_tree', 'ball_tree']:
            raise ValueError("unrecognized algorithm: '%s'" % algorithm)

        if algorithm == 'auto':
            alg_check = 'ball_tree'
        else:
            alg_check = algorithm

        if callable(metric):
            if algorithm == 'kd_tree':
                # callable metric is only valid for brute force and ball_tree
                raise ValueError(
                    "kd_tree algorithm does not support callable metric '%s'"
                    % metric)
        elif metric not in VALID_METRICS[alg_check]:
            raise ValueError("Metric '%s' not valid for algorithm '%s'"
                             % (metric, algorithm))

        if self.metric_params is not None and 'p' in self.metric_params:
            warnings.warn("Parameter p is found in metric_params. "
                          "The corresponding parameter from __init__ "
                          "is ignored.", SyntaxWarning, stacklevel=3)
            effective_p = metric_params['p']
        else:
            effective_p = self.p

        if self.metric in ['wminkowski', 'minkowski'] and effective_p < 1:
            raise ValueError("p must be greater than one for minkowski metric")

        self._fit_X = None
        self._tree = None
        self._fit_method = None

    def _fit(self, X):
        if self.metric_params is None:
            self.effective_metric_params_ = {}
        else:
            self.effective_metric_params_ = self.metric_params.copy()

        effective_p = self.effective_metric_params_.get('p', self.p)
        if self.metric in ['wminkowski', 'minkowski']:
            self.effective_metric_params_['p'] = effective_p

        self.effective_metric_ = self.metric
        # For minkowski distance, use more efficient methods where available
        if self.metric == 'minkowski':
            p = self.effective_metric_params_.pop('p', 2)
            if p < 1:
                raise ValueError("p must be greater than one "
                                 "for minkowski metric")
            elif p == 1:
                self.effective_metric_ = 'manhattan'
            elif p == 2:
                self.effective_metric_ = 'euclidean'
            elif p == np.inf:
                self.effective_metric_ = 'chebyshev'
            else:
                self.effective_metric_params_['p'] = p

        if isinstance(X, NeighborsBase):
            self._fit_X = X._fit_X
            self._tree = X._tree
            self._fit_method = X._fit_method
            return self

        elif isinstance(X, BallTree):
            self._fit_X = X.data
            self._tree = X
            self._fit_method = 'ball_tree'
            return self

        elif isinstance(X, KDTree):
            self._fit_X = X.data
            self._tree = X
            self._fit_method = 'kd_tree'
            return self

        X = check_array(X, accept_sparse='csr')

        n_samples = X.shape[0]
        if n_samples == 0:
            raise ValueError("n_samples must be greater than 0")

        if issparse(X):
            if self.algorithm not in ('auto', 'brute'):
                warnings.warn("cannot use tree with sparse input: "
                              "using brute force")
            if self.effective_metric_ not in VALID_METRICS_SPARSE['brute']:
                raise ValueError("metric '%s' not valid for sparse input"
                                 % self.effective_metric_)
            self._fit_X = X.copy()
            self._tree = None
            self._fit_method = 'brute'
            return self

        self._fit_method = self.algorithm
        self._fit_X = X

        if self._fit_method == 'auto':
            # A tree approach is better for small number of neighbors,
            # and KDTree is generally faster when available
            if (self.n_neighbors is None
                    or self.n_neighbors < self._fit_X.shape[0] // 2):
                if self.effective_metric_ in VALID_METRICS['kd_tree']:
                    self._fit_method = 'kd_tree'
                else:
                    self._fit_method = 'ball_tree'
            else:
                self._fit_method = 'brute'

        if self._fit_method == 'ball_tree':
            self._tree = BallTree(X, self.leaf_size,
                                  metric=self.effective_metric_,
                                  **self.effective_metric_params_)
        elif self._fit_method == 'kd_tree':
            self._tree = KDTree(X, self.leaf_size,
                                metric=self.effective_metric_,
                                **self.effective_metric_params_)
        elif self._fit_method == 'brute':
            self._tree = None
        else:
            raise ValueError("algorithm = '%s' not recognized"
                             % self.algorithm)
        return self


class KNeighborsMixin(object):
    """Mixin for k-neighbors searches"""

    def kneighbors(self, X, n_neighbors=None, return_distance=True):
        """Finds the K-neighbors of a point.

        Returns distance

        Parameters
        ----------
        X : array-like, last dimension same as that of fit data
            The new point.

        n_neighbors : int
            Number of neighbors to get (default is the value
            passed to the constructor).

        return_distance : boolean, optional. Defaults to True.
            If False, distances will not be returned

        Returns
        -------
        dist : array
            Array representing the lengths to point, only present if
            return_distance=True

        ind : array
            Indices of the nearest points in the population matrix.

        Examples
        --------
        In the following example, we construct a NeighborsClassifier
        class from an array representing our data set and ask who's
        the closest point to [1,1,1]

        >>> samples = [[0., 0., 0.], [0., .5, 0.], [1., 1., .5]]
        >>> from sklearn.neighbors import NearestNeighbors
        >>> neigh = NearestNeighbors(n_neighbors=1)
        >>> neigh.fit(samples) # doctest: +ELLIPSIS
        NearestNeighbors(algorithm='auto', leaf_size=30, ...)
        >>> print(neigh.kneighbors([1., 1., 1.])) # doctest: +ELLIPSIS
        (array([[ 0.5]]), array([[2]]...))

        As you can see, it returns [[0.5]], and [[2]], which means that the
        element is at distance 0.5 and is the third element of samples
        (indexes start at 0). You can also query for multiple points:

        >>> X = [[0., 1., 0.], [1., 0., 1.]]
        >>> neigh.kneighbors(X, return_distance=False) # doctest: +ELLIPSIS
        array([[1],
               [2]]...)

        """
        if self._fit_method is None:
            raise ValueError("must fit neighbors before querying")

        X = check_array(X, accept_sparse='csr')

        if n_neighbors is None:
            n_neighbors = self.n_neighbors

        if self._fit_method == 'brute':
            # for efficiency, use squared euclidean distances
            if self.effective_metric_ == 'euclidean':
                dist = pairwise_distances(X, self._fit_X, 'euclidean',
                                          squared=True)
            else:
                dist = pairwise_distances(X, self._fit_X,
                                          self.effective_metric_,
                                          **self.effective_metric_params_)

            neigh_ind = argpartition(dist, n_neighbors - 1, axis=1)
            neigh_ind = neigh_ind[:, :n_neighbors]
            # argpartition doesn't guarantee sorted order, so we sort again
            j = np.arange(neigh_ind.shape[0])[:, None]
            neigh_ind = neigh_ind[j, np.argsort(dist[j, neigh_ind])]
            if return_distance:
                if self.effective_metric_ == 'euclidean':
                    return np.sqrt(dist[j, neigh_ind]), neigh_ind
                else:
                    return dist[j, neigh_ind], neigh_ind
            else:
                return neigh_ind
        elif self._fit_method in ['ball_tree', 'kd_tree']:
            result = self._tree.query(X, n_neighbors,
                                      return_distance=return_distance)
            return result
        else:
            raise ValueError("internal: _fit_method not recognized")

    def kneighbors_graph(self, X, n_neighbors=None,
                         mode='connectivity'):
        """Computes the (weighted) graph of k-Neighbors for points in X

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Sample data

        n_neighbors : int
            Number of neighbors for each sample.
            (default is value passed to the constructor).

        mode : {'connectivity', 'distance'}, optional
            Type of returned matrix: 'connectivity' will return the
            connectivity matrix with ones and zeros, in 'distance' the
            edges are Euclidean distance between points.

        Returns
        -------
        A : sparse matrix in CSR format, shape = [n_samples, n_samples_fit]
            n_samples_fit is the number of samples in the fitted data
            A[i, j] is assigned the weight of edge that connects i to j.

        Examples
        --------
        >>> X = [[0], [3], [1]]
        >>> from sklearn.neighbors import NearestNeighbors
        >>> neigh = NearestNeighbors(n_neighbors=2)
        >>> neigh.fit(X) # doctest: +ELLIPSIS
        NearestNeighbors(algorithm='auto', leaf_size=30, ...)
        >>> A = neigh.kneighbors_graph(X)
        >>> A.toarray()
        array([[ 1.,  0.,  1.],
               [ 0.,  1.,  1.],
               [ 1.,  0.,  1.]])

        See also
        --------
        NearestNeighbors.radius_neighbors_graph
        """
        X = check_array(X, accept_sparse=['csr', 'csc', 'coo'])

        if n_neighbors is None:
            n_neighbors = self.n_neighbors

        n_samples1 = X.shape[0]
        n_samples2 = self._fit_X.shape[0]
        n_nonzero = n_samples1 * n_neighbors
        A_indptr = np.arange(0, n_nonzero + 1, n_neighbors)

        # construct CSR matrix representation of the k-NN graph
        if mode == 'connectivity':
            A_data = np.ones((n_samples1, n_neighbors))
            A_ind = self.kneighbors(X, n_neighbors, return_distance=False)

        elif mode == 'distance':
            data, ind = self.kneighbors(X, n_neighbors + 1,
                                        return_distance=True)
            A_data, A_ind = data[:, 1:], ind[:, 1:]

        else:
            raise ValueError(
                'Unsupported mode, must be one of "connectivity" '
                'or "distance" but got "%s" instead' % mode)

        return csr_matrix((A_data.ravel(), A_ind.ravel(), A_indptr),
                          shape=(n_samples1, n_samples2))


class RadiusNeighborsMixin(object):
    """Mixin for radius-based neighbors searches"""

    def radius_neighbors(self, X, radius=None, return_distance=True):
        """Finds the neighbors within a given radius of a point or points.

        Returns indices of and distances to the neighbors of each point.

        Parameters
        ----------
        X : array-like, last dimension same as that of fit data
            The new point or points

        radius : float
            Limiting distance of neighbors to return.
            (default is the value passed to the constructor).

        return_distance : boolean, optional. Defaults to True.
            If False, distances will not be returned

        Returns
        -------
        dist : array
            Array representing the euclidean distances to each point,
            only present if return_distance=True.

        ind : array
            Indices of the nearest points in the population matrix.

        Examples
        --------
        In the following example, we construct a NeighborsClassifier
        class from an array representing our data set and ask who's
        the closest point to [1,1,1]

        >>> samples = [[0., 0., 0.], [0., .5, 0.], [1., 1., .5]]
        >>> from sklearn.neighbors import NearestNeighbors
        >>> neigh = NearestNeighbors(radius=1.6)
        >>> neigh.fit(samples) # doctest: +ELLIPSIS
        NearestNeighbors(algorithm='auto', leaf_size=30, ...)
        >>> print(neigh.radius_neighbors([1., 1., 1.])) # doctest: +ELLIPSIS
        (array([[ 1.5,  0.5]]...), array([[1, 2]]...)

        The first array returned contains the distances to all points which
        are closer than 1.6, while the second array returned contains their
        indices.  In general, multiple points can be queried at the same time.

        Notes
        -----
        Because the number of neighbors of each point is not necessarily
        equal, the results for multiple query points cannot be fit in a
        standard data array.
        For efficiency, `radius_neighbors` returns arrays of objects, where
        each object is a 1D array of indices or distances.
        """

        if self._fit_method is None:
            raise ValueError("must fit neighbors before querying")

        X = check_array(X, accept_sparse='csr')

        if radius is None:
            radius = self.radius

        if self._fit_method == 'brute':
            # for efficiency, use squared euclidean distances
            if self.effective_metric_ == 'euclidean':
                dist = pairwise_distances(X, self._fit_X, 'euclidean',
                                          squared=True)
                radius *= radius
            else:
                dist = pairwise_distances(X, self._fit_X,
                                          self.effective_metric_,
                                          **self.effective_metric_params_)
            neigh_ind = [np.where(d < radius)[0] for d in dist]

            # if there are the same number of neighbors for each point,
            # we can do a normal array.  Otherwise, we return an object
            # array with elements that are numpy arrays
            try:
                neigh_ind = np.asarray(neigh_ind, dtype=int)
                dtype_F = float
            except ValueError:
                neigh_ind = np.asarray(neigh_ind, dtype='object')
                dtype_F = object

            if return_distance:
                if self.effective_metric_ == 'euclidean':
                    dist = np.array([np.sqrt(d[neigh_ind[i]])
                                     for i, d in enumerate(dist)],
                                    dtype=dtype_F)
                else:
                    dist = np.array([d[neigh_ind[i]]
                                     for i, d in enumerate(dist)],
                                    dtype=dtype_F)
                return dist, neigh_ind
            else:
                return neigh_ind
        elif self._fit_method in ['ball_tree', 'kd_tree']:
            results = self._tree.query_radius(X, radius,
                                              return_distance=return_distance)
            if return_distance:
                ind, dist = results
                return dist, ind
            else:
                return results
        else:
            raise ValueError("internal: _fit_method not recognized")

    def radius_neighbors_graph(self, X, radius=None, mode='connectivity'):
        """Computes the (weighted) graph of Neighbors for points in X

        Neighborhoods are restricted the points at a distance lower than
        radius.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Sample data

        radius : float
            Radius of neighborhoods.
            (default is the value passed to the constructor).

        mode : {'connectivity', 'distance'}, optional
            Type of returned matrix: 'connectivity' will return the
            connectivity matrix with ones and zeros, in 'distance' the
            edges are Euclidean distance between points.

        Returns
        -------
        A : sparse matrix in CSR format, shape = [n_samples, n_samples]
            A[i, j] is assigned the weight of edge that connects i to j.

        Examples
        --------
        >>> X = [[0], [3], [1]]
        >>> from sklearn.neighbors import NearestNeighbors
        >>> neigh = NearestNeighbors(radius=1.5)
        >>> neigh.fit(X) # doctest: +ELLIPSIS
        NearestNeighbors(algorithm='auto', leaf_size=30, ...)
        >>> A = neigh.radius_neighbors_graph(X)
        >>> A.toarray()
        array([[ 1.,  0.,  1.],
               [ 0.,  1.,  0.],
               [ 1.,  0.,  1.]])

        See also
        --------
        kneighbors_graph
        """
        X = check_array(X, accept_sparse=['csr', 'csc', 'coo'])

        if radius is None:
            radius = self.radius

        n_samples1 = X.shape[0]
        n_samples2 = self._fit_X.shape[0]

        # construct CSR matrix representation of the NN graph
        if mode == 'connectivity':
            A_ind = self.radius_neighbors(X, radius,
                                          return_distance=False)
            A_data = None
        elif mode == 'distance':
            dist, A_ind = self.radius_neighbors(X, radius,
                                                return_distance=True)
            A_data = np.concatenate(list(dist))
        else:
            raise ValueError(
                'Unsupported mode, must be one of "connectivity", '
                'or "distance" but got %s instead' % mode)

        n_neighbors = np.array([len(a) for a in A_ind])
        n_nonzero = np.sum(n_neighbors)
        if A_data is None:
            A_data = np.ones(n_nonzero)
        A_ind = np.concatenate(list(A_ind))
        A_indptr = np.concatenate((np.zeros(1, dtype=int),
                                   np.cumsum(n_neighbors)))

        return csr_matrix((A_data, A_ind, A_indptr),
                          shape=(n_samples1, n_samples2))


class SupervisedFloatMixin(object):
    def fit(self, X, y):
        """Fit the model using X as training data and y as target values

        Parameters
        ----------
        X : {array-like, sparse matrix, BallTree, KDTree}
            Training data. If array or matrix, shape = [n_samples, n_features]

        y : {array-like, sparse matrix}
            Target values, array of float values, shape = [n_samples]
             or [n_samples, n_outputs]
        """
        if not isinstance(X, (KDTree, BallTree)):
            X, y = check_X_y(X, y, "csr", multi_output=True)
        self._y = y
        return self._fit(X)


class SupervisedIntegerMixin(object):
    def fit(self, X, y):
        """Fit the model using X as training data and y as target values

        Parameters
        ----------
        X : {array-like, sparse matrix, BallTree, KDTree}
            Training data. If array or matrix, shape = [n_samples, n_features]

        y : {array-like, sparse matrix}
            Target values of shape = [n_samples] or [n_samples, n_outputs]

        """
        if not isinstance(X, (KDTree, BallTree)):
            X, y = check_X_y(X, y, "csr", multi_output=True)

        if y.ndim == 1 or y.ndim == 2 and y.shape[1] == 1:
            if y.ndim != 1:
                warnings.warn("A column-vector y was passed when a 1d array "
                              "was expected. Please change the shape of y to "
                              "(n_samples, ), for example using ravel().",
                              DataConversionWarning, stacklevel=2)

            self.outputs_2d_ = False
            y = y.reshape((-1, 1))
        else:
            self.outputs_2d_ = True

        self.classes_ = []
        self._y = np.empty(y.shape, dtype=np.int)
        for k in range(self._y.shape[1]):
            classes, self._y[:, k] = np.unique(y[:, k], return_inverse=True)
            self.classes_.append(classes)

        if not self.outputs_2d_:
            self.classes_ = self.classes_[0]
            self._y = self._y.ravel()

        return self._fit(X)


class UnsupervisedMixin(object):
    def fit(self, X, y=None):
        """Fit the model using X as training data

        Parameters
        ----------
        X : {array-like, sparse matrix, BallTree, KDTree}
            Training data. If array or matrix, shape = [n_samples, n_features]
        """
        return self._fit(X)
