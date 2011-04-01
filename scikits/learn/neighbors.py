"""Nearest Neighbor related algorithms"""

# Author: Fabian Pedregosa <fabian.pedregosa@inria.fr>
#         Alexandre Gramfort <alexandre.gramfort@inria.fr>
#
# License: BSD, (C) INRIA

import numpy as np

from .base import BaseEstimator, ClassifierMixin, RegressorMixin
from .ball_tree import BallTree
from .metrics import euclidean_distances


class NeighborsClassifier(BaseEstimator, ClassifierMixin):
    """Classifier implementing k-Nearest Neighbor Algorithm.

    Parameters
    ----------
    n_neighbors : int, optional
        Default number of neighbors. Defaults to 5.

    window_size : int, optional
        Window size passed to BallTree

    algorithm : {'auto', 'ball_tree', 'brute', 'brute_inplace'}, optional
        Algorithm used to compute the nearest neighbors. 'ball_tree'
        will construct a BallTree, 'brute' and 'brute_inplace' will
        perform brute-force search.'auto' will guess the most
        appropriate based on current dataset.

    Examples
    --------
    >>> samples = [[0, 0, 1], [1, 0, 0]]
    >>> labels = [0, 1]
    >>> from scikits.learn.neighbors import NeighborsClassifier
    >>> neigh = NeighborsClassifier(n_neighbors=1)
    >>> neigh.fit(samples, labels)
    NeighborsClassifier(n_neighbors=1, window_size=1, algorithm='auto')
    >>> print neigh.predict([[0,0,0]])
    [1]

    See also
    --------
    BallTree

    References
    ----------
    http://en.wikipedia.org/wiki/K-nearest_neighbor_algorithm
    """

    def __init__(self, n_neighbors=5, algorithm='auto', window_size=1):
        self.n_neighbors = n_neighbors
        self.window_size = window_size
        self.algorithm = algorithm

    def fit(self, X, y, **params):
        """Fit the model using X, y as training data

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Training data.

        y : array-like, shape = [n_samples]
            Target values, array of integer values.

        params : list of keyword, optional
            Overwrite keywords from __init__
        """
        X = np.asanyarray(X)
        if y is None:
            raise ValueError("y must not be None")
        self._y = np.asanyarray(y)
        self._set_params(**params)

        if self.algorithm == 'ball_tree' or \
           (self.algorithm == 'auto' and X.shape[1] < 20):
            self.ball_tree = BallTree(X, self.window_size)
        else:
            self.ball_tree = None
            self._fit_X = X
        return self

    def kneighbors(self, X, return_distance=True, **params):
        """Finds the K-neighbors of a point.

        Returns distance

        Parameters
        ----------
        point : array-like
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
        In the following example, we construnct a NeighborsClassifier
        class from an array representing our data set and ask who's
        the closest point to [1,1,1]

        >>> samples = [[0., 0., 0.], [0., .5, 0.], [1., 1., .5]]
        >>> labels = [0, 0, 1]
        >>> from scikits.learn.neighbors import NeighborsClassifier
        >>> neigh = NeighborsClassifier(n_neighbors=1)
        >>> neigh.fit(samples, labels)
        NeighborsClassifier(n_neighbors=1, window_size=1, algorithm='auto')
        >>> print neigh.kneighbors([1., 1., 1.]) # doctest: +ELLIPSIS
        (array([[ 0.5]]), array([[2]]...))

        As you can see, it returns [[0.5]], and [[2]], which means that the
        element is at distance 0.5 and is the third element of samples
        (indexes start at 0). You can also query for multiple points:

        >>> X = [[0., 1., 0.], [1., 0., 1.]]
        >>> neigh.kneighbors(X, return_distance=False) # doctest: +ELLIPSIS
        array([[1],
               [2]]...)

        """
        self._set_params(**params)
        X = np.atleast_2d(X)
        if self.ball_tree is None:
            if self.algorithm == 'brute_inplace' and not return_distance:
                return knn_brute(self._fit_X, X, self.n_neighbors)
            else:
                dist = euclidean_distances(X, self._fit_X, squared=True)
                # XXX: should be implemented with a partial sort
                neigh_ind = dist.argsort(axis=1)[:, :self.n_neighbors]
            if not return_distance:
                return neigh_ind
            else:
                return dist.T[neigh_ind], neigh_ind
        else:
            return self.ball_tree.query(X, self.n_neighbors,
                                        return_distance=return_distance)

    def predict(self, X, **params):
        """Predict the class labels for the provided data

        Parameters
        ----------
        X: array
            A 2-D array representing the test point.

        n_neighbors : int
            Number of neighbors to get (default is the value
            passed to the constructor).

        Returns
        -------
        labels: array
            List of class labels (one for each data sample).
        """
        X = np.atleast_2d(X)
        self._set_params(**params)

        # get neighbors
        neigh_ind = self.kneighbors(X, return_distance=False)

        # compute the most popular label
        pred_labels = self._y[neigh_ind]
        from scipy import stats
        mode, _ = stats.mode(pred_labels, axis=1)
        return mode.flatten().astype(np.int)


###############################################################################
# NeighborsRegressor class for regression problems

class NeighborsRegressor(NeighborsClassifier, RegressorMixin):
    """Regression based on k-Nearest Neighbor Algorithm

    The target is predicted by local interpolation of the targets
    associated of the k-Nearest Neighbors in the training set.

    Different modes for estimating the result can be set via parameter
    mode. 'barycenter' will apply the weights that best reconstruct
    the point from its neighbors while 'mean' will apply constant
    weights to each point.

    Parameters
    ----------
    n_neighbors : int, optional
        Default number of neighbors. Defaults to 5.

    window_size : int, optional
        Window size passed to BallTree

    mode : {'mean', 'barycenter'}, optional
        Weights to apply to labels.

    algorithm : {'auto', 'ball_tree', 'brute', 'brute_inplace'}, optional
        Algorithm used to compute the nearest neighbors. 'ball_tree'
        will construct a BallTree, 'brute' and 'brute_inplace' will
        perform brute-force search.'auto' will guess the most
        appropriate based on current dataset.

    Examples
    --------
    >>> X = [[0], [1], [2], [3]]
    >>> y = [0, 0, 1, 1]
    >>> from scikits.learn.neighbors import NeighborsRegressor
    >>> neigh = NeighborsRegressor(n_neighbors=2)
    >>> neigh.fit(X, y)
    NeighborsRegressor(n_neighbors=2, window_size=1, mode='mean',
              algorithm='auto')
    >>> print neigh.predict([[1.5]])
    [ 0.5]

    Notes
    -----
    http://en.wikipedia.org/wiki/K-nearest_neighbor_algorithm
    """

    def __init__(self, n_neighbors=5, mode='mean', algorithm='auto',
                 window_size=1):
        self.n_neighbors = n_neighbors
        self.window_size = window_size
        self.mode = mode
        self.algorithm = algorithm

    def predict(self, X, **params):
        """Predict the target for the provided data

        Parameters
        ----------
        X : array
            A 2-D array representing the test data.

        n_neighbors : int, optional
            Number of neighbors to get (default is the value
            passed to the constructor).

        Returns
        -------
        y: array
            List of target values (one for each data sample).
        """
        X = np.atleast_2d(np.asanyarray(X))
        self._set_params(**params)

        # compute nearest neighbors
        neigh_ind = self.kneighbors(X, return_distance=False)
        if self.ball_tree is None:
            neigh = self._fit_X[neigh_ind]
        else:
            neigh = self.ball_tree.data[neigh_ind]

        # compute interpolation on y
        if self.mode == 'barycenter':
            W = barycenter_weights(X, neigh)
            return (W * self._y[neigh_ind]).sum(axis=1)

        elif self.mode == 'mean':
            return np.mean(self._y[neigh_ind], axis=1)

        else:
            raise ValueError(
                'Unsupported mode, must be one of "barycenter" or '
                '"mean" but got %s instead' % self.mode)


###############################################################################
# Utils k-NN based Functions

def barycenter_weights(X, Z, cond=None):
    """Compute barycenter weights of X from Y along the first axis

    We estimate the weights to assign to each point in Y[i] to recover
    the point X[i]. The barycenter weights sum to 1.

    Parameters
    ----------
    X : array-like, shape (n_samples, n_dim)

    Z : array-like, shape (n_samples, n_neighbors, n_dim)

    cond: float, optional
        Cutoff for small singular values; used to determine effective
        rank of Z[i]. Singular values smaller than ``rcond *
        largest_singular_value`` are considered zero.

    Returns
    -------
    B : array-like, shape (n_samples, n_neighbors)

    Notes
    -----
    See developers note for more information.
    """
    from scipy import linalg
    X, Z = map(np.asanyarray, (X, Z))
    n_samples, n_neighbors = X.shape[0], Z.shape[1]
    if X.dtype.kind == 'i':
        X = X.astype(np.float)
    B = np.empty((n_samples, n_neighbors), dtype=X.dtype)
    v = np.ones(n_neighbors, dtype=X.dtype)
    rank_update, = linalg.get_blas_funcs(('ger',), (X,))

    # constrained least squares
    v[0] -= np.sqrt(n_neighbors)
    B[:, 0] = 1. / np.sqrt(n_neighbors)
    if n_neighbors <= 1:
        return B
    alpha = - 1. / (n_neighbors - np.sqrt(n_neighbors))
    for i, A in enumerate(Z.transpose(0, 2, 1)):
        C = rank_update(alpha, np.dot(A, v), v, a=A)
        B[i, 1:] = linalg.lstsq(
            C[:, 1:], X[i] - C[:, 0] / np.sqrt(n_neighbors), cond=cond,
            overwrite_a=True, overwrite_b=True)[0].ravel()
        B[i] = rank_update(alpha, v, np.dot(v.T, B[i]), a=B[i])

    return B


def kneighbors_graph(X, n_neighbors, mode='connectivity'):
    """Computes the (weighted) graph of k-Neighbors for points in X

    Parameters
    ----------
    X : array-like, shape = [n_samples, n_features]
        Coordinates of samples. One sample per row.

    n_neighbors : int
        Number of neighbors for each sample.

    mode : {'connectivity', 'distance', 'barycenter'}
        Type of returned matrix: 'connectivity' will return the
        connectivity matrix with ones and zeros, in 'distance' the
        edges are euclidian distance between points. In 'barycenter'
        they are barycenter weights estimated by solving a linear
        system for each point.

    Returns
    -------
    A : CSR sparse matrix, shape = [n_samples, n_samples]
        A[i,j] is assigned the weight of edge that connects i to j.

    Examples
    --------
    >>> X = [[0], [3], [1]]
    >>> from scikits.learn.neighbors import kneighbors_graph
    >>> A = kneighbors_graph(X, 2)
    >>> A.todense()
    matrix([[ 1.,  0.,  1.],
            [ 0.,  1.,  1.],
            [ 1.,  0.,  1.]])
    """
    from scipy import sparse
    X = np.asanyarray(X)
    n_samples = X.shape[0]
    ball_tree = BallTree(X)
    n_nonzero = n_neighbors * n_samples
    A_indptr = np.arange(0, n_nonzero + 1, n_neighbors)

    # construct CSR matrix representation of the k-NN graph
    if mode is 'connectivity':
        A_data = np.ones((n_samples, n_neighbors))
        A_ind = ball_tree.query(
            X, k=n_neighbors, return_distance=False)

    elif mode is 'distance':
        data, ind = ball_tree.query(X, k=n_neighbors + 1)
        A_data, A_ind = data[:, 1:], ind[:, 1:]

    elif mode is 'barycenter':
        ind = ball_tree.query(
            X, k=n_neighbors + 1, return_distance=False)
        A_ind = ind[:, 1:]
        A_data = barycenter_weights(X, X[A_ind])

    else:
        raise ValueError(
            'Unsupported mode, must be one of "connectivity", '
            '"distance" or "barycenter" but got %s instead' % mode)

    A = sparse.csr_matrix(
        (A_data.reshape(-1), A_ind.reshape(-1), A_indptr),
        shape=(n_samples, n_samples))

    return A
