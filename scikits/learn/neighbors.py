"""
k-Nearest Neighbor Algorithm.

Uses BallTree algorithm, which is an efficient way to perform fast
neighbor searches in high dimensionality.
"""
import numpy as np
from scipy import stats
from scipy import linalg

from .base import BaseEstimator, ClassifierMixin, RegressorMixin
from .ball_tree import BallTree


class Neighbors(BaseEstimator, ClassifierMixin):
    """Classifier implementing k-Nearest Neighbor Algorithm.

    Parameters
    ----------
    data : array-like, shape (n, k)
        The data points to be indexed. This array is not copied, and so
        modifying this data will result in bogus results.
    labels : array
        An array representing labels for the data (only arrays of
        integers are supported).
    n_neighbors : int
        default number of neighbors.
    window_size : int
        Window size passed to BallTree

    Examples
    --------
    >>> samples = [[0.,0.,1.], [1.,0.,0.], [2.,2.,2.], [2.,5.,4.]]
    >>> labels = [0,0,1,1]
    >>> from scikits.learn.neighbors import Neighbors
    >>> neigh = Neighbors(n_neighbors=3)
    >>> neigh.fit(samples, labels)
    Neighbors(n_neighbors=3, window_size=1)
    >>> print neigh.predict([[0,0,0]])
    [ 0.]

    Notes
    -----
    http://en.wikipedia.org/wiki/K-nearest_neighbor_algorithm
    """

    def __init__(self, n_neighbors=5, window_size=1):
        """Internally uses the ball tree datastructure and algorithm for fast
        neighbors lookups on high dimensional datasets.
        """
        self.n_neighbors = n_neighbors
        self.window_size = window_size

    def fit(self, X, Y=()):
        # we need Y to be an integer, because after we'll use it an index
        self.Y = np.asanyarray(Y, dtype=np.int)
        self.ball_tree = BallTree(X, self.window_size)
        return self

    def kneighbors(self, data, n_neighbors=None):
        """Finds the K-neighbors of a point.

        Parameters
        ----------
        point : array-like
            The new point.
        n_neighbors : int
            Number of neighbors to get (default is the value
            passed to the constructor).

        Returns
        -------
        dist : array
            Array representing the lengths to point.
        ind : array
            Array representing the indices of the nearest points in the
            population matrix.

        Examples
        --------
        In the following example, we construnct a Neighbors class from an
        array representing our data set and ask who's the closest point to
        [1,1,1]

        >>> samples = [[0., 0., 0.], [0., .5, 0.], [1., 1., .5]]
        >>> labels = [0, 0, 1]
        >>> from scikits.learn.neighbors import Neighbors
        >>> neigh = Neighbors(n_neighbors=1)
        >>> neigh.fit(samples, labels)
        Neighbors(n_neighbors=1, window_size=1)
        >>> print neigh.kneighbors([1., 1., 1.])
        (array(0.5), array(2))

        As you can see, it returns [0.5], and [2], which means that the
        element is at distance 0.5 and is the third element of samples
        (indexes start at 0). You can also query for multiple points:

        >>> print neigh.kneighbors([[0., 1., 0.], [1., 0., 1.]])
        (array([ 0.5       ,  1.11803399]), array([1, 2]))

        """
        if n_neighbors is None:
            n_neighbors = self.n_neighbors
        return self.ball_tree.query(data, k=n_neighbors)

    def predict(self, T, n_neighbors=None):
        """Predict the class labels for the provided data.

        Parameters
        ----------
        test: array
            A 2-D array representing the test point.
        n_neighbors : int
            Number of neighbors to get (default is the value
            passed to the constructor).

        Returns
        -------
        labels: array
            List of class labels (one for each data sample).

        Examples
        --------
        >>> samples = [[0., 0., 0.], [0., .5, 0.], [1., 1., .5]]
        >>> labels = [0, 0, 1]
        >>> from scikits.learn.neighbors import Neighbors
        >>> neigh = Neighbors(n_neighbors=1)
        >>> neigh.fit(samples, labels)
        Neighbors(n_neighbors=1, window_size=1)
        >>> print neigh.predict([.2, .1, .2])
        0
        >>> print neigh.predict([[0., -1., 0.], [3., 2., 0.]])
        [0 1]
        """
        T = np.asanyarray(T)
        if n_neighbors is None:
            n_neighbors = self.n_neighbors
        return _predict_from_BallTree(self.ball_tree, self.Y, T, n_neighbors)


def _predict_from_BallTree(ball_tree, Y, test, n_neighbors):
    """Predict target from BallTree object containing the data points.

    This is a helper method, not meant to be used directly. It will
    not check that input is of the correct type.
    """
    Y_ = Y[ball_tree.query(test, k=n_neighbors, return_distance=False)]
    if n_neighbors == 1:
        return Y_
    return (stats.mode(Y_, axis=1)[0]).ravel()

###############################################################################
# Neighbors Barycenter class for regression problems

class NeighborsBarycenter(BaseEstimator, RegressorMixin):
    """Regression based on k-Nearest Neighbor Algorithm.

    The target is predicted by local interpolation of the targets
    associated of the k-Nearest Neighbors in the training set.
    The interpolation weights correspond to barycenter weights.

    Parameters
    ----------
    X : array-like, shape (n_samples, n_features)
        The data points to be indexed. This array is not copied, and so
        modifying this data will result in bogus results.
    y : array
        An array representing labels for the data (only arrays of
        integers are supported).
    n_neighbors : int
        default number of neighbors.
    window_size : int
        Window size passed to BallTree

    Examples
    --------
    >>> X = [[0], [1], [2], [3]]
    >>> y = [0, 0, 1, 1]
    >>> from scikits.learn.neighbors import NeighborsBarycenter
    >>> neigh = NeighborsBarycenter(n_neighbors=2)
    >>> neigh.fit(X, y)
    NeighborsBarycenter(n_neighbors=2, window_size=1)
    >>> print neigh.predict([[1.5]])
    [ 0.5]

    Notes
    -----
    http://en.wikipedia.org/wiki/K-nearest_neighbor_algorithm
    """

    def __init__(self, n_neighbors=5, window_size=1):
        """Internally uses the ball tree datastructure and algorithm for fast
        neighbors lookups on high dimensional datasets.
        """
        self.n_neighbors = n_neighbors
        self.window_size = window_size

    def fit(self, X, y, copy=True):
        self._y = np.array(y, copy=copy)
        self.ball_tree = BallTree(X, self.window_size)
        return self

    def predict(self, T, n_neighbors=None):
        """Predict the target for the provided data.

        Parameters
        ----------
        T : array
            A 2-D array representing the test data.
        n_neighbors : int
            Number of neighbors to get (default is the value
            passed to the constructor).

        Returns
        -------
        y: array
            List of target values (one for each data sample).

        Examples
        --------
        >>> X = [[0], [1], [2]]
        >>> y = [0, 0, 1]
        >>> from scikits.learn.neighbors import NeighborsBarycenter
        >>> neigh = NeighborsBarycenter(n_neighbors=2)
        >>> neigh.fit(X, y)
        NeighborsBarycenter(n_neighbors=2, window_size=1)
        >>> print neigh.predict([[.5], [1.5]])
        [ 0.   0.5]
        """
        T = np.asanyarray(T)
        if T.ndim == 1:
            T = T[:,None]
        if n_neighbors is None:
            n_neighbors = self.n_neighbors
        A = kneighbors_graph(T, n_neighbors=n_neighbors, weight="barycenter",
                                  ball_tree=self.ball_tree).tocsr()
        return A * self._y

###############################################################################
# Utils k-NN based Functions

def barycenter_weights(x, X_neighbors, tol=1e-3):
    """Computes barycenter weights

    We estimate the weights to assign to each point in X_neighbors
    to recover the point x. The barycenter weights sum to 1.
    If x do not belong to the span of X_neighbors, it's the
    projection of x onto the span that is recovered.

    Parameters
    ----------
    x : array
        a 1D array

    X_neighbors : array
        a 2D array containing samples

    tol : float
        tolerance

    Returns
    -------
    array of barycenter weights that sum to 1

    Examples
    --------
    >>> X_neighbors = [[0], [2]]
    >>> x = [0.5]
    >>> from scikits.learn.neighbors import barycenter_weights
    >>> print barycenter_weights(x, X_neighbors)
    [ 0.74968789  0.25031211]
    """
    x = np.asanyarray(x)
    X_neighbors = np.asanyarray(X_neighbors)
    if x.ndim == 1:
        x = x[None,:]
    if X_neighbors.ndim == 1:
        X_neighbors = X_neighbors[:,None]
    z = x - X_neighbors
    gram = np.dot(z, z.T)
    # Add constant on diagonal to avoid singular matrices
    diag_stride = gram.shape[0] + 1
    gram.flat[::diag_stride] += tol * np.trace(gram)
    w = linalg.solve(gram, np.ones(len(X_neighbors)))
    w /= np.sum(w)
    return w


def kneighbors_graph(X, n_neighbors, weight=None, ball_tree=None,
                     window_size=1):
    """Computes the (weighted) graph of k-Neighbors

    Parameters
    ----------
    X : array-like, shape = [n_samples, n_features]
        Coordinates of samples. One sample per row.

    n_neighbors : int
        Number of neighbors for each sample.

    weight : None (default)
        Weights to apply on graph edges. If weight is None
        then no weighting is applied (1 for each edge).
        If weight equals "distance" the edge weight is the
        euclidian distance. If weight equals "barycenter"
        the weights are barycenter weights estimated by
        solving a linear system for each point.

    ball_tree : None or instance of precomputed BallTree

    window_size : int
        Window size pass to the BallTree

    Returns
    -------
    A : sparse matrix, shape = [n_samples, n_samples]
        A is returned as LInked List Sparse matrix
        A[i,j] = weight of edge that connects i to j

    Examples
    --------
    >>> X = [[0], [2], [1]]
    >>> A = kneighbors_graph(X, n_neighbors=2, weight=None)
    >>> print A.todense()
    [[ 1.  0.  1.]
     [ 0.  1.  1.]
     [ 0.  1.  1.]]
    """
    from scipy import sparse
    X = np.asanyarray(X)
    n_samples = X.shape[0]
    if ball_tree is None:
        ball_tree = BallTree(X, window_size)
    A = sparse.lil_matrix((n_samples, ball_tree.size))
    dist, ind = ball_tree.query(X, k=n_neighbors)
    if weight is None:
        for i, li in enumerate(ind):
            if n_neighbors > 1:
                A[i, list(li)] = np.ones(n_neighbors)
            else:
                A[i, li] = 1.0
    elif weight is "distance":
        for i, li in enumerate(ind):
            if n_neighbors > 1:
                A[i, list(li)] = dist[i, :]
            else:
                A[i, li] = dist[i, 0]
    elif weight is "barycenter":
        # XXX : the next loop could be done in parallel
        # by parallelizing groups of indices
        for i, li in enumerate(ind):
            if n_neighbors > 1:
                X_i = ball_tree.data[li]
                A[i, list(li)] = barycenter_weights(X[i], X_i)
            else:
                A[i, li] = 1.0
    else:
        raise ValueError("Unknown weight type")
    return A
