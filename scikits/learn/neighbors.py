"""Nearest Neighbor related algorithms"""

# Author: Fabian Pedregosa <fabian.pedregosa@inria.fr>
#         Alexandre Gramfort <alexandre.gramfort@inria.fr>
#         Sparseness support by Lars Buitinck <L.J.Buitinck@uva.nl>
#
# License: BSD, (C) INRIA, University of Amsterdam

import numpy as np
from scipy import linalg
from scipy.sparse import csr_matrix, issparse

from .base import BaseEstimator, ClassifierMixin, RegressorMixin
from .ball_tree import BallTree
from .metrics import euclidean_distances
from .utils import safe_asanyarray, atleast2d_or_csr


class NeighborsClassifier(BaseEstimator, ClassifierMixin):
    """Classifier implementing the k-nearest neighbors (k-NN) algorithm.

    Parameters
    ----------
    n_neighbors : int, optional
        Default number of neighbors. Defaults to 5.

    leaf_size : int, optional
        Leaf size passed to BallTree.  This can affect the speed of
        the nearest neighbors query.  The optimal value depends on the
        nature of the problem (see Discussion below).  Defaults to 20.

    algorithm : {'auto', 'ball_tree', 'brute'}, optional
        Algorithm used to compute the nearest neighbors. 'ball_tree' will
        construct a BallTree while 'brute'will perform brute-force
        search. 'auto' will guess the most appropriate based on current
        dataset.  See Discussion below for notes on choosing the optimal
        method. Fitting on sparse input will override the setting of this
        parameter.

    Examples
    --------
    >>> samples = [[0, 0, 2], [1, 0, 0]]
    >>> labels = [0, 1]
    >>> from scikits.learn.neighbors import NeighborsClassifier
    >>> neigh = NeighborsClassifier(n_neighbors=1)
    >>> neigh.fit(samples, labels)
    NeighborsClassifier(algorithm='auto', leaf_size=20, n_neighbors=1)
    >>> print neigh.predict([[0,0,0]])
    [1]

    See also
    --------
    BallTree

    Discussion
    ----------
    The optimal algorithm for a given dataset is a complicated choice, and
    depends on a number of factors:
    * number of samples N (n_samples).
      * 'brute' query time grows as O[N], while
      * 'ball_tree' query time grows as O[log(N)]
    * dimensionality D (n_features)
      Intrinsic dimensionality refers to the dimension of a manifold which
      is linearly or nonlinearly embedded within the parameter space
      * 'brute' query time grows as O[D], and is unaffected by the value of d.
      * 'ball_tree' query time may grow faster or slower than this, depending
        on the structure of the data.
    * data structure: intrinsic dimensionality of the data and/or sparsity
      of the data. Intrinsic dimensionality refers to the dimension d<=D
      of a manifold on which the data lies, which can be linearly or
      nonlinearly embedded in the parameter space. Sparsity refers to the
      degree to which the data fills the parameter space (this is to be
      distinguished from the concept as used in "sparse" matrices.  The data
      matrix may have no zero entries, but the structure can still be
      "sparse" in this sense).
      * 'brute' query time is unchanged by data structure.
      * 'ball_tree' query time is greatly influenced by data structure.
        In general, the more sparse the data is, and the smaller the intrinsic
        dimension d, the faster the ball_tree algorithm will be compared to
        a brute-force search.  For data which densely fills the parameter
      space, brute-force is generally a better choice.
    * number of neighbors k requested for a query point.
      * 'brute' query time is unaffected by the value of k
      * 'ball_tree' query time is slower for k>1, mainly due to the internal
        queueing and sorting that takes place during the query.
    * leaf_size of the ball_tree
      The leaf_size parameter controls the point at which the ball_tree
      algorithm switches from a tree-based query to a brute-force search.
      As the number of points 'n' in a node becomes smaller, the difference
      between a brute-force query time O[n] and a ball-tree query time
      O[log(n)] becomes less significant.  At some point, depending on all
      the above-mentioned factors, brute-force becomes more efficient.
      The parameter leaf_size sets this threshold.
    * number of query points.
      The ball_tree algorithm requires a one-time building phase.  For very
      few queries, the time to build the tree may overwhelm any gain during
      the query time.

    Currently, the algorithm='auto' option chooses between 'brute' and
    'ball_tree' using an unsophisticated rubric.  In practice, it is
    suggested that the user experiment with different choices of
    `algorithm` and `leaf_size` to determine the optimal configuration.

    References
    ----------
    http://en.wikipedia.org/wiki/K-nearest_neighbor_algorithm
    """

    def __init__(self, n_neighbors=5, algorithm='auto', leaf_size=20):
        self.n_neighbors = n_neighbors
        self.leaf_size = leaf_size
        self.algorithm = algorithm

    def fit(self, X, y):
        """Fit the model using X, y as training data

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Training data.

        y : {array-like, sparse matrix}, shape = [n_samples]
            Target values, array of integer values.

        """
        X = safe_asanyarray(X)
        if y is None:
            raise ValueError("y must not be None")
        self._y = np.asanyarray(y)

        if issparse(X):
            self.ball_tree = None
            self._fit_X = X.tocsr()
        elif self.algorithm == 'ball_tree' or \
           (self.algorithm == 'auto' and X.shape[1] < 20):
            self.ball_tree = BallTree(X, self.leaf_size)
        else:
            self.ball_tree = None
            self._fit_X = X
        return self

    def kneighbors(self, X, n_neighbors=None, return_distance=True):
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
        NeighborsClassifier(algorithm='auto', leaf_size=20, n_neighbors=1)
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
        X = atleast2d_or_csr(X)
        if n_neighbors is None:
            n_neighbors = self.n_neighbors
        if self.ball_tree is None:
            dist = euclidean_distances(X, self._fit_X, squared=True)
            # XXX: should be implemented with a partial sort
            neigh_ind = dist.argsort(axis=1)[:, :n_neighbors]
            if not return_distance:
                return neigh_ind
            else:
                return dist.T[neigh_ind], neigh_ind
        else:
            return self.ball_tree.query(X, n_neighbors,
                                        return_distance=return_distance)

    def predict(self, X):
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
        X = atleast2d_or_csr(X)

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

    leaf_size : int, optional
        Leaf size passed to BallTree.  This can affect the speed of
        the nearest neighbors query.  The optimal value depends on the
        nature of the problem.  Defaults to 20.

    mode : {'mean', 'barycenter'}, optional
        Weights to apply to labels.

    algorithm : {'auto', 'ball_tree', 'brute'}, optional
        Algorithm used to compute the nearest neighbors. 'ball_tree' will
        construct a BallTree, while 'brute' will perform brute-force
        search. 'auto' will guess the most appropriate based on current
        dataset.

    Examples
    --------
    >>> X = [[0], [1], [2], [3]]
    >>> y = [0, 0, 1, 1]
    >>> from scikits.learn.neighbors import NeighborsRegressor
    >>> neigh = NeighborsRegressor(n_neighbors=2)
    >>> neigh.fit(X, y)
    NeighborsRegressor(algorithm='auto', leaf_size=20, mode='mean', n_neighbors=2)
    >>> print neigh.predict([[1.5]])
    [ 0.5]

    Notes
    -----
    See Discussion in NeighborsClassifier docstring regarding the optimal
    choice of algorithm and leaf_size.
    http://en.wikipedia.org/wiki/K-nearest_neighbor_algorithm
    """

    def __init__(self, n_neighbors=5, mode='mean', algorithm='auto',
                 leaf_size=20):
        self.n_neighbors = n_neighbors
        self.leaf_size = leaf_size
        self.mode = mode
        self.algorithm = algorithm

    def predict(self, X):
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
        X = atleast2d_or_csr(X)

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

def barycenter_weights(X, Z, reg=1e-3):
    """Compute barycenter weights of X from Y along the first axis

    We estimate the weights to assign to each point in Y[i] to recover
    the point X[i]. The barycenter weights sum to 1.

    Parameters
    ----------
    X : array-like, shape (n_samples, n_dim)

    Z : array-like, shape (n_samples, n_neighbors, n_dim)

    reg: float, optional
        amount of regularization to add for the problem to be
        well-posed in the case of n_neighbors > n_dim

    Returns
    -------
    B : array-like, shape (n_samples, n_neighbors)

    Notes
    -----
    See developers note for more information.
    """
    X, Z = map(np.asanyarray, (X, Z))
    n_samples, n_neighbors = X.shape[0], Z.shape[1]
    if X.dtype.kind == 'i':
        X = X.astype(np.float)
    if Z.dtype.kind == 'i':
        Z = Z.astype(np.float)
    B = np.empty((n_samples, n_neighbors), dtype=X.dtype)
    v = np.ones(n_neighbors, dtype=X.dtype)

    # this might raise a LinalgError if G is singular and has trace
    # zero
    for i, A in enumerate(Z.transpose(0, 2, 1)):
        C = A.T - X[i]  # broadcasting
        G = np.dot(C, C.T)
        trace = np.trace(G)
        if trace > 0:
            R = reg * trace
        else:
            R = reg
        G.flat[::Z.shape[1] + 1] += R
        w = linalg.solve(G, v, sym_pos=True)
        B[i, :] = w / np.sum(w)
    return B


def kneighbors_graph(X, n_neighbors, mode='connectivity', reg=1e-3):
    """Computes the (weighted) graph of k-Neighbors for points in X

    Parameters
    ----------
    X : array-like or BallTree, shape = [n_samples, n_features]
        Sample data, in the form of a numpy array or a precomputed
        :class:`BallTree`.

    n_neighbors : int
        Number of neighbors for each sample.

    mode : {'connectivity', 'distance', 'barycenter'}, optional
        Type of returned matrix: 'connectivity' will return the
        connectivity matrix with ones and zeros, in 'distance' the
        edges are Euclidean distance between points. In 'barycenter'
        they are the weights that best reconstruncts the point from
        its nearest neighbors.

    reg : float, optional
        Amount of regularization when solving the least-squares
        problem. Only relevant if mode='barycenter'. If None, use the
        default.

    Returns
    -------
    A : sparse matrix in CSR format, shape = [n_samples, n_samples]
        A[i, j] is assigned the weight of edge that connects i to j.

    Examples
    --------
    >>> X = [[0], [3], [1]]
    >>> from scikits.learn.neighbors import kneighbors_graph
    >>> A = kneighbors_graph(X, 2)
    >>> A.todense()
    matrix([[ 1.,  0.,  1.],
            [ 0.,  1.,  1.],
            [ 1.,  0.,  1.]])

    See also
    --------
    radius_neighbors_graph
    """
    if isinstance(X, BallTree):
        ball_tree = X
        X = ball_tree.data
    else:
        X = np.asanyarray(X)
        ball_tree = BallTree(X)

    n_samples = X.shape[0]
    n_nonzero = n_neighbors * n_samples
    A_indptr = np.arange(0, n_nonzero + 1, n_neighbors)

    # construct CSR matrix representation of the k-NN graph
    if mode == 'connectivity':
        A_data = np.ones((n_samples, n_neighbors))
        A_ind = ball_tree.query(
            X, k=n_neighbors, return_distance=False)

    elif mode == 'distance':
        data, ind = ball_tree.query(X, k=n_neighbors + 1)
        A_data, A_ind = data[:, 1:], ind[:, 1:]

    elif mode == 'barycenter':
        ind = ball_tree.query(
            X, k=n_neighbors + 1, return_distance=False)
        A_ind = ind[:, 1:]
        A_data = barycenter_weights(X, X[A_ind], reg=reg)

    else:
        raise ValueError(
            'Unsupported mode, must be one of "connectivity", '
            '"distance" or "barycenter" but got %s instead' % mode)

    return csr_matrix((A_data.ravel(), A_ind.ravel(), A_indptr),
                      shape=(n_samples, n_samples))


def radius_neighbors_graph(X, radius, mode='connectivity'):
    """Computes the (weighted) graph of Neighbors for points in X

    Neighborhoods are restricted the points at a distance lower than
    radius.

    Parameters
    ----------
    X : array-like or BallTree, shape = [n_samples, n_features]
        Sample data, in the form of a numpy array or a precomputed
        :class:`BallTree`.

    radius : float
        Radius of neighborhoods.

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
    >>> from scikits.learn.neighbors import radius_neighbors_graph
    >>> A = radius_neighbors_graph(X, 1.5)
    >>> A.todense()
    matrix([[ 1.,  0.,  1.],
            [ 0.,  1.,  0.],
            [ 1.,  0.,  1.]])

    See also
    --------
    kneighbors_graph
    """
    if isinstance(X, BallTree):
        ball_tree = X
        X = ball_tree.data
    else:
        X = np.asanyarray(X)
        ball_tree = BallTree(X)

    n_samples = X.shape[0]

    # construct CSR matrix representation of the NN graph
    if mode == 'connectivity':
        A_ind = ball_tree.query_radius(X, radius, return_distance=False)
        A_data = None
    elif mode == 'distance':
        A_ind, dist = ball_tree.query_radius(X, radius, return_distance=True)
        A_data = np.concatenate(list(dist))
    else:
        raise ValueError(
            'Unsupported mode, must be one of "connectivity", '
            'or "distance" but got %s instead' % mode)

    n_neighbors = np.array([len(a) for a in A_ind])
    n_nonzero = np.sum(n_neighbors)
    A_ind = np.concatenate(list(A_ind))
    A_indptr = np.concatenate((np.zeros(1, dtype=int), np.cumsum(n_neighbors)))

    if A_data is None:
        A_data = np.ones(n_nonzero)

    return csr_matrix((A_data, A_ind, A_indptr), shape=(n_samples, n_samples))
