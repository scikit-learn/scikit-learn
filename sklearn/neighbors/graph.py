"""Nearest Neighbors graph functions"""

# Author: Jake Vanderplas <vanderplas@astro.washington.edu>
#
# License: BSD, (C) INRIA, University of Amsterdam

from .base import KNeighborsMixin, RadiusNeighborsMixin
from .unsupervised import NearestNeighbors


def kneighbors_graph(X, n_neighbors, mode='connectivity'):
    """Computes the (weighted) graph of k-Neighbors for points in X

    Parameters
    ----------
    X : array-like or BallTree, shape = [n_samples, n_features]
        Sample data, in the form of a numpy array or a precomputed
        :class:`BallTree`.

    n_neighbors : int
        Number of neighbors for each sample.

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
    >>> from sklearn.neighbors import kneighbors_graph
    >>> A = kneighbors_graph(X, 2)
    >>> A.todense()
    matrix([[ 1.,  0.,  1.],
            [ 0.,  1.,  1.],
            [ 1.,  0.,  1.]])

    See also
    --------
    radius_neighbors_graph
    """
    if not isinstance(X, KNeighborsMixin):
        X = NearestNeighbors(n_neighbors).fit(X)
    return X.kneighbors_graph(X._fit_X, n_neighbors, mode=mode)


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
    >>> from sklearn.neighbors import radius_neighbors_graph
    >>> A = radius_neighbors_graph(X, 1.5)
    >>> A.todense()
    matrix([[ 1.,  0.,  1.],
            [ 0.,  1.,  0.],
            [ 1.,  0.,  1.]])

    See also
    --------
    kneighbors_graph
    """
    if not isinstance(X, RadiusNeighborsMixin):
        X = NearestNeighbors(radius=radius).fit(X)
    return X.radius_neighbors_graph(X._fit_X, radius, mode)
