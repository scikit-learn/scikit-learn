"""Isomap for manifold learning"""

# Author: Jake Vanderplas  -- <vanderplas@astro.washington.edu>
# License: BSD, (C) JV 2011

import numpy as np
from scipy.linalg import eigh
from scipy.sparse import linalg
from ..base import BaseEstimator
from ..utils.arpack import eigsh
from ..neighbors import kneighbors_graph


def isomap(X, n_neighbors, out_dim, eigen_solver='dense',
           tol=0, max_iter=None):
    """Perform an Isomap analysis of the data

    Parameters
    ----------
    X : array-like or BallTree, shape [n_samples, n_features]
        Input data, in the form of a numpy array or a precomputed
        :class:`BallTree`.

    n_neighbors : integer
        number of neighbors to consider for each point.

    out_dim : integer
        number of coordinates for the manifold.

    eigen_solver : {'arpack', 'dense'}
        arpack can handle both dense and sparse data efficiently

    tol : float
        convergence tolerance passed to arpack or lobpcg.
        not used if eigen_solver == 'dense'

    max_iter : integer
        maximum number of iterations for the arpack solver.
        only used if eigen_solver == 'arpack'

    Returns
    -------
    Y : array-like, shape [n_samples, out_dim]
        Embedding vectors.

    References
    ----------
    [1] Tenenbaum, J.B.; De Silva, V.; & Langford, J.C. A global geometric
        framework for nonlinear dimensionality reduction.
        Science 290 (5500)
    """
    # Construct the distance graph of the data
    G = kneighbors_graph(X, n_neighbors, mode='distance')
    G = G.toarray()
    G[np.where(G == 0)] = np.inf
    # note: this could be made faster by converting to coo and
    #  calling the utility routine coo_todense with B a matrix
    #  of infinities rather than a matrix of zeros.

    n_samples = G.shape[0]

    # make sure G is symmetric.
    G = np.minimum(G, G.T)
    # note: there's probably a more efficient way to do this

    # find the matrix of minimum paths between points
    for k in xrange(n_samples):
        for i in xrange(n_samples):
            G[i] = np.minimum(G[i], G[i, k] + G[k])
    # note: as written this is an N^3 process with slow python for-loops.
    #  We should use Dijkstra's algorithm or some variant to speed this up.

    # now compute tau = -0.5 * H.(G^2).H where H = (I - 1/N)
    G **= 2
    HG = G - G.mean(0)
    HGH = HG - HG.mean(1)[:,None]
    tau = -0.5 * HGH

    # compute the out_dim largest eigenvalues and vectors of tau
    if eigen_solver == 'arpack':
        eigen_values, eigen_vectors = eigsh(tau, out_dim, which='LM',
                                            tol=tol, maxiter=max_iter)
    elif eigen_solver == 'dense':
        eigen_values, eigen_vectors = eigh(
            tau, eigvals=(n_samples - out_dim, n_samples - 1),
            overwrite_a=True)
        index = np.argsort(eigen_values)
        eigen_values = eigen_values[index]
        eigen_vectors = eigen_vectors[:, index]
    else:
        raise ValueError("Unrecognized eigen_solver '%s'" % eigen_solver)

    print eigen_values
    return eigen_vectors * np.sqrt(eigen_values)


class Isomap(BaseEstimator):
    """Isomap

    Parameters
    ----------
    n_neighbors : integer
        number of neighbors to consider for each point.

    out_dim : integer
        number of coordinates for the manifold

    Attributes
    ----------
    `embedding_vectors_` : array-like, shape [out_dim, n_samples]
        Stores the embedding vectors
    """

    def __init__(self, n_neighbors=5, out_dim=2, random_state=None):
        self.n_neighbors = n_neighbors
        self.out_dim = out_dim
        self.random_state = random_state

    def fit(self, X, Y=None, eigen_solver='dense', tol=0,
            max_iter=None, **params):
        """Compute the embedding vectors for data X

        Parameters
        ----------
        X : array-like of shape [n_samples, n_features]
            training set.

        out_dim : integer
            number of coordinates for the manifold

        eigen_solver : {'arpack', 'dense'}
            use the arpack eigensolver or a dense eigensolver based on LAPACK
            routines.

        max_iter : integer
            maximum number of iterations for the lobpcg solver.

        Returns
        -------
        self : returns an instance of self.
        """
        self.random_state = check_random_state(self.random_state)
        self._set_params(**params)
        self.embedding_, = \
            isomap(X, self.n_neighbors, self.out_dim,
                   eigen_solver=eigen_solver, tol=tol,
                   max_iter=max_iter)
        return self
