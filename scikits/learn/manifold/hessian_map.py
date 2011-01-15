# -*- coding: utf-8 -*-

"""
Computes coordinates based on the similarities given as parameters
"""

__all__ = ['HessianMap']

from math import sqrt

import numpy as np
from scipy import linalg

from .base_embedding import BaseEmbedding
from ..neighbors import kneighbors_graph


def hessian_map(X, n_coords, ball_tree, n_neighbors=9):
    """Computes a Hessian eigenmap for a manifold

    Parameters
    ----------
    X : array_like
        the samples that will be embedded
    n_coords : int
        the number of coordinates in the manifold
    ball_tree : callable
        the neighborer used (optional, default BallTree)
    n_neighbors : int
        the number of neighbors (optional, default 9)

    Returns
    -------
    E : array of shape [n_samples, n_coords]
        The coordinates of all samples in the lower dimensional space
    """
    graph = kneighbors_graph(X, ball_tree=ball_tree,
                                      n_neighbors=n_neighbors)
    graph = graph.tolil()
    n_samples = X.shape[0]
    dp = n_coords * (n_coords + 1) / 2
    W = np.zeros((n_samples * dp, n_samples))

    for i in range(len(graph.rows)):
        neighs = graph.rows[i]
        neighborhood = X[neighs] - np.mean(X[neighs], axis=0)
        U, s, Vh = linalg.svd(neighborhood.T, full_matrices=False)
        tangent = Vh.T[:, :n_coords]

        Yi = np.zeros((len(tangent), dp))
        ct = 0
        for j in range(n_coords):
            startp = tangent[:, j]
            for k in range(j, n_coords):
                Yi[:, ct + k - j] = startp * tangent[:, k]
            ct += n_coords - j

        Yi = np.hstack((np.ones((len(neighs), 1)), tangent, Yi))

        Yt = mgs(Yi)
        Pii = Yt[:, n_coords + 1:]
        means = np.mean(Pii, axis=0)[:, None]
        means[np.where(means < 0.0001)[0]] = 1
        W[i * dp:(i + 1) * dp, neighs] = Pii.T / means

    G = np.dot(W.T, W)
    w, V = linalg.eigh(G)

    index = np.argsort(w)
    ws = w[index]
    too_small = np.sum(ws < 10 * np.finfo(np.float).eps)

    index = index[too_small:too_small + n_coords]

    return sqrt(n_samples) * V[:, index]


class HessianMap(BaseEmbedding):
    """Hessian Map embedding object

    Parameters
    ----------
    n_coords : int
      The dimension of the embedding space

    n_neighbors : int
      The number of K-neighboors to use (optional, default 9) if neigh is not
      given.

    ball_tree : BallTree
      A neighboorer (optional). By default, a K-Neighbor research is done.
      If provided, neigh must be a functor class . `neigh_alternate_arguments`
      will be passed to this class constructor.

    Attributes
    ----------
    embedding_ : array_like
        BaseEmbedding of the learning data

    X_ : array_like
        Original data that is embedded

    See also
    --------
    See examples/plot_swissroll.py and examples/plot_manifold_embeddings.py
    for an example.

    .. [1] David L. Donoho and Carrie Grimes,
           "Hessian Eigenmaps: new locally linear embedding techniques for
            high-dimensional data",

    Examples
    --------
    >>> from scikits.learn.manifold import HessianMap
    >>> import numpy as np
    >>> samples = np.array((0., 0., 0., \
      1., 0., 0., \
      0., 1., 0., \
      1., 1., 0., \
      0., .5, 0., \
      .5, 0., 0., \
      1., 1., 0.5, \
      )).reshape((-1,3))
    >>> hessian = HessianMap(n_coords=2, n_neighbors=4)
    >>> hessian = hessian.fit(samples)
    """

    def fit(self, X):
        """Compute the embedding

        Parameters
        ----------
        X : array_like
            The learning dataset

        Returns
        -------
        self
        """
        X_ = np.asanyarray(X)
        self.embedding_ = hessian_map(X_, n_coords=self.n_coords,
            ball_tree=self.ball_tree, n_neighbors=self.n_neighbors)
        return self


def mgs(A):
    """Computes a Gram-Schmidt orthogonalization

    The matrix A is modified inplace

    Parameters
    ----------
    A : array_like

    Returns
    -------
    A : array_like
        The orthogonalization of A
    """
    A = np.asanyarray(A)
    m, n = A.shape
    R = np.zeros((n, n))

    for i in range(0, n):
        R[i, i] = linalg.norm(A[:, i])
        A[:, i] /= R[i, i]
        for j in range(i + 1, n):
            R[i, j] = np.dot(A[:, i].T, A[:, j])
            A[:, j] -= R[i, j] * A[:, i]
    return A
