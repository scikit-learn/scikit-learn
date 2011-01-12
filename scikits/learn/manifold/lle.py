# -*- coding: utf-8 -*-

"""
Computes coordinates based on the similarities given as parameters
"""

__all__ = ['LLE']

import numpy as np
from scipy import linalg

from .base_embedding import BaseEmbedding
from ..neighbors import kneighbors_graph


class LLE(BaseEmbedding):
    """
    LLE embedding object

    Parameters
    ----------
    n_coords : int
      The dimension of the embedding space

    n_neighbors : int
      The number of K-neighboors to use (optional, default 9) if neigh is not
      given.

    ball_tree : BallTree
      A neighboorer (optional). By default, a K-Neighbor research is done.
      If provided, neigh must provide a query method.

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

    .. [1] Sam T. Roweis  and Lawrence K. Saul,
           "Nonlinear Dimensionality Reduction by Locally Linear BaseEmbedding",
           Science, Vol. 290. no. 5500, pp. 2323 -- 2326, 22 December 2000

    Examples
    --------
    >>> from scikits.learn.manifold import LLE
    >>> import numpy as np
    >>> samples = np.array((0., 0., 0., \
      1., 0., 0., \
      0., 1., 0., \
      1., 1., 0., \
      0., .5, 0., \
      .5, 0., 0., \
      1., 1., 0.5, \
      )).reshape((-1,3))
    >>> lle = LLE(n_coords=2, n_neighbors=3)
    >>> lle = lle.fit(samples)
    """
    def __init__(self, n_coords, n_neighbors=None, ball_tree=None):
        BaseEmbedding.__init__(self, n_coords, n_neighbors, ball_tree)

    def fit(self, X):
        """
        Parameters
        ----------
        X : array_like
        The learning dataset

        Returns
        -------
        Self
        """
        self.X_ = np.asanyarray(X)
        W = kneighbors_graph(self.X_, ball_tree=self.ball_tree, drop_first=True,
            n_neighbors=self.n_neighbors, weight="barycenter")

        t = np.eye(len(self.X_), len(self.X_)) - W
        M = np.asarray(np.dot(t.T, t))

        w, vectors = linalg.eigh(M)
        index = np.argsort(w)[1:1+self.n_coords]

        self.embedding_ = np.sqrt(len(self.X_)) * vectors[:,index]
        return self
