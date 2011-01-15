# -*- coding: utf-8 -*-

"""
Computes coordinates based on the similarities given as parameters
"""

__all__ = ['LaplacianEigenmap']

from math import sqrt
import numpy as np
from scipy import sparse

from ..metrics.pairwise import euclidian_distances
from ..neighbors import kneighbors_graph

from .base_embedding import BaseEmbedding

try:
    from pyamg import smoothed_aggregation_solver
    from scipy.sparse.linalg import lobpcg
    raise RuntimeError("pyamg is currently not supported")
    pyamg_loaded = True
except:
    from ..utils.fixes import arpack_eigsh
    pyamg_loaded = False


def find_largest_eigenvectors(L, n_coords, largest=True):
    """Finds the first eigenvalues/eigenvectors of a matrix

    It returns both eigenvalues and eigenvectors

    Parameters
    ----------
    L : ndarray or sparse matrix

    n_coords : int
        The number of eigenvalues to compute.

    largest : bool
        False if the lowest eigenvalues must be computed

    Returns
    -------
    w : array_like
        XXX ?

    vectors : array_like
        XXX ?
    """
    if pyamg_loaded and sparse.isspmatrix(L):
        L = L.tocoo()
        diag_idx = (L.row == L.col)
        L.data[diag_idx] = -1
        L = L.tocsr()
        ml = smoothed_aggregation_solver(L)

        X = np.random.rand(L.shape[0], n_coords + 1)

        # preconditioner based on ml
        M = ml.aspreconditioner()

        # compute eigenvalues and eigenvectors with LOBPCG
        return lobpcg(L, X, M=M, tol=1e-3, largest=not largest)
    else:
        return arpack_eigsh(L, k=n_coords+1, which='LM' if largest else 'SM')


def _laplacian_maps(W, n_coords):
    """Computes a Laplacian eigenmap given a similarity matrix

    Parameters
    ----------
    W : ndarray or sparse matrix of shape [n_samples, n_samples]
        The similarity (or weight) matrix
    n_coords : int
        the number of coordinates in the manifold

    Returns
    -------
    E : array of shape [n_samples, n_coords]
        The coordinates of all samples in the lower dimensional space
    """
    D = np.sqrt(W.sum(axis=0))
    Di = 1. / D
    dia = sparse.dia_matrix((Di, (0, )), shape=W.shape)
    L = dia * W * dia

    w, vectors = find_largest_eigenvectors(L, n_coords)
    vectors = np.asarray(vectors)

    Di = np.asarray(Di).squeeze()
    index = np.argsort(w)[-2:-2 - n_coords:-1]

    return sqrt(W.shape[0]*np.sum(D)) * Di[:, np.newaxis] * vectors[:, index]


def sparse_heat_kernel(samples, kernel_width=.5, **kwargs):
    """Uses a heat kernel for computing similarities in a neighborhood

    Parameters
    ----------
    samples : array_like
        the samples that will be embedded
    method : function
        the method to create the similarity matrix
    ball_tree : function
        the neighborer used (optional, default BallTree)
    n_neighbors : int
        the number of neighbors (optional, default 9)

    Returns
    -------
    graph : sparse matrix
        The sparse distance matrix
    """
    graph = kneighbors_graph(samples, weight="distance", drop_first = True,
                             **kwargs)
    tri = graph.T.tocsr()
    for i, ind in enumerate(zip(*tri.nonzero())):
        graph[ind] = tri.data[i]

    graph.data = np.exp(-graph.data / kernel_width)
    return graph


def heat_kernel(X, kernel_width=.5, **kwargs):
    """Uses a heat kernel for computing similarities in the whole array

    Parameters
    ----------
    samples : array_like
        the samples that will be embedded
    kernel_width : array_like
        the size of the heat kernel

    Returns
    -------
    D : ndarray
        The similarity matrix.
    """
    distances = euclidian_distances(X, X) ** 2
    return np.exp(-distances / kernel_width)


def normalized_heat_kernel(samples, **kwargs):
    """Uses a heat kernel for computing similarities in the whole array

    Parameters
    ----------
    samples : array_like
           the samples that will be embedded

    Returns
    -------
    D : ndarray
        The similarity matrix.
    """
    similarities = heat_kernel(samples, **kwargs)
    p1 = 1. / np.sqrt(np.sum(similarities, axis=0))
    return p1[:, np.newaxis] * similarities * p1


class LaplacianEigenmap(BaseEmbedding):
    """Laplacian Eigenmap embedding object

    Parameters
    ----------
    n_coords : int
        The dimension of the embedding space
    n_neighbors : int
        The number of K-neighboors to use (optional, default 9)
        if neigh is not given.
    ball_tree : BallTree
        A neighboorer (optional). By default, a K-Neighbor research is done.
        If provided, neigh must provide a query method.
    kernel_width : float
        Width of the heat kernel
    method : function
        the method to create the similarity matrix
    ball_tree : function
        the neighborer used (optional, default BallTree)
    n_neighbors : int
        the number of neighbors (optional, default 9)

    Attributes
    ----------
    embedding_ : array_like
        BaseEmbedding of the learning data

    X_ : array_like
        Original data that is embedded

    Notes
    -----
    See examples/plot_swissroll.py and examples/plot_manifold_embeddings.py
    for an example.

    .. [1] Partha Niyogi, andMikhail Belkin,
           "Laplacian Eigenmaps for Dimensionality Reduction and Data
           Representation",
           Neural Computation, Vol. 15, No. 6, Pages 1373-1396
           doi: 10.1162/089976603321780317

    Examples
    --------
    >>> from scikits.learn.manifold import LaplacianEigenmap
    >>> import numpy as np
    >>> samples = np.array((0., 0., 0., \
      1., 0., 0., \
      0., 1., 0., \
      1., 1., 0., \
      0., .5, 0., \
      .5, 0., 0., \
      1., 1., 0.5, \
      )).reshape((-1,3))
    >>> laplacian = LaplacianEigenmap(n_coords=2, n_neighbors=3, \
                                      kernel_width=.5)
    >>> laplacian = laplacian.fit(samples)
    """

    def __init__(self, n_coords, n_neighbors=None, ball_tree=None,
        kernel_width=.5):
        BaseEmbedding.__init__(self, n_coords, n_neighbors, ball_tree)
        self.kernel_width = kernel_width

    def fit(self, X):
        """Compute the embedding

        Parameters
        ----------
        X : array_like
        The learning dataset

        Returns
        -------
        Self
        """
        X_ = np.asanyarray(X)
        W = sparse_heat_kernel(X_, ball_tree=self.ball_tree,
                               n_neighbors=self.n_neighbors,
                               kernel_width=self.kernel_width)
        self.embedding_ = _laplacian_maps(W, self.n_coords)
        return self


class DiffusionMap(BaseEmbedding):
    """Diffusion Map embedding object

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

    kernel_width : float
      Width of the heat kernel

    Attributes
    ----------
    embedding_ : array_like
        BaseEmbedding of the learning data

    X_ : array_like
        Original data that is embedded

    Notes
    -----
    See examples/plot_swissroll.py and examples/plot_manifold_embeddings.py
    for an example.

    .. [1] R.R. Coifman, and S. Lafon, "Diffusion maps",
           Applied and Computational Harmonic Analysis,
           Vol 21, July 2006, pp 5-30

    Examples
    --------
    >>> from scikits.learn.manifold import DiffusionMap
    >>> import numpy as np
    >>> samples = np.array((0., 0., 0., \
      1., 0., 0., \
      0., 1., 0., \
      1., 1., 0., \
      0., .5, 0., \
      .5, 0., 0., \
      1., 1., 0.5, \
      )).reshape((-1,3))
    >>> diffusion = DiffusionMap(n_coords=2, n_neighbors=3, kernel_width=.5)
    >>> diffusion = diffusion.fit(samples)
    """

    def __init__(self, n_coords, n_neighbors=None, ball_tree=None,
        kernel_width=.5):
        BaseEmbedding.__init__(self, n_coords, n_neighbors, ball_tree)
        self.kernel_width = kernel_width

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
        X_ = X_ - np.mean(X_, axis=0) # center
        X_ /= np.std(X_, axis=0) # normalize
        W = normalized_heat_kernel(X_, ball_tree=self.ball_tree,
                                   n_neighbors=self.n_neighbors,
                                   kernel_width=self.kernel_width)
        self.embedding_ = _laplacian_maps(W, self.n_coords)
        return self
