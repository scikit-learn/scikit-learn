# -*- coding: utf-8 -*-

"""
Computes coordinates based on the similarities given as parameters
"""

__all__ = ['LLE', 'HessianMap']

import numpy as np
from scipy import linalg
import scipy.sparse
import math

from .base_embedding import BaseEmbedding
from ..mapping import builder as mapping_builder

from .tools import create_graph, create_sym_graph
from ...utils.euclidian_distances import euclidian_distances
from ...preprocessing import Scaler

try:
    from pyamg import smoothed_aggregation_solver
    from scipy.sparse.linalg import lobpcg
    raise RuntimeError("pyamg is currently not supported")
    pyamg_loaded = True
except:
    from scipy.sparse.linalg.eigen.arpack import eigen_symmetric
    pyamg_loaded = False

def laplacian_maps(samples, n_coords, method, **kwargs):
    """
    Computes a Laplacian eigenmap for a manifold

    Parameters
    ----------
      samples : array_like
           the samples that will be embedded
      n_coords : int
           the number of coordinates in the manifold
      method : function
           the method to create the similarity matrix
      neigh : function
           the neighborer used (optional, default Kneighbor)
      n_neighbors : int
           the number of neighbors (optional, default 9)

    Returns
    -------
    The embedding
    """
    W = method(samples, **kwargs)

    if scipy.sparse.issparse(W):
        D = np.sqrt(W.sum(axis=0))
        Di = 1./D
        dia = scipy.sparse.dia_matrix((Di, (0,)), shape=W.shape)
        L = dia * W * dia
        if pyamg_loaded:
            L = L.tocoo()
            diag_idx = (L.row == L.col)
            L.data[diag_idx] = -1
            L = L.tocsr()
            ml = smoothed_aggregation_solver(L)

            X = scipy.rand(L.shape[0], n_coords+1)

            # preconditioner based on ml
            M = ml.aspreconditioner()

            # compute eigenvalues and eigenvectors with LOBPCG
            w, vectors = lobpcg(L, X, M=M, tol=1e-3, largest=False)
        else:
            w, vectors = eigen_symmetric(L, k=n_coords+1)
        vectors = np.asarray(vectors)

    else:
        D = np.sqrt(np.sum(W, axis=0))
        Di = 1./D
        L = Di * W * Di[:,np.newaxis]
        w, vectors = scipy.linalg.eigh(L)

    Di = np.asarray(Di).squeeze()
    index = np.argsort(w)[-2:-2-n_coords:-1]
    #print index, w[index]

    return np.sqrt(len(samples)) * Di[:,np.newaxis] * vectors[:,index] * \
        math.sqrt(np.sum(D))


def sparse_heat_kernel(samples, kernel_width=.5, **kwargs):
    """
    Uses a heat kernel for computing similarities in a neighborhood

    Parameters
    ----------
      samples : array_like
           the samples that will be embedded
      method : function
           the method to create the similarity matrix
      neigh : function
           the neighborer used (optional, default Kneighbor)
      n_neighbors : int
           the number of neighbors (optional, default 9)

    Returns
    -------
    The sparse distance matrix
    """
    graph = create_sym_graph(samples, **kwargs)

    W = []
    indices=[]
    indptr=[0]
    for i in range(len(samples)):
        neighs = graph[i]
        z = samples[i] - samples[neighs]
        wi = np.sum(z ** 2, axis = 1) / kernel_width
        W.extend(np.exp(-wi))
        indices.extend(neighs)
        indptr.append(indptr[-1] + len(neighs))

    W = np.asarray(W)
    indices = np.asarray(indices, dtype=np.intc)
    indptr = np.asarray(indptr, dtype=np.intc)
    return scipy.sparse.csr_matrix((W, indices, indptr), shape=(len(samples), len(samples)))


def heat_kernel(samples, kernel_width=.5, **kwargs):
    """
    Uses a heat kernel for computing similarities in the whole array

    Parameters
    ----------
      samples : array_like
           the samples that will be embedded
      kernel_width : array_like
           the size of the heat kernel

    Returns
    -------
    The new distance
    """
    distances = euclidian_distances(samples, samples)**2

    return np.exp(-distances/kernel_width)


def normalized_heat_kernel(samples, **kwargs):
    """
    Uses a heat kernel for computing similarities in the whole array

    Parameters
    ----------
      samples : array_like
           the samples that will be embedded

    Returns
    -------
    The embedding
    """
    similarities = heat_kernel(samples, **kwargs)
    p1 = 1./np.sqrt(np.sum(similarities, axis=0))
    return p1[:, np.newaxis] * similarities * p1


def centered_normalized(samples):
    """
    Returns a set of samples that are centered and of variance 1

    >>> import numpy as np
    >>> from  scikits.learn.manifold.embedding.laplacian_map import centered_normalized
    >>> samples = np.array((0., 0., 0., \
      1., 0., 0., \
      0., 1., 0., \
      1., 1., 0., \
      0., .5, 0., \
      .5, 0., 0., \
      1., 1., 0.5, \
      )).reshape((-1,3))
    >>> centered_normalized(samples)
    array([[-1.08012345, -1.08012345, -0.40824829],
           [ 1.08012345, -1.08012345, -0.40824829],
           [-1.08012345,  1.08012345, -0.40824829],
           [ 1.08012345,  1.08012345, -0.40824829],
           [-1.08012345,  0.        , -0.40824829],
           [ 0.        , -1.08012345, -0.40824829],
           [ 1.08012345,  1.08012345,  2.44948974]])
    """
    scaler = Scaler(with_std=True)
    scaler.fit(samples)
    return scaler.transform(samples)

class LaplacianEigenmap(BaseEmbedding):
    """
    Laplacian Eigenmap embedding object

    Parameters
    ----------
    n_coords : int
      The dimension of the embedding space

    n_neighbors : int
      The number of K-neighboors to use (optional, default 9) if neigh is not
      given.

    neigh : Neighbors
      A neighboorer (optional). By default, a K-Neighbor research is done.
      If provided, neigh must be a functor class . `neigh_alternate_arguments`
      will be passed to this class constructor.

    neigh_alternate_arguments : dictionary
      Dictionary of arguments that will be passed to the `neigh` constructor

    kernel_width : float
      Width of the heat kernel

    mapping_kind : object
      The type of mapper to use. Can be:
          * None : no mapping built
          * "Barycenter" (default) : Barycenter mapping
          * a class object : a class that will be instantiated with the
              arguments of this function
          * an instance : an instance that will be fit() and then
              transform()ed

    Attributes
    ----------
    embedding_ : array_like
        BaseEmbedding of the learning data

    X_ : array_like
        Original data that is embedded

    Notes
    -----
    See also examples/plot_swissroll.py

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
    >>> laplacian = LaplacianEigenmap(n_coords=2, mapping_kind=None,\
          n_neighbors=3, kernel_width=.5)
    >>> laplacian = laplacian.fit(samples)
    """
    def __init__(self, n_coords, n_neighbors=None, neigh=None,
        neigh_alternate_arguments=None, mapping_kind="Barycenter",
        kernel_width = .5):
        BaseEmbedding.__init__(self, n_coords, n_neighbors,
            neigh,neigh_alternate_arguments, mapping_kind)
        self.kernel_width = kernel_width

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
        self.embedding_ = laplacian_maps(self.X_, n_coords=self.n_coords,
            neigh=self.neigh, n_neighbors=self.n_neighbors,
            neigh_alternate_arguments=self.neigh_alternate_arguments,
            method=sparse_heat_kernel, kernel_width=self.kernel_width)
        self.mapping = mapping_builder(self, self.mapping_kind,
            neigh=self.neigh, n_neighbors=self.n_neighbors - 1,
            neigh_alternate_arguments=self.neigh_alternate_arguments)
        return self

class DiffusionMap(BaseEmbedding):
    """
    Diffusion Map embedding object

    Parameters
    ----------
    n_coords : int
      The dimension of the embedding space

    n_neighbors : int
      The number of K-neighboors to use (optional, default 9) if neigh is not
      given.

    neigh : Neighbors
      A neighboorer (optional). By default, a K-Neighbor research is done.
      If provided, neigh must be a functor class . `neigh_alternate_arguments`
      will be passed to this class constructor.

    neigh_alternate_arguments : dictionary
      Dictionary of arguments that will be passed to the `neigh` constructor

    kernel_width : float
      Width of the heat kernel

    mapping_kind : object
      The type of mapper to use. Can be:
          * None : no mapping built
          * "Barycenter" (default) : Barycenter mapping
          * a class object : a class that will be instantiated with the
              arguments of this function
          * an instance : an instance that will be fit() and then
              transform()ed

    Attributes
    ----------
    embedding_ : array_like
        BaseEmbedding of the learning data

    X_ : array_like
        Original data that is embedded

    Notes
    -----
    See also examples/plot_swissroll.py

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
    >>> diffusion = DiffusionMap(n_coords=2, mapping_kind=None,\
          n_neighbors=3, kernel_width=.5)
    >>> diffusion = diffusion.fit(samples)
    """
    def __init__(self, n_coords, n_neighbors=None, neigh=None,
        neigh_alternate_arguments=None, mapping_kind="Barycenter",
        kernel_width=.5):
        BaseEmbedding.__init__(self, n_coords, n_neighbors,
            neigh,neigh_alternate_arguments, mapping_kind)
        self.kernel_width = kernel_width

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
        self.embedding_ = laplacian_maps(centered_normalized(self.X_),
            n_coords = self.n_coords,
            neigh = self.neigh, n_neighbors = self.n_neighbors,
            neigh_alternate_arguments = self.neigh_alternate_arguments,
            method=normalized_heat_kernel, kernel_width=self.kernel_width)
        self.mapping = mapping_builder(self, self.mapping_kind,
            neigh=self.neigh, n_neighbors=self.n_neighbors - 1,
            neigh_alternate_arguments=self.neigh_alternate_arguments)
        return self
