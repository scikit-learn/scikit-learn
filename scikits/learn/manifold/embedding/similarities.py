# -*- coding: utf-8 -*-

"""
Computes coordinates based on the similarities given as parameters
"""

__all__ = ['LLE', 'HessianMap']

from .base_embedding import BaseEmbedding
from ..mapping import builder as mapping_builder

from barycenters import barycenters

import numpy as np
from scipy import linalg
import scipy.sparse
from scipy.sparse.linalg.eigen.arpack import eigen_symmetric
import math
from .tools import create_graph, create_sym_graph
from ...cluster.spectral import spectral_embedding

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

    neigh : Neighbors
      A neighboorer (optional). By default, a K-Neighbor research is done.
      If provided, neigh must be a functor class . `neigh_alternate_arguments`
      will be passed to this class constructor.

    neigh_alternate_arguments : dictionary
      Dictionary of arguments that will be passed to the `neigh` constructor

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

    See Also
    --------


    Notes
    -----

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
    >>> lle = LLE(n_coords=2, mapping_kind=None, n_neighbors=3)
    >>> lle = lle.fit(samples)
    """
    def __init__(self, n_coords, n_neighbors=None, neigh=None,
        neigh_alternate_arguments=None, mapping_kind="Barycenter"):
        BaseEmbedding.__init__(self, n_coords, n_neighbors,
            neigh,neigh_alternate_arguments, mapping_kind)

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
        W = barycenters(self.X_, neigh=self.neigh,
            n_neighbors=self.n_neighbors,
            neigh_alternate_arguments=self.neigh_alternate_arguments)
        t = np.eye(len(self.X_), len(self.X_)) - W
        M = np.asarray(np.dot(t.T, t))

        w, vectors = linalg.eigh(M)
        index = np.argsort(w)[1:1+self.n_coords]

        t = scipy.sparse.eye(len(self.X_), len(self.X_)) - W
        M = t.T * t

        self.embedding_ = np.sqrt(len(self.X_)) * vectors[:,index]
        self.mapping = mapping_builder(self, self.mapping_kind,
            neigh=self.neigh, n_neighbors=self.n_neighbors - 1,
            neigh_alternate_arguments=self.neigh_alternate_arguments)
        return self

def laplacian_maps(samples, n_coords, method, mode=None, **kwargs):
    """
    Computes a Laplacian eigenmap for a manifold
    Parameters:
      - samples are the samples that will be reduced
      - n_coords is the number of coordinates in the manifold
      - method is the method to create the similarity matrix
      - neigh is the neighborer used (optional, default Kneighbor)
      - neighbor is the number of neighbors (optional, default 9)
    """
    W = method(samples, **kwargs)

    if scipy.sparse.issparse(W):
        D = np.sqrt(W.sum(axis=0))
        Di = 1./D
        dia = scipy.sparse.dia_matrix((Di, (0,)), shape=W.shape)
        L = dia * W * dia

        w, vectors = eigen_symmetric(L, k=n_coords+1)
        vectors = np.asarray(vectors)
        D = np.asarray(D)
        Di = np.asarray(Di).squeeze()

    else:
        D = np.sqrt(np.sum(W, axis=0))
        Di = 1./D
        L = Di * W * Di[:,np.newaxis]
        w, vectors = scipy.linalg.eigh(L)

    index = np.argsort(w)[-2:-2-n_coords:-1]

    return np.sqrt(len(samples)) * Di[:,np.newaxis] * vectors[:,index] * \
        math.sqrt(np.sum(D))


def sparse_heat_kernel(samples, kernel_width=.5, **kwargs):
    """
    Uses a heat kernel for computing similarities in a neighborhood
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
    """
    from tools import dist2hd
    distances = dist2hd(samples, samples)**2

    return np.exp(-distances/kernel_width)


def normalized_heat_kernel(samples, **kwargs):
    """
    Uses a heat kernel for computing similarities in the whole array
    """
    similarities = heat_kernel(samples, **kwargs)
    p1 = 1./np.sqrt(np.sum(similarities, axis=0))
    return p1[:, np.newaxis] * similarities * p1


def hessian_map(samples, n_coords, **kwargs):
    """
    Computes a Hessian eigenmap for a manifold
    Parameters:
    - samples are the samples that will be reduced
    - n_coords is the number of coordinates in the manifold
    - neigh is the neighborer used (optional, default Kneighbor)
    - neighbor is the number of neighbors (optional, default 9)
    """
    graph = create_graph(samples, **kwargs)
    dp = n_coords * (n_coords + 1) / 2
    W = np.zeros((len(samples) * dp, len(samples)))

    for i in range(len(samples)):
        neighs = graph[i][1]
        neighborhood = samples[neighs] - np.mean(samples[neighs], axis=0)
        u, s, vh = linalg.svd(neighborhood.T, full_matrices=False)
        tangent = vh.T[:,:n_coords]

        Yi = np.zeros((len(tangent), dp))
        ct = 0
        for j in range(n_coords):
            startp = tangent[:,j]
            for k in range(j, n_coords):
                Yi[:, ct + k - j] = startp * tangent[:,k]
            ct = ct + n_coords - j

        Yi = np.hstack((np.ones((len(neighs), 1)), tangent, Yi))

        Yt = mgs(Yi)
        Pii = Yt[:, n_coords + 1:]
        means = np.mean(Pii, axis=0)[:,None]
        means[np.where(means < 0.0001)[0]] = 1
        W[i * dp:(i+1) * dp, neighs] = Pii.T / means

    G = np.dot(W.T, W)
    w, v = linalg.eigh(G)

    index = np.argsort(w)
    ws = w[index]
    too_small = np.sum(ws < 10 * np.finfo(np.float).eps)

    index = index[too_small:too_small+n_coords]

    return np.sqrt(len(samples)) * v[:,index]


class HessianMap(BaseEmbedding):
    """
    Hessian Map embedding object

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

    See Also
    --------


    Notes
    -----

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
    >>> hessian = HessianMap(n_coords = 2, mapping_kind = None, n_neighbors = 4)
    >>> hessian = hessian.fit(samples)
    """
    def __init__(self, n_coords, n_neighbors=None, neigh=None,
        neigh_alternate_arguments=None, mapping_kind="Barycenter"):
        BaseEmbedding.__init__(self, n_coords, n_neighbors,
            neigh,neigh_alternate_arguments, mapping_kind)

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
        self.embedding_ = hessian_map(self.X_, n_coords=self.n_coords,
            neigh=self.neigh, n_neighbors=self.n_neighbors,
            neigh_alternate_arguments=self.neigh_alternate_arguments)
        self.mapping = mapping_builder(self, self.mapping_kind,
            neigh=self.neigh, n_neighbors=self.n_neighbors - 1,
            neigh_alternate_arguments=self.neigh_alternate_arguments)
        return self

def mgs(A):
    """
    Computes a Gram-Schmidt orthogonalization
    """
    V = np.array(A)
    m, n = V.shape
    R = np.zeros((n, n))

    for i in range(0, n):
        R[i, i] = linalg.norm(V[:, i])
        V[:, i] /= R[i, i]
        for j in range(i+1, n):
            R[i, j] = np.dot(V[:, i].T, V[:, j])
            V[:, j] -= R[i, j] * V[:, i]
    return V
