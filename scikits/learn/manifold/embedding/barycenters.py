
"""
Computes barycenters weights from a graph and saves it in a sparse matrix
"""

import math

import scipy.sparse
from numpy import asarray, dot, eye, ones, sum, trace, zeros, intc
from numpy.linalg import solve

from .tools import create_graph, create_neighborer

__all__ = ['barycenters', ]

def barycenters(samples, neigh = None, n_neighbors = None,
    neigh_alternate_arguments = None):
    """
    Computes the barycenters of samples given as parameters and returns them.
    
    Parameters
    ----------
    samples : matrix
        The points to consider.

    n_neighbors : int
      The number of K-neighboors to use (optional, default 9) if neigh is not
      given.

    neigh : Neighbors
      A neighboorer (optional). By default, a K-Neighbor research is done.
      If provided, neigh must be a functor class . `neigh_alternate_arguments` 
      will be passed to this class constructor.

    neigh_alternate_arguments : dictionary
      Dictionary of arguments that will be passed to the `neigh` constructor

    Returns
    -------
    A CSR sparse matrix containing the barycenters weights
    """
    bary = zeros((len(samples), len(samples)))

    graph = create_graph(samples, neigh, n_neighbors,
        neigh_alternate_arguments)

    W = []
    indices=[]
    indptr=[0]
    for i in range(len(samples)):
        neighs, ind = graph[i]
        wi = barycenter(samples[i], samples[ind])
        W.extend(wi)
        indices.extend(neighs)
        indptr.append(indptr[-1] + len(neighs))

    W = asarray(W)
    indices = asarray(indices, dtype=intc)
    indptr = asarray(indptr, dtype=intc)
    return scipy.sparse.csr_matrix((W, indices, indptr), shape=(len(samples),
        len(samples)))

def barycenter(point, point_neighbors, tol = 1e-3, **kwargs):
    """
    Computes barycenter weights so that point may be reconstructed from its
    neighbors
    
    Parameters
    ----------
    point : array
        a 1D array
    
    point_neighbors : array
        a 2D array containing samples
    
    tol : float
        tolerance

    Returns
    -------
    Barycenter weights
    """
    z = point - point_neighbors
    Gram = dot(z, z.T)
    Gram += eye(len(point_neighbors), len(point_neighbors)) * tol * trace(Gram)
    wi = solve(Gram, ones(len(point_neighbors)))
    wi /= sum(wi)
    return wi
