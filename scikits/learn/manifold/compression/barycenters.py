
"""
Computes barycenters weights from a graph and saves it in a sparse matrix
"""

import math

import scipy.sparse
from numpy import asarray, dot, eye, ones, sum, trace, zeros, intc
from numpy.linalg import solve

from tools import create_graph

__all__ = ['barycenters', ]

def barycenters(samples, **kwargs):
  """
  Computes the barycenters of samples given as parameters and returns them.
  
  Parameters
  ----------
  samples : matrix
    The points to consider.

  neigh : Neighbors
    A neighboorer (optional). By default, a K-Neighbor research is done. If provided, neigh must be a functor

  neighbors : int
    The number of K-neighboors to use (optional, default 9) if neigh is not given.

  Returns
  -------
  A CSR sparse matrix containing the barycenters weights
  """
  bary = zeros((len(samples), len(samples)))

  graph = create_graph(samples, **kwargs)

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
  return scipy.sparse.csr_matrix((W, indices, indptr), shape=(len(samples), len(samples)))

def barycenter(point, neighbors, tol = 1e-3, **kwargs):
  """
  Computes barycenter weights so that point may be reconstructed from its neighbors
  
  Parameters
  ----------
  point : array
    a 1D array
  
  neighbors : array
    a 2D array containing samples
  
  Returns
  -------
  Barycenter weights
  """
  z = point - neighbors
  Gram = dot(z, z.T)
  Gram += eye(len(neighbors), len(neighbors)) * tol * trace(Gram)
  wi = solve(Gram, ones(len(neighbors)))
  wi /= sum(wi)
  return wi
