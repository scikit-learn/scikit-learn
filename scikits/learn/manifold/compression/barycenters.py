
"""
Computes barycenters weights from a graph and saves it in a sparse matrix
"""

# Matthieu Brucher
# Last Change : 2008-02-28 14:06

import math

import scipy.sparse
from numpy import asarray, dot, eye, ones, sum, trace, zeros, intc
from numpy.linalg import solve

from tools import create_graph

__all__ = ['barycenters', ]

def barycenters(samples, **kwargs):
  """
  Computes the barycenters of samples given as parameters and returns them.
  """
  bary = zeros((len(samples), len(samples)))

  graph = create_graph(samples, **kwargs)

  tol = 1e-3 #math.sqrt(finfo(samples.dtype).eps)

  W = []
  indices=[]
  indptr=[0]
  for i in range(len(samples)):
    neighs, ind = graph[i]
    z = samples[i] - samples[ind]
    Gram = dot(z, z.T)
    Gram += eye(len(neighs), len(neighs)) * tol * trace(Gram)
    wi = solve(Gram, ones(len(neighs)))
    wi /= sum(wi)
    W.extend(wi)
    indices.extend(neighs)
    indptr.append(indptr[-1] + len(neighs))

  W = asarray(W)
  indices = asarray(indices, dtype=intc)
  indptr = asarray(indptr, dtype=intc)
  return scipy.sparse.csr_matrix((W, indices, indptr), shape=(len(samples), len(samples)))
