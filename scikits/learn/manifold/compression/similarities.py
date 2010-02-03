
"""
Computes coordinates based on the similarities given as parameters
"""

# Matthieu Brucher
# Last Change : 2008-04-11 14:42

__all__ = ['LLE', 'laplacian_maps', 'hessianMap']

from barycenters import barycenters

import numpy
import numpy.linalg as linalg
import scipy.sparse
import scipy.linalg
import scipy.sparse.linalg.dsolve
import math

def LLE(samples, nb_coords, **kwargs):
  """
  Computes the LLE reduction for a manifold
  Parameters :
    - samples are the samples that will be reduced
    - nb_coords is the number of coordinates in the manifold
    - neigh is the neighborer used (optional, default Kneighbor)
    - neighbor is the number of neighbors (optional, default 9)
  """
  W = barycenters(samples, **kwargs)
  t = numpy.eye(len(samples), len(samples)) - W
  M = numpy.asarray(numpy.dot(t.T, t))

  w, vectors = numpy.linalg.eigh(M)
  print w
  index = numpy.argsort(w)[1:1+nb_coords]

  t = scipy.sparse.eye(len(samples), len(samples)) - W
  M = t.T * t

  #sigma_solve = scipy.sparse.linalg.dsolve.splu(M).solve
  #w, vectors = scipy.sparse.linalg.eigen_symmetric(M, k=nb_coords+1, which='LR')
  #w, vectors = scipy.sparse.linalg.speigs.ARPACK_gen_eigs(M.matvec, sigma_solve, n=M.shape[0], sigma = 0, nev=nb_coords+1, which='SM')
  #vectors = numpy.asarray(vectors)
  #print w
  #index = numpy.argsort(w)[1:1+nb_coords]

  return numpy.sqrt(len(samples)) * vectors[:,index]

def laplacian_maps(samples, nb_coords, method, **kwargs):
  """
  Computes a Laplacian eigenmap for a manifold
  Parameters:
    - samples are the samples that will be reduced
    - nb_coords is the number of coordinates in the manifold
    - method is the method to create the similarity matrix
    - neigh is the neighborer used (optional, default Kneighbor)
    - neighbor is the number of neighbors (optional, default 9)
  """
  W = method(samples, **kwargs)

  if scipy.sparse.issparse(W):
    D = numpy.sqrt(W.sum(axis=0))
    Di = 1./D
    dia = scipy.sparse.dia_matrix((Di, (0,)), shape=W.shape)
    L = dia * W * dia

    w, vectors = scipy.sparse.linalg.eigen_symmetric(L, k=nb_coords+1)
    vectors = numpy.asarray(vectors)
    D = numpy.asarray(D)
    Di = numpy.asarray(Di).squeeze()

  else:
    D = numpy.sqrt(numpy.sum(W, axis=0))
    Di = 1./D
    L = Di * W * Di[:,numpy.newaxis]
    w, vectors = scipy.linalg.eigh(L)

  index = numpy.argsort(w)[-2:-2-nb_coords:-1]

  return numpy.sqrt(len(samples)) * Di[:,numpy.newaxis] * vectors[:,index] * math.sqrt(numpy.sum(D))

def laplacian_maps2(samples, nb_coords, method, **kwargs):
  """
  Computes a Laplacian eigenmap for a manifold
  Parameters:
    - samples are the samples that will be reduced
    - nb_coords is the number of coordinates in the manifold
    - method is the method to create the similarity matrix
    - neigh is the neighborer used (optional, default Kneighbor)
    - neighbor is the number of neighbors (optional, default 9)
  """
  W = method(samples, **kwargs)
  D = numpy.sum(W, axis=0)
  L = 1/D[:, numpy.newaxis] * (numpy.diag(D) - W)
  w, vectors = scipy.linalg.eig(L)
  index = numpy.argsort(w)[1:1+nb_coords]
  return numpy.sqrt(len(samples)) * vectors[:,index]

def sparse_heat_kernel(samples, kernel_width = .5, **kwargs):
  """
  Uses a heat kernel for computing similarities in a neighborhood
  """
  from tools import create_sym_graph

  graph = create_sym_graph(samples, **kwargs)

  W = []
  indices=[]
  indptr=[0]
  for i in range(len(samples)):
    neighs = graph[i]
    z = samples[i] - samples[neighs]
    wi = numpy.sum(z ** 2, axis = 1) / parameter
    W.extend(numpy.exp(-wi))
    indices.extend(neighs)
    indptr.append(indptr[-1] + len(neighs))

  W = numpy.asarray(W)
  indices = numpy.asarray(indices, dtype=numpy.intc)
  indptr = numpy.asarray(indptr, dtype=numpy.intc)
  return scipy.sparse.csr_matrix((W, indices, indptr), shape=(len(samples), len(samples)))

def heat_kernel(samples, kernel_width = .5, **kwargs):
  """
  Uses a heat kernel for computing similarities in the whole array
  """
  from tools import dist2hd
  distances = dist2hd(samples, samples)**2

  return numpy.exp(-distances/parameter)

def normalized_heat_kernel(samples, **kwargs):
  """
  Uses a heat kernel for computing similarities in the whole array
  """
  similarities = heat_kernel(samples, **kwargs)
  p1 = 1./numpy.sqrt(numpy.sum(similarities, axis=0))

  return p1[:, numpy.newaxis] * similarities * p1

def hessianMap(samples, nb_coords, **kwargs):
  """
  Computes a Hessian eigenmap for a manifold
  Parameters:
  - samples are the samples that will be reduced
  - nb_coords is the number of coordinates in the manifold
  - neigh is the neighborer used (optional, default Kneighbor)
  - neighbor is the number of neighbors (optional, default 9)
  """
  from tools import create_graph
  from numpy import linalg

  graph = create_graph(samples, **kwargs)
  dp = nb_coords * (nb_coords + 1) / 2
  W = numpy.zeros((len(samples) * dp, len(samples)))

  for i in range(len(samples)):
    neighs = graph[i]
    neighborhood = samples[neighs] - numpy.mean(samples[neighs], axis=0)
    u, s, vh = linalg.svd(neighborhood.T, full_matrices=False)
    tangent = vh.T[:,:nb_coords]

    Yi = numpy.zeros((len(tangent), dp))
    ct = 0
    for j in range(nb_coords):
      startp = tangent[:,j]
      for k in range(j, nb_coords):
        Yi[:, ct + k - j] = startp * tangent[:,k]
      ct = ct + nb_coords - j

    Yi = numpy.hstack((numpy.ones((len(neighs), 1)), tangent, Yi))

    Yt = mgs(Yi)
    Pii = Yt[:, nb_coords + 1:]
    means = numpy.mean(Pii, axis=0)[:,None]
    means[numpy.where(means < 0.0001)[0]] = 1
    W[i * dp:(i+1) * dp, neighs] = Pii.T / means

  G = numpy.dot(W.T, W)
  w, v = linalg.eigh(G)

  index = numpy.argsort(w)
  ws = w[index]
  too_small = numpy.sum(ws < 10 * numpy.finfo(numpy.float).eps)

  index = index[too_small:too_small+nb_coords]

  return numpy.sqrt(len(samples)) * v[:,index]

def mgs(A):
  """
  Computes a Gram-Schmidt orthogonalization
  """
  V = numpy.array(A)
  m, n = V.shape
  R = numpy.zeros((n, n))

  for i in range(0, n):
    R[i, i] = linalg.norm(V[:, i])
    V[:, i] /= R[i, i]
    for j in range(i+1, n):
       R[i, j] = numpy.dot(V[:, i].T, V[:, j])
       V[:, j] -= R[i, j] * V[:, i]
  return V
