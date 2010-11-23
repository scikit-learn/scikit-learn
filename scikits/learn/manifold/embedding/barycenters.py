# -*- coding: utf-8 -*-

"""
Computes barycenters weights from a graph and saves it in a sparse matrix
"""

import numpy as np
from scipy import linalg
from scipy import sparse

from .tools import create_graph

__all__ = ['barycenters', ]

def barycenters(samples, neigh=None, n_neighbors=None,
    neigh_alternate_arguments=None):
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
    graph = create_graph(samples, neigh, n_neighbors,
        neigh_alternate_arguments)

    W = []
    indices=[]
    indptr=[0]
    for i in range(len(samples)):
        neighs, ind = graph[i]
        ind = ind[1:]
        wi = barycenter_weights(samples[i], samples[ind])
        W.extend(wi)
        indices.extend(ind)
        indptr.append(indptr[-1] + len(ind))

    W = np.asarray(W)
    indices = np.asarray(indices, dtype=np.intc)
    indptr = np.asarray(indptr, dtype=np.intc)
    return sparse.csr_matrix((W, indices, indptr), shape=(len(samples),
        len(samples)))


def barycenter_weights(point, point_neighbors, tol=1e-3, **kwargs):
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
    Gram = np.dot(z, z.T)
    Gram += np.eye(len(point_neighbors), len(point_neighbors)) * tol *\
       np.trace(Gram)
    wi = linalg.solve(Gram, np.ones(len(point_neighbors)))
    wi /= np.sum(wi)
    return wi
