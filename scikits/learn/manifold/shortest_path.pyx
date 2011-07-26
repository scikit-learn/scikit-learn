"""
Routines for performing shortest-path graph searches
"""

# Author: Jake Vanderplas  -- <vanderplas@astro.washington.edu>
# License: BSD, (C) 2011

import numpy as np
cimport numpy as np

cimport cython

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

ITYPE = np.int32
ctypedef np.int32_t ITYPE_t

cdef inline DTYPE_t fmin(DTYPE_t a, DTYPE_t b):
   if a <= b:
       return a
   else:
       return b

def shortest_path(neighbors, distances, method='best'):
    """shortest_path(N, neighbors, distances)
    
    Perform a shortest-path graph search on data

    Parameters
    ----------
    neighbors : array-like, integer, shape = [N, k]
        indices of edges.  The vector neighbors[i,:] gives the indices of
        the k neighbors which are connected to point i
    distances : array-like, float, shape = [N, k]
        distances along edges.  distances[i,j] gives the distance from
        point j to point neighbors[i,j]
    method : string
        method to use.  Options are
	'best' : attempt to choose the best method
	'FW' : Floyd-Warshall algorithm.  O[N^3]
	'D' : Dijkstra's algorithm with Fibonacci stacks.  O[(k+log(N))N^2]

    Returns
    -------
    G : np.ndarray, float, shape = [N,N]
        G[i,j] gives the shortest distance from point i to point j
        along the graph.
    """
    neighbors = np.asarray(neighbors, dtype=ITYPE, order='C')
    distances = np.asarray(distances, dtype=DTYPE, order='C')

    assert neighbors.shape == distances.shape
    N, k = neighbors.shape

    graph = np.empty((N,N), dtype=DTYPE, order='C')

    if method=='best':
        if k < N/2:
            method = 'FW'
        else:
            method = 'D'

    if method == 'FW':
        FloydWarshall(neighbors, distances, graph)

    elif method == 'D':
        Dijkstra(neighbors, distances, graph)

    else:
        raise ValueError("unrecognized method '%s'" % method)

    return graph


@cython.boundscheck(False)
cdef void FloydWarshall(np.ndarray[ITYPE_t, ndim=2, mode='c'] neighbors,
                        np.ndarray[DTYPE_t, ndim=2, mode='c'] distances,
                        np.ndarray[DTYPE_t, ndim=2, mode='c'] graph):
    cdef int N = graph.shape[0]
    cdef int n_neighbors = neighbors.shape[1]
    cdef unsigned int i, j, k, m

    cdef DTYPE_t infinity = np.inf
    cdef DTYPE_t sum_ijk

    #initialize all distances to infinity
    graph[:,:] = np.inf

    #graph[i,i] should be zero
    graph.flat[::N+1] = 0

    #first populate the graph with the given distances
    for i from 0 <= i < N:
        for j from 0 <= j < n_neighbors:
            m = neighbors[i,j]
            graph[i, m] = distances[i, j]
            graph[m, i] = distances[i, j]

    #now perform the Floyd-Warshall algorithm
    for k from 0 <= k < N:
        for i from 0 <= i < N:
            if graph[i,k] == infinity:
                continue
            for j from 0 <= j < N:
                sum_ijk = graph[i, k] + graph[k, j]
                if sum_ijk < graph[i,j]:
                    graph[i,j] = sum_ijk


@cython.boundscheck(False)
cdef void Dijkstra(np.ndarray[ITYPE_t, ndim=2] neighbors,
                   np.ndarray[DTYPE_t, ndim=2] distances,
                   np.ndarray[DTYPE_t, ndim=2] graph):
    #place-holder for right now: code this up later
    FloydWarshall(neighbors, distances, graph)
