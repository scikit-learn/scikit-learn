"""
Routines for performing shortest-path graph searches
"""

# Author: Jake Vanderplas  -- <vanderplas@astro.washington.edu>
# License: BSD, (C) 2011

import numpy as np
cimport numpy as np

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

ITYPE = np.int32
ctypedef np.int32_t ITYPE_t


def graph_search(neighbors, distances, method='best'):
    """graph_search(N, neighbors, distances)
    
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
    neighbors = np.asarray(neighbors, dtype=ITYPE)
    distances = np.asarray(distances, dtype=DTYPE)

    assert neighbors.shape == distances.shape
    N, k = neighbors.shape

    if method=='best':
        if k < N/2:
            method = 'FW'
        else:
            method = 'D'

    graph = np.empty((N,N), dtype=DTYPE)

    if method == 'FW':
        FloydWarshall(neighbors, distances, graph)

    elif method == 'D':
        Dijkstra(neighbors, distances, graph)

    else:
        raise ValueError("unrecognized method '%s'" % method)

    return graph


cdef void FloydWarshall(np.ndarray[ITYPE_t, ndim=2] neighbors,
                        np.ndarray[DTYPE_t, ndim=2] distances,
                        np.ndarray[DTYPE_t, ndim=2] graph):
    cdef int N = graph.shape[0]
    cdef int k = neighbors.shape[1]
    cdef int i, j, m

    #initialize all distances to infinity
    graph[:,:] = np.inf

    #first populate the graph with the given distances
    for i from 0 <= i < N:
        for j from 0 <= j < k:
            m = neighbors[i,j]
            graph[i, m] = distances[i, j]
            graph[m, i] = distances[i, j]

    #now perform the Floyd-Warshall algorithm
    for k from 0 <= k < N:
        for i from 0 <= i < N:
            for j from 0 <= j < N:
                graph[i, j] = min(graph[i, j], graph[i, k] + graph[k, j])
                
cdef void Dijkstra(np.ndarray[ITYPE_t, ndim=2] neighbors,
                   np.ndarray[DTYPE_t, ndim=2] distances,
                   np.ndarray[DTYPE_t, ndim=2] graph):
    #place-holder for right now: code this up later
    FloydWarshall(neighbors, distances, graph)
