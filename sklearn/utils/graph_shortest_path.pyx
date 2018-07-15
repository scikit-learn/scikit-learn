"""
Routines for performing shortest-path graph searches

The main interface is in the function `graph_shortest_path`.  This
calls cython routines that compute the shortest path using either
the Floyd-Warshall algorithm, or Dijkstra's algorithm with Fibonacci Heaps.
"""

# Author: Jake Vanderplas  -- <vanderplas@astro.washington.edu>
# License: BSD 3 clause, (C) 2011

import numpy as np
cimport numpy as np

from scipy.sparse import csr_matrix, isspmatrix, isspmatrix_csr

cimport cython

from libc.stdlib cimport malloc, free

np.import_array()

cimport fibonacci
from fibonacci cimport (FibonacciHeap, FibonacciNode, initialize_node,
                        rightmost_sibling, leftmost_sibling, add_child,
                        add_sibling, remove, insert_node, decrease_val,
                        link, remove_min)

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

ITYPE = np.int32
ctypedef np.int32_t ITYPE_t


def graph_shortest_path(dist_matrix, directed=True, method='auto'):
    """
    Perform a shortest-path graph search on a positive directed or
    undirected graph.

    Parameters
    ----------
    dist_matrix : arraylike or sparse matrix, shape = (N,N)
        Array of positive distances.
        If vertex i is connected to vertex j, then dist_matrix[i,j] gives
        the distance between the vertices.
        If vertex i is not connected to vertex j, then dist_matrix[i,j] = 0
    directed : boolean
        if True, then find the shortest path on a directed graph: only
        progress from a point to its neighbors, not the other way around.
        if False, then find the shortest path on an undirected graph: the
        algorithm can progress from a point to its neighbors and vice versa.
    method : string ['auto'|'FW'|'D']
        method to use.  Options are
        'auto' : attempt to choose the best method for the current problem
        'FW' : Floyd-Warshall algorithm.  O[N^3]
        'D' : Dijkstra's algorithm with Fibonacci stacks.  O[(k+log(N))N^2]

    Returns
    -------
    G : np.ndarray, float, shape = [N,N]
        G[i,j] gives the shortest distance from point i to point j
        along the graph.

    Notes
    -----
    As currently implemented, Dijkstra's algorithm does not work for
    graphs with direction-dependent distances when directed == False.
    i.e., if dist_matrix[i,j] and dist_matrix[j,i] are not equal and
    both are nonzero, method='D' will not necessarily yield the correct
    result.

    Also, these routines have not been tested for graphs with negative
    distances.  Negative distances can lead to infinite cycles that must
    be handled by specialized algorithms.
    """
    if not isspmatrix_csr(dist_matrix):
        dist_matrix = csr_matrix(dist_matrix)

    N = dist_matrix.shape[0]
    Nk = len(dist_matrix.data)

    if method == 'auto':
        if Nk < N * N / 4:
            method = 'D'
        else:
            method = 'FW'

    if method == 'FW':
        graph = np.asarray(dist_matrix.toarray(), dtype=DTYPE, order='C')
        floyd_warshall(graph, directed)
    elif method == 'D':
        graph = np.zeros((N, N), dtype=DTYPE, order='C')
        dijkstra(dist_matrix, graph, directed)
    else:
        raise ValueError("unrecognized method '%s'" % method)

    return graph


@cython.boundscheck(False)
cdef np.ndarray floyd_warshall(np.ndarray[DTYPE_t, ndim=2, mode='c'] graph,
                              int directed=0):
    """
    FloydWarshall algorithm

    Parameters
    ----------
    graph : ndarray
        on input, graph is the matrix of distances between connected points.
        unconnected points have distance=0
        on exit, graph is overwritten with the output
    directed : bool, default = False
        if True, then the algorithm will only traverse from a point to
        its neighbors when finding the shortest path.
        if False, then the algorithm will traverse all paths in both
        directions.

    Returns
    -------
    graph : ndarray
        the matrix of shortest paths between points.
        If no path exists, the path length is zero
    """
    cdef int N = graph.shape[0]
    assert graph.shape[1] == N

    cdef unsigned int i, j, k, m

    cdef DTYPE_t infinity = np.inf
    cdef DTYPE_t sum_ijk

    #initialize all distances to infinity
    graph[np.where(graph == 0)] = infinity

    #graph[i,i] should be zero
    graph.flat[::N + 1] = 0

    # for a non-directed graph, we need to symmetrize the distances
    if not directed:
        for i from 0 <= i < N:
            for j from i + 1 <= j < N:
                if graph[j, i] <= graph[i, j]:
                    graph[i, j] = graph[j, i]
                else:
                    graph[j, i] = graph[i, j]

    #now perform the Floyd-Warshall algorithm
    for k from 0 <= k < N:
        for i from 0 <= i < N:
            if graph[i, k] == infinity:
                continue
            for j from 0 <= j < N:
                sum_ijk = graph[i, k] + graph[k, j]
                if sum_ijk < graph[i, j]:
                    graph[i, j] = sum_ijk

    graph[np.where(np.isinf(graph))] = 0

    return graph


@cython.boundscheck(False)
cdef np.ndarray dijkstra(dist_matrix,
                         np.ndarray[DTYPE_t, ndim=2] graph,
                         int directed=0):
    """
    Dijkstra algorithm using Fibonacci Heaps

    Parameters
    ----------
    graph : array or sparse matrix
        dist_matrix is the matrix of distances between connected points.
        unconnected points have distance=0.  It will be converted to
        a csr_matrix internally
    indptr :
        These arrays encode a distance matrix in compressed-sparse-row
        format.
    graph : ndarray
        on input, graph is the matrix of distances between connected points.
        unconnected points have distance=0
        on exit, graph is overwritten with the output
    directed : bool, default = False
        if True, then the algorithm will only traverse from a point to
        its neighbors when finding the shortest path.
        if False, then the algorithm will traverse all paths in both
        directions.

    Returns
    -------
    graph : array
        the matrix of shortest paths between points.
        If no path exists, the path length is zero
    """
    cdef unsigned int N = graph.shape[0]
    cdef unsigned int i

    cdef FibonacciHeap heap

    cdef FibonacciNode* nodes = <FibonacciNode*> malloc(N *
                                                        sizeof(FibonacciNode))

    cdef np.ndarray distances, neighbors, indptr
    cdef np.ndarray distances2, neighbors2, indptr2

    if not isspmatrix_csr(dist_matrix):
        dist_matrix = csr_matrix(dist_matrix)

    distances = np.asarray(dist_matrix.data, dtype=DTYPE, order='C')
    neighbors = np.asarray(dist_matrix.indices, dtype=ITYPE, order='C')
    indptr = np.asarray(dist_matrix.indptr, dtype=ITYPE, order='C')

    for i from 0 <= i < N:
        initialize_node(&nodes[i], i)

    heap.min_node = NULL

    if directed:
        for i from 0 <= i < N:
            dijkstra_directed_one_row(i, neighbors, distances, indptr,
                                      graph, &heap, nodes)
    else:
        #use the csr -> csc sparse matrix conversion to quickly get
        # both directions of neigbors
        dist_matrix_T = dist_matrix.T.tocsr()

        distances2 = np.asarray(dist_matrix_T.data,
                                dtype=DTYPE, order='C')
        neighbors2 = np.asarray(dist_matrix_T.indices,
                                dtype=ITYPE, order='C')
        indptr2 = np.asarray(dist_matrix_T.indptr,
                             dtype=ITYPE, order='C')

        for i from 0 <= i < N:
            dijkstra_one_row(i, neighbors, distances, indptr,
                             neighbors2, distances2, indptr2,
                             graph, &heap, nodes)

    free(nodes)

    return graph


######################################################################
# Debugging: Functions for printing the fibonacci heap
#
#cdef void print_node(FibonacciNode* node, int level=0):
#    print '%s(%i,%i) %i' % (level*'   ', node.index, node.val, node.rank)
#    if node.children:
#        print_node(leftmost_sibling(node.children), level+1)
#    if node.right_sibling:
#        print_node(node.right_sibling, level)
#
#
#cdef void print_heap(FibonacciHeap* heap):
#    print "---------------------------------"
#    if heap.min_node:
#        print_node(leftmost_sibling(heap.min_node))
#    else:
#        print "[empty heap]"


@cython.boundscheck(False)
cdef void dijkstra_directed_one_row(
                    unsigned int i_node,
                    np.ndarray[ITYPE_t, ndim=1, mode='c'] neighbors,
                    np.ndarray[DTYPE_t, ndim=1, mode='c'] distances,
                    np.ndarray[ITYPE_t, ndim=1, mode='c'] indptr,
                    np.ndarray[DTYPE_t, ndim=2, mode='c'] graph,
                    FibonacciHeap* heap,
                    FibonacciNode* nodes):
    """
    Calculate distances from a single point to all targets using a
    directed graph.

    Parameters
    ----------
    i_node : index of source point
    neighbors : array, shape = [N,]
        indices of neighbors for each point
    distances : array, shape = [N,]
        lengths of edges to each neighbor
    indptr : array, shape = (N+1,)
        the neighbors of point i are given by
        neighbors[indptr[i]:indptr[i+1]]
    graph : array, shape = (N,N)
        on return, graph[i_node] contains the path lengths from
        i_node to each target
    heap : the Fibonacci heap object to use
    nodes : the array of nodes to use
    """
    cdef unsigned int N = graph.shape[0]
    cdef unsigned int i
    cdef FibonacciNode *v
    cdef FibonacciNode *current_neighbor
    cdef DTYPE_t dist

    # initialize nodes
    for i from 0 <= i < N:
        initialize_node(&nodes[i], i)

    heap.min_node = NULL
    insert_node(heap, &nodes[i_node])

    while heap.min_node:
        v = remove_min(heap)
        v.state = 2  # 2 -> SCANNED

        for i from indptr[v.index] <= i < indptr[v.index + 1]:
            current_neighbor = &nodes[neighbors[i]]
            if current_neighbor.state != 2:      # 2 -> SCANNED
                dist = distances[i]
                if current_neighbor.state == 0:  # 0 -> NOT_IN_HEAP
                    current_neighbor.state = 1   # 1 -> IN_HEAP
                    current_neighbor.val = v.val + dist
                    insert_node(heap, current_neighbor)
                elif current_neighbor.val > v.val + dist:
                    decrease_val(heap, current_neighbor,
                                 v.val + dist)

        #v has now been scanned: add the distance to the results
        graph[i_node, v.index] = v.val


@cython.boundscheck(False)
cdef void dijkstra_one_row(unsigned int i_node,
                    np.ndarray[ITYPE_t, ndim=1, mode='c'] neighbors1,
                    np.ndarray[DTYPE_t, ndim=1, mode='c'] distances1,
                    np.ndarray[ITYPE_t, ndim=1, mode='c'] indptr1,
                    np.ndarray[ITYPE_t, ndim=1, mode='c'] neighbors2,
                    np.ndarray[DTYPE_t, ndim=1, mode='c'] distances2,
                    np.ndarray[ITYPE_t, ndim=1, mode='c'] indptr2,
                    np.ndarray[DTYPE_t, ndim=2, mode='c'] graph,
                    FibonacciHeap* heap,
                    FibonacciNode* nodes):
    """
    Calculate distances from a single point to all targets using an
    undirected graph.

    Parameters
    ----------
    i_node : index of source point
    neighbors[1,2] : array, shape = [N,]
        indices of neighbors for each point
    distances[1,2] : array, shape = [N,]
        lengths of edges to each neighbor
    indptr[1,2] : array, shape = (N+1,)
        the neighbors of point i are given by
        neighbors1[indptr1[i]:indptr1[i+1]] and
        neighbors2[indptr2[i]:indptr2[i+1]]
    graph : array, shape = (N,)
        on return, graph[i_node] contains the path lengths from
        i_node to each target
    heap : the Fibonacci heap object to use
    nodes : the array of nodes to use
    """
    cdef unsigned int N = graph.shape[0]
    cdef unsigned int i
    cdef FibonacciNode *v
    cdef FibonacciNode *current_neighbor
    cdef DTYPE_t dist

    # re-initialize nodes
    # children, parent, left_sibling, right_sibling should already be NULL
    # rank should already be 0, index will already be set
    # we just need to re-set state and val
    for i from 0 <= i < N:
        nodes[i].state = 0  # 0 -> NOT_IN_HEAP
        nodes[i].val = 0

    insert_node(heap, &nodes[i_node])

    while heap.min_node:
        v = remove_min(heap)
        v.state = 2  # 2 -> SCANNED

        for i from indptr1[v.index] <= i < indptr1[v.index + 1]:
            current_neighbor = &nodes[neighbors1[i]]
            if current_neighbor.state != 2:      # 2 -> SCANNED
                dist = distances1[i]
                if current_neighbor.state == 0:  # 0 -> NOT_IN_HEAP
                    current_neighbor.state = 1   # 1 -> IN_HEAP
                    current_neighbor.val = v.val + dist
                    insert_node(heap, current_neighbor)
                elif current_neighbor.val > v.val + dist:
                    decrease_val(heap, current_neighbor,
                                 v.val + dist)

        for i from indptr2[v.index] <= i < indptr2[v.index + 1]:
            current_neighbor = &nodes[neighbors2[i]]
            if current_neighbor.state != 2:      # 2 -> SCANNED
                dist = distances2[i]
                if current_neighbor.state == 0:  # 0 -> NOT_IN_HEAP
                    current_neighbor.state = 1   # 1 -> IN_HEAP
                    current_neighbor.val = v.val + dist
                    insert_node(heap, current_neighbor)
                elif current_neighbor.val > v.val + dist:
                    decrease_val(heap, current_neighbor,
                                 v.val + dist)

        #v has now been scanned: add the distance to the results
        graph[i_node, v.index] = v.val
