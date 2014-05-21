"""
Routines for performing shortest-path graph searches

The main interface is in the function `graph_shortest_path`.  This
calls cython routines that compute the shortest path using either
the Floyd-Warshall algorithm, or Dykstra's algorithm with Fibonacci Heaps.
"""

# Author: Jake Vanderplas  -- <vanderplas@astro.washington.edu>
# License: BSD 3 clause, (C) 2011
# Modified by Peter Fischer -- <peter.fischer@fau.de>
# - also return shortest paths (predecessor matrix)
# - modified dijkstra for partial shortest paths

import numpy as np
cimport numpy as np
cdef extern from "numpy/npy_math.h":
    double NPY_INFINITY

from scipy.sparse import csr_matrix, isspmatrix, isspmatrix_csr

cimport cython

from libc.stdlib cimport malloc, free

np.import_array()

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

ITYPE = np.int32
ctypedef np.int32_t ITYPE_t


def graph_shortest_path(dist_matrix not None, directed=True, method='auto'):
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

    P : np.ndarray, int, shape = [N,N]
        Predecessor matrix which stores the shortest paths.
        P[i,j] gives the index k of the predecessor of node j on the shortest
        path from node i to node j

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

    Addition (peterf):
    Dijkstra is not optimally implemented for undirected graphs. If G[i,j] is
    known, it can also be used for G[j,i], which is not done here. This would
    reduce the number of computations by a factor of 2.
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

    predecessors = np.zeros(dist_matrix.shape, dtype=ITYPE, order='C')
    if method == 'FW':
        graph = np.asarray(dist_matrix.toarray(), dtype=DTYPE, order='C')
        floyd_warshall(graph, predecessors, directed)
    elif method == 'D':
        graph = np.zeros((N, N), dtype=DTYPE, order='C')
        dijkstra(dist_matrix, graph, predecessors, directed)
    else:
        raise ValueError("unrecognized method '%s'" % method)

    return graph, predecessors


@cython.boundscheck(False)
cdef void floyd_warshall(DTYPE_t[:, ::1] graph,
                         ITYPE_t[:, ::1] predecessors,
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
    cdef unsigned int N = graph.shape[0]
    assert graph.shape[1] == N

    cdef unsigned int i, j, k, m

    cdef DTYPE_t sum_ijk

    #initialize all distances to infinity
    for i in xrange(N):
        for j in xrange(N):
            if graph[i, j] == 0:
                graph[i, j] = NPY_INFINITY

    #graph[i,i] should be zero
    for i in xrange(N):
        graph[i, i] = 0

    # for a non-directed graph, we need to symmetrize the distances
    if not directed:
        for i from 0 <= i < N:
            for j from i + 1 <= j < N:
                if graph[j, i] <= graph[i, j]:
                    graph[i, j] = graph[j, i]
                else:
                    graph[j, i] = graph[i, j]

    # initialize predecessors
    for i in xrange(N):
        for j in xrange(N):
            if graph[i, j] != 0 and graph[i, j] != NPY_INFINITY:
                predecessors[i, j] = i
            else:
                predecessors[i, j] = -1

    #now perform the Floyd-Warshall algorithm
    for k from 0 <= k < N:
        for i from 0 <= i < N:
            if graph[i, k] == NPY_INFINITY:
                continue
            for j from 0 <= j < N:
                sum_ijk = graph[i, k] + graph[k, j]
                if sum_ijk < graph[i, j]:
                    graph[i, j] = sum_ijk
                    predecessors[i, j] = predecessors[k, j]

    # Set unknown distances to 0 instead of inf
    for i in xrange(N):
        for j in xrange(N):
            if graph[i, j] == NPY_INFINITY:
                graph[i, j] = 0

    return


@cython.boundscheck(False)
cdef void dijkstra(dist_matrix,
                   DTYPE_t[:, ::1] graph,
                   ITYPE_t[:, ::1] predecessors,
                   int directed=0):
    """
    Dijkstra algorithm using Fibonacci Heaps

    Parameters
    ----------
    graph : array or sparse matrix
        dist_matrix is the matrix of distances betweeen connected points.
        unconnected points have distance=0.  It will be converted to
        a csr_matrix internally
    indptr :
        These arrays encode a distance matrix in compressed-sparse-row
        format.
    graph : ndarray
        on input, graph is the matrix of distances betweeen connected points.
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

    cdef DTYPE_t[::1] distances, distances2
    cdef ITYPE_t[::1] neighbors, indptr, neighbors2, indptr2
    cdef ITYPE_t[::1] predecessor_buffer

    if not isspmatrix_csr(dist_matrix):
        dist_matrix = csr_matrix(dist_matrix)

    distances = dist_matrix.data
    neighbors = dist_matrix.indices
    indptr = dist_matrix.indptr
    predecessor_buffer = np.zeros(N, dtype=ITYPE, order='C')

    for i from 0 <= i < N:
        initialize_node(&nodes[i], i)

    heap.min_node = NULL

    if directed:
        for i from 0 <= i < N:
            predecessor_buffer[:] = -1
            dijkstra_directed_one_row(i, neighbors, distances, indptr,
                                      graph, &heap, nodes, predecessor_buffer)
            predecessors[i, :] = predecessor_buffer
    else:
        #use the csr -> csc sparse matrix conversion to quickly get
        # both directions of neigbors
        # Assumption is: dist_matrix is NOT symmetric
        dist_matrix_T = dist_matrix.T.tocsr()

        distances2 = dist_matrix_T.data
        neighbors2 = dist_matrix_T.indices
        indptr2 = dist_matrix_T.indptr

        for i from 0 <= i < N:
            predecessor_buffer[:] = -1
            dijkstra_one_row(i, neighbors, distances, indptr,
                             neighbors2, distances2, indptr2,
                             graph, &heap, nodes, predecessor_buffer)
            predecessors[i, :] = predecessor_buffer

    free(nodes)

    return


######################################################################
# FibonacciNode structure
#  This structure and the operations on it are the nodes of the
#  Fibonacci heap.
#

cdef struct FibonacciNode:
    unsigned int index
    unsigned int rank
    unsigned int state
    DTYPE_t val
    FibonacciNode* parent
    FibonacciNode* left_sibling
    FibonacciNode* right_sibling
    FibonacciNode* children


cdef void initialize_node(FibonacciNode* node,
                          unsigned int index,
                          DTYPE_t val=0):
    # Assumptions: - node is a valid pointer
    #              - node is not currently part of a heap
    node.index = index
    node.val = val
    node.rank = 0
    node.state = 0  # 0 -> NOT_IN_HEAP

    node.parent = NULL
    node.left_sibling = NULL
    node.right_sibling = NULL
    node.children = NULL


cdef FibonacciNode* rightmost_sibling(FibonacciNode* node):
    # Assumptions: - node is a valid pointer
    cdef FibonacciNode* temp = node
    while(temp.right_sibling):
        temp = temp.right_sibling
    return temp


cdef FibonacciNode* leftmost_sibling(FibonacciNode* node):
    # Assumptions: - node is a valid pointer
    cdef FibonacciNode* temp = node
    while(temp.left_sibling):
        temp = temp.left_sibling
    return temp


cdef void add_child(FibonacciNode* node, FibonacciNode* new_child):
    # Assumptions: - node is a valid pointer
    #              - new_child is a valid pointer
    #              - new_child is not the sibling or child of another node
    new_child.parent = node

    if node.children:
        add_sibling(node.children, new_child)
    else:
        node.children = new_child
        new_child.right_sibling = NULL
        new_child.left_sibling = NULL
        node.rank = 1


cdef void add_sibling(FibonacciNode* node, FibonacciNode* new_sibling):
    # Assumptions: - node is a valid pointer
    #              - new_sibling is a valid pointer
    #              - new_sibling is not the child or sibling of another node
    cdef FibonacciNode* temp = rightmost_sibling(node)
    temp.right_sibling = new_sibling
    new_sibling.left_sibling = temp
    new_sibling.right_sibling = NULL
    new_sibling.parent = node.parent
    if new_sibling.parent:
        new_sibling.parent.rank += 1


cdef void remove(FibonacciNode* node):
    # Assumptions: - node is a valid pointer
    if node.parent:
        node.parent.rank -= 1
        if node.left_sibling:
            node.parent.children = node.left_sibling
        elif node.right_sibling:
            node.parent.children = node.right_sibling
        else:
            node.parent.children = NULL

    if node.left_sibling:
        node.left_sibling.right_sibling = node.right_sibling
    if node.right_sibling:
        node.right_sibling.left_sibling = node.left_sibling

    node.left_sibling = NULL
    node.right_sibling = NULL
    node.parent = NULL


######################################################################
# FibonacciHeap structure
#  This structure and operations on it use the FibonacciNode
#  routines to implement a Fibonacci heap

ctypedef FibonacciNode* pFibonacciNode


cdef struct FibonacciHeap:
    FibonacciNode* min_node
    pFibonacciNode[100] roots_by_rank  # maximum number of nodes is ~2^100.


cdef void insert_node(FibonacciHeap* heap,
                      FibonacciNode* node):
    # Assumptions: - heap is a valid pointer
    #              - node is a valid pointer
    #              - node is not the child or sibling of another node
    if heap.min_node:
        add_sibling(heap.min_node, node)
        if node.val < heap.min_node.val:
            heap.min_node = node
    else:
        heap.min_node = node


cdef void decrease_val(FibonacciHeap* heap,
                       FibonacciNode* node,
                       DTYPE_t newval):
    # Assumptions: - heap is a valid pointer
    #              - newval <= node.val
    #              - node is a valid pointer
    #              - node is not the child or sibling of another node
    node.val = newval
    if node.parent and (node.parent.val >= newval):
        remove(node)
        insert_node(heap, node)
    elif heap.min_node.val > node.val:
        heap.min_node = node


cdef void link(FibonacciHeap* heap, FibonacciNode* node):
    # Assumptions: - heap is a valid pointer
    #              - node is a valid pointer
    #              - node is already within heap

    cdef FibonacciNode *linknode, *parent, *child

    if heap.roots_by_rank[node.rank] == NULL:
        heap.roots_by_rank[node.rank] = node
    else:
        linknode = heap.roots_by_rank[node.rank]
        heap.roots_by_rank[node.rank] = NULL

        if node.val < linknode.val or node == heap.min_node:
            remove(linknode)
            add_child(node, linknode)
            link(heap, node)
        else:
            remove(node)
            add_child(linknode, node)
            link(heap, linknode)


cdef FibonacciNode* remove_min(FibonacciHeap* heap):
    # Assumptions: - heap is a valid pointer
    #              - heap.min_node is a valid pointer
    cdef FibonacciNode *temp, *temp_right, *out
    cdef unsigned int i

    # make all min_node children into root nodes
    if heap.min_node.children:
        temp = leftmost_sibling(heap.min_node.children)
        temp_right = NULL

        while temp:
            temp_right = temp.right_sibling
            remove(temp)
            add_sibling(heap.min_node, temp)
            temp = temp_right

        heap.min_node.children = NULL

    # choose a root node other than min_node
    temp = leftmost_sibling(heap.min_node)
    if temp == heap.min_node:
        if heap.min_node.right_sibling:
            temp = heap.min_node.right_sibling
        else:
            out = heap.min_node
            heap.min_node = NULL
            return out

    # remove min_node, and point heap to the new min
    out = heap.min_node
    remove(heap.min_node)
    heap.min_node = temp

    # re-link the heap
    for i from 0 <= i < 100:
        heap.roots_by_rank[i] = NULL

    while temp:
        if temp.val < heap.min_node.val:
            heap.min_node = temp
        temp_right = temp.right_sibling
        link(heap, temp)
        temp = temp_right

    return out


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


cdef void dijkstra_directed_one_row(
                    unsigned int i_node,
                    ITYPE_t[::1] neighbors,
                    DTYPE_t[::1] distances,
                    ITYPE_t[::1] indptr,
                    DTYPE_t[:, ::1] graph,
                    FibonacciHeap* heap,
                    FibonacciNode* nodes,
                    ITYPE_t[::1] predecessors):
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
    heap: the Fibonacci heap object to use
    nodes : the array of nodes to use
    """
    cdef unsigned int N = graph.shape[0]
    cdef unsigned int i
    cdef FibonacciNode *v, *current_neighbor
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
                    predecessors[current_neighbor.index] = v.index
                elif current_neighbor.val > v.val + dist:
                    decrease_val(heap, current_neighbor,
                                 v.val + dist)
                    predecessors[current_neighbor.index] = v.index

        #v has now been scanned: add the distance to the results
        graph[i_node, v.index] = v.val
    return


cdef void dijkstra_one_row(unsigned int i_node,
                    ITYPE_t[::1] neighbors1,
                    DTYPE_t[::1] distances1,
                    ITYPE_t[::1] indptr1,
                    ITYPE_t[::1] neighbors2,
                    DTYPE_t[::1] distances2,
                    ITYPE_t[::1] indptr2,
                    DTYPE_t[:, ::1] graph,
                    FibonacciHeap* heap,
                    FibonacciNode* nodes,
                    ITYPE_t[::1] predecessors):
    """
    Calculate distances from a single point to all targets using an
    undirected graph.

    TODO:
    update both directions, G[i,j] and G[j,i]

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
    heap: the Fibonacci heap object to use
    nodes : the array of nodes to use
    """
    cdef unsigned int N = graph.shape[0]
    cdef unsigned int i
    cdef FibonacciNode *v, *current_neighbor
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
                    predecessors[current_neighbor.index] = v.index
                elif current_neighbor.val > v.val + dist:
                    decrease_val(heap, current_neighbor,
                                 v.val + dist)
                    predecessors[current_neighbor.index] = v.index

        for i from indptr2[v.index] <= i < indptr2[v.index + 1]:
            current_neighbor = &nodes[neighbors2[i]]
            if current_neighbor.state != 2:      # 2 -> SCANNED
                dist = distances2[i]
                if current_neighbor.state == 0:  # 0 -> NOT_IN_HEAP
                    current_neighbor.state = 1   # 1 -> IN_HEAP
                    current_neighbor.val = v.val + dist
                    insert_node(heap, current_neighbor)
                    predecessors[current_neighbor.index] = v.index
                elif current_neighbor.val > v.val + dist:
                    decrease_val(heap, current_neighbor,
                                 v.val + dist)
                    predecessors[current_neighbor.index] = v.index

        #v has now been scanned: add the distance to the results
        graph[i_node, v.index] = v.val
    return


def modified_graph_shortest_path(i_node,
                                 N,
                                 np.ndarray targets not None,
                                 kng not None,
                                 np.ndarray graph not None,
                                 np.ndarray predecessors not None,
                                 directed=True):
    """
    Perform a shortest-path graph search on a positive directed or
    undirected graph between a specific set of points. Only Dijkstra method allowed
 
    Parameters
    ----------
    kng : arraylike or sparse matrix, shape = (N,N)
        Array of positive distances.
        If vertex i is connected to vertex j, then kng[i,j] gives
        the distance between the vertices.
        If vertex i is not connected to vertex j, then kng[i,j] = 0
    directed : boolean
        if True, then find the shortest path on a directed graph: only
        progress from a point to its neighbors, not the other way around.
        if False, then find the shortest path on an undirected graph: the
        algorithm can progress from a point to its neighbors and vice versa.

    Returns
    -------
    G : np.ndarray, float, shape = [N,N]
        G[i,j] gives the shortest distance from point i to point j
        along the graph.

    P : np.ndarray, int, shape = [N,N]
        Predecessor matrix which stores the shortest paths.
        P[i,j] gives the index k of the predecessor of node j on the shortest
        path from node i to node j

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
    if not isspmatrix_csr(kng):
        neighbor_matrix = csr_matrix(kng)
    else:
        neighbor_matrix = kng

    modified_dijkstra(i_node, N, targets, neighbor_matrix, 
                      graph, predecessors, directed)

    return graph, predecessors


@cython.boundscheck(False)
cdef void modified_dijkstra(unsigned int i_node,
                         unsigned int N,
                         ITYPE_t[::1] targets,
                         kng,
                         DTYPE_t[:, ::1] graph,
                         ITYPE_t[:, ::1] predecessors,
                         int directed=0):
    """
    Dijkstra algorithm using Fibonacci Heaps

    Parameters
    ----------
    graph : array or sparse matrix
        dist_matrix is the matrix of distances betweeen connected points.
        unconnected points have distance=0.  It will be converted to
        a csr_matrix internally
    indptr :
        These arrays encode a distance matrix in compressed-sparse-row
        format.
    graph : ndarray
        on input, graph is the matrix of distances betweeen connected points.
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
    cdef unsigned int i
    cdef FibonacciHeap heap
    cdef FibonacciNode* nodes = <FibonacciNode*> malloc(N *
                                                        sizeof(FibonacciNode))
    cdef DTYPE_t[::1] distances
    cdef ITYPE_t[::1] neighbors, indptr

    distances = kng.data
    neighbors = kng.indices
    indptr = kng.indptr

    for i from 0 <= i < N:
        initialize_node(&nodes[i], i)
    heap.min_node = NULL

    if directed:
        modified_dijkstra_directed_one_row(i_node, N, targets, neighbors,
                                           distances, indptr, graph, &heap,
                                           nodes, predecessors[i_node, :])
    else:
        modified_dijkstra_one_row(i_node, N, targets, neighbors, distances,
                                  indptr, graph, &heap, nodes,
                                  predecessors[i_node, :],
                                  predecessors[:, i_node])

    free(nodes)
    return


@cython.boundscheck(False)
cdef void modified_dijkstra_directed_one_row(
                    unsigned int i_node,
                    unsigned int N,
                    ITYPE_t[::1] targets,
                    ITYPE_t[::1] neighbors,
                    DTYPE_t[::1] distances,
                    ITYPE_t[::1] indptr,
                    DTYPE_t[:, ::1] graph,
                    FibonacciHeap* heap,
                    FibonacciNode* nodes,
                    ITYPE_t[::1] predecessors):
    """
    Calculate distances from a single point to a specific set of targets using a
    directed graph.

    Parameters
    ----------
    i_node : index of source point
    targets : array, shape = [K,]
        indices targets to which shortest path is computed
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
    predecessors : array, shape = [N,]
        on return, predecessors contains the predecessor on the shortest
        path from i_node to each target
    heap: the Fibonacci heap object to use
    nodes : the array of nodes to use
    """
    cdef unsigned int K = targets.shape[0]
    cdef unsigned int i, j, k, l, m
    cdef FibonacciNode *v, *current_neighbor
    cdef DTYPE_t dist
    cdef DTYPE_t delta

    # initialize nodes
    for i from 0 <= i < N:
        initialize_node(&nodes[i], i)

    heap.min_node = NULL

    for i in range(K):  # avoid processing of targets
        j = targets[i]
        v = &nodes[j]
        v.state = 2  # should not be processed

    for i in range(K):
        j = targets[i]

        delta = NPY_INFINITY
        for l from indptr[j] <= l < indptr[j + 1]:
            k = neighbors[l]
            current_neighbor = &nodes[k]
            if current_neighbor.state != 2:
                dist = graph[i_node, k]
                for m from indptr[j] <= m < indptr[j + 1]:
                    if neighbors[m] == k:
                        dist += distances[m]
                if dist < delta:
                    delta = dist
                    predecessors[j] = k
        v = &nodes[j]
        v.val = delta
        insert_node(heap, &nodes[j])

    while heap.min_node:
        v = remove_min(heap)
        v.state = 3  # 3 -> DISTANCE VALID, can be used again
        k = v.index
        graph[i_node, k] = v.val

        for i from indptr[k] <= i < indptr[k + 1]:
            j = neighbors[i]
            current_neighbor = &nodes[j]
            if current_neighbor.state == 2:
                dist = graph[i_node, k]
                for m from indptr[j] <= m < indptr[j + 1]:
                    if neighbors[m] == k:
                        dist += distances[m]
                for m from indptr[k] <= m < indptr[k + 1]:
                    if neighbors[m] == j:
                        current_neighbor = &nodes[neighbors[m]]

                if dist < current_neighbor.val:
                    decrease_val(heap, current_neighbor, dist)
                    predecessors[current_neighbor.index] = k
    return


@cython.boundscheck(False)
cdef void modified_dijkstra_one_row(
                    unsigned int i_node,
                    unsigned int N,
                    ITYPE_t[::1] targets,
                    ITYPE_t[::1] neighbors1,
                    DTYPE_t[::1] distances1,
                    ITYPE_t[::1] indptr1,
                    DTYPE_t[:, ::1] graph,
                    FibonacciHeap* heap,
                    FibonacciNode* nodes,
                    ITYPE_t[::1] predecessors1,
                    ITYPE_t[:] predecessors2):
    """
    Calculate distances from a single point to a specific set of targets using
    an undirected graph.

    At the same time, calculates the inverse distances and predecessors
    (from targets to single point)
    Careful: this semantic is different from the unmodified versions!

    Parameters
    ----------
    i_node : index of source point
    targets : array, shape = [K,]
        indices targets to which shortest path is computed
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
    predecessors : array, shape = [N,]
        on return, predecessors contains the predecessor on the shortest
        path from i_node to each target
    heap: the Fibonacci heap object to use
    nodes : the array of nodes to use
    """
    cdef unsigned int K = targets.shape[0]
    cdef int i, j, k, l, m
    cdef FibonacciNode *v, *current_neighbor
    cdef DTYPE_t dist
    cdef DTYPE_t delta
 
    # initialize nodes
    for i from 0 <= i < N:
        initialize_node(&nodes[i], i)
    heap.min_node = NULL

    # avoid processing of targets
    for i in range(K):
        j = targets[i]
        v = &nodes[j]
        v.state =  2 # should not be processed

    for i in range(K):
        j = targets[i]
        delta = NPY_INFINITY
        for l from indptr1[j] <= l < indptr1[j + 1]:
            k = neighbors1[l]
            current_neighbor = &nodes[k]
            if current_neighbor.state != 2:
                dist = graph[i_node, k]
                for m from indptr1[j] <= m < indptr1[j + 1]:
                    if neighbors1[m] == k:
                        # if element was not removed from sparse matrix
                        if distances1[m] < 0.00001:
                            dist = NPY_INFINITY
                        else:
                            dist += distances1[m]
                # More efficient version of:
                # dist = graph[i_node, k] + distances1[indptr1[j]:indptr1[j+1]][neighbors1[indptr1[j]:indptr1[j+1]] == k][0]
                if dist < delta:
                    delta = dist
                    predecessors1[j] = k
                    # successor on path from i_node to j is the same as
                    # from i_node to k
                    # if i_node is equal to k, the successor is j itself
                    if i_node == k:
                        predecessors2[j] = j
                    else:
                        predecessors2[j] = predecessors2[k]

        v = &nodes[j]
        v.val = delta
        insert_node(heap, &nodes[j])

    while heap.min_node:
        v = remove_min(heap)
        v.state = 3     # 3 -> DISTANCE VALID, can be used again
        k = v.index
        graph[i_node, k] = v.val
        # in undirected dijkstra, assume distance in both directions is the same
        graph[k, i_node] = v.val

        for i from indptr1[k] <= i < indptr1[k + 1]:
            j = neighbors1[i]
            current_neighbor = &nodes[j]
            if current_neighbor.state == 2:
                # Should not assume that if j in N(k), then also k in N(j), where N:Neighborhood
                # What is assumed here: kng[i,j] should be identical to kng[j,i], if it is not 0
                dist = graph[i_node, k]
                for m from indptr1[k] <= m < indptr1[k + 1]:
                    if neighbors1[m] == j:
                        if distances1[m] < 0.00001:
                            dist = NPY_INFINITY
                        else:
                            dist += distances1[m]
                            current_neighbor = &nodes[neighbors1[m]]

                if dist < current_neighbor.val:
                    decrease_val(heap, current_neighbor, dist)
                    predecessors1[current_neighbor.index] = k
                    if i_node == k:
                        predecessors2[current_neighbor.index] = \
                                current_neighbor.index
                    else:
                        predecessors2[current_neighbor.index] = \
                                predecessors2[k]
    return
