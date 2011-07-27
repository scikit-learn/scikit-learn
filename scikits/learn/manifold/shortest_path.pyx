"""
Routines for performing shortest-path graph searches

The main interface is in the function `shortest_path`.  This
calls cython routines that compute the shortest path using either
the Floyd-Warshall algorithm, or Dykstra's algorithm with Fibonacci Heaps.
"""

# Author: Jake Vanderplas  -- <vanderplas@astro.washington.edu>
# License: BSD, (C) 2011

import numpy as np
cimport numpy as np

from scipy.sparse import csr_matrix

cimport cython

from libc.stdlib cimport malloc, free

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

ITYPE = np.int32
ctypedef np.int32_t ITYPE_t

cdef inline DTYPE_t fmin(DTYPE_t a, DTYPE_t b):
   if a <= b:
       return a
   else:
       return b

def shortest_path(neighbors, distances, directed=False, method='best'):
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
    directed : boolean
        if True, then find the shortest path on a directed graph: only
        progress from a point to its neighbors, not the other way around.
        if False, then find the shortest path on an undirected graph: the
        algorithm can progress from a point to its neighbors and vice versa.
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
        if k < N/4:
            method = 'D'
        else:
            method = 'FW'

    if method == 'FW':
        if directed:
            FloydWarshallDirected(neighbors, distances, graph)
        else:
            FloydWarshall(neighbors, distances, graph)
    elif method == 'D':
        if directed:
            DijkstraDirected(neighbors, distances, graph)
        else:
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
cdef void FloydWarshallDirected(np.ndarray[ITYPE_t, ndim=2, mode='c'] neighbors,
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

    #now perform the Floyd-Warshall algorithm
    for k from 0 <= k < N:
        for i from 0 <= i < N:
            if graph[i,k] == infinity:
                continue
            for j from 0 <= j < N:
                sum_ijk = graph[i, k] + graph[k, j]
                if sum_ijk < graph[i,j]:
                    graph[i,j] = sum_ijk


######################################################################
# FibonacciNode structure
#  This structure and the operations on it are the nodes of the
#  Fibonacci heap.

cdef struct FibonacciNode:
    unsigned int index
    unsigned int rank
    unsigned int state
    DTYPE_t val
    FibonacciNode* parent
    FibonacciNode* left_sibling
    FibonacciNode* right_sibling
    FibonacciNode* children


cdef FibonacciNode* initialize_node(FibonacciNode* node,
                                    unsigned int index,
                                    DTYPE_t val=0):
     global UNLABELED

     node.index = index
     node.val = val
     node.rank = 0
     node.state = 0
     
     node.parent = NULL
     node.left_sibling = NULL
     node.right_sibling = NULL
     node.children = NULL
     
     return node
     

cdef FibonacciNode* rightmost_sibling(FibonacciNode* node):
    cdef FibonacciNode* temp = node
    while(temp.right_sibling):
        temp = temp.right_sibling
    return temp


cdef FibonacciNode* leftmost_sibling(FibonacciNode* node):
    cdef FibonacciNode* temp = node
    while(temp.left_sibling):
        temp = temp.left_sibling
    return temp


cdef FibonacciNode* add_child(FibonacciNode* node, FibonacciNode* child):
    child.right_sibling = NULL
    child.parent = node

    if node.children:
        add_sibling(node.children, child)
    else:
        node.children = child
        child.left_sibling = NULL
        node.rank = 1

    return child


cdef FibonacciNode* add_sibling(FibonacciNode* node, FibonacciNode* sibling):
    cdef FibonacciNode* temp = rightmost_sibling(node)
    temp.right_sibling = sibling
    sibling.left_sibling = temp
    sibling.right_sibling = NULL
    sibling.parent = node.parent
    if sibling.parent:
        sibling.parent.rank += 1

    return sibling


cdef FibonacciNode* remove(FibonacciNode* node):
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

    return node


######################################################################
# FibonacciHeap structure
#  This structure and operations on it use the FibonacciNode 
#  routines to implement a Fibonacci heap

ctypedef FibonacciNode* pFibonacciNode


cdef struct FibonacciHeap:
    FibonacciNode* min_node
    pFibonacciNode[100] roots_by_rank


cdef FibonacciHeap* initialize_heap(FibonacciHeap* heap):
    heap.min_node = NULL
    cdef unsigned int i
    for i from 0 <= i < 100:
        heap.roots_by_rank[i] = NULL


cdef FibonacciNode* insert_node(FibonacciHeap* heap,
                                FibonacciNode* node):
    if heap.min_node:
        add_sibling(heap.min_node, node)
        if node.val < heap.min_node.val:
            heap.min_node = node
    else:
        heap.min_node = node
    
    return node


cdef FibonacciNode* decrease_val(FibonacciHeap* heap,
                                 FibonacciNode* node,
                                 DTYPE_t newval):
    node.val = newval
    if node.parent:
        if node.parent.val >= newval:
            remove(node)
            add_sibling(heap.min_node, node)
            if node.val < heap.min_node.val:
                heap.min_node = node
    return node


cdef FibonacciNode* find_min(FibonacciHeap* heap):
    return heap.min_node


cdef FibonacciNode* link(FibonacciHeap* heap, FibonacciNode* node):
    cdef FibonacciNode *linknode, *tmp_parent, *tmp_child
    if heap.roots_by_rank[node.rank] == NULL:
        heap.roots_by_rank[node.rank] = node
    else:
        linknode = heap.roots_by_rank[node.rank]
        heap.roots_by_rank[node.rank] = NULL
        if node.val < linknode.val:
            tmp_parent = node
            tmp_child = linknode
        else:
            tmp_parent = linknode
            tmp_child = node
        
        remove(tmp_child)
        add_child(tmp_parent, tmp_child)
        if heap.roots_by_rank[tmp_parent.rank]:
            link(heap, tmp_parent)
        else:
            heap.roots_by_rank[tmp_parent.rank] = tmp_parent

    return node


cdef FibonacciNode* delete_min(FibonacciHeap* heap):
    cdef FibonacciNode *temp, *next_temp, *out
    cdef unsigned int i
    
    if heap.min_node == NULL:
        return NULL
    
    if heap.min_node.children:
        temp = leftmost_sibling(heap.min_node.children)
        next_temp = NULL

        while temp:
            next_temp = temp.right_sibling
            remove(temp)
            add_sibling(heap.min_node, temp)
            temp = next_temp
        
    temp = leftmost_sibling(heap.min_node)
    
    if temp == heap.min_node:
        if heap.min_node.right_sibling:
            temp = heap.min_node.right_sibling
        else:
            out = heap.min_node
            heap.min_node = NULL
            return out

    out = heap.min_node
    remove(heap.min_node)
    heap.min_node = temp
    
    for i from 0 <= i < 100:
        heap.roots_by_rank[i] = NULL

    while temp:
        if temp.val < heap.min_node.val:
            heap.min_node = temp
                
        next_temp = temp.right_sibling
        link(heap, temp)
        temp = next_temp
        
    return out


######################################################################
# Debugging: Functions for printing the fibonacci heap

cdef void print_node(FibonacciNode* node, int level=0):
    print '%s(%i,%i) %i' % (level*'   ', node.index, node.val, node.rank)

    if node.children:
        print_node(leftmost_sibling(node.children), level+1)
    
    if node.right_sibling:
        print_node(node.right_sibling, level)


cdef void print_heap(FibonacciHeap* heap):
    print "---------------------------------"
    if heap.min_node:
        print_node(leftmost_sibling(heap.min_node))
    else:
        print "[empty heap]"


@cython.boundscheck(False)
cdef void DijkstraDirectedOneRow(
    unsigned int i_node,
    np.ndarray[ITYPE_t, ndim=2, mode='c'] neighbors,
    np.ndarray[DTYPE_t, ndim=2, mode='c'] distances,
    np.ndarray[DTYPE_t, ndim=2, mode='c'] graph,
    FibonacciHeap* heap,
    FibonacciNode* nodes):
    """
    Calculate distances from a single point to all targets using a
    directed graph.

    Parameters
    ----------
    i_node : index of source point
    neighbors : array, shape = [N,k]
        indices of neighbors for each point
    distances : array, shape = [N,k]
        lengths of edges to each neighbor
    graph : array, shape = (N,N)
        on return, graph[i_node] contains the path lengths from 
        i_node to each target
    heap: the Fibonacci heap object to use
    nodes : the array of nodes to use
    """
    cdef int UNLABELED = 0
    cdef int LABELED = 1
    cdef int SCANNED = 2

    cdef unsigned int N = neighbors.shape[0]
    cdef unsigned int k = neighbors.shape[1]
    cdef unsigned int i
    cdef FibonacciNode *v, *current_neighbor
    cdef DTYPE_t dist

    # initialize nodes
    for i from 0 <= i < N:
        nodes[i].state = UNLABELED
        nodes[i].val = 0

    insert_node(heap, &nodes[i_node])

    while True:
        v = delete_min(heap)
        v.state = SCANNED

        for i from 0 <= i < k:
            current_neighbor = &nodes[neighbors[v.index, i]]
            dist = distances[v.index, i]
            if current_neighbor.state != SCANNED:
                if current_neighbor.state == UNLABELED:
                    current_neighbor.state = LABELED
                    current_neighbor.val = v.val + dist
                    insert_node(heap, current_neighbor)
                elif current_neighbor.val > v.val + dist:
                    decrease_val(heap, current_neighbor,
                                 v.val + dist)
        
        if heap.min_node == NULL:
            break

    for i from 0 <= i < N:
        graph[i_node, i] = nodes[i].val
    
        
@cython.boundscheck(False)
cdef void DijkstraDirected(np.ndarray[ITYPE_t, ndim=2, mode='c'] neighbors,
                           np.ndarray[DTYPE_t, ndim=2, mode='c'] distances,
                           np.ndarray[DTYPE_t, ndim=2, mode='c'] graph):
    cdef unsigned int N = neighbors.shape[0]
    cdef unsigned int k = neighbors.shape[1]
    cdef unsigned int i

    cdef FibonacciHeap heap

    cdef FibonacciNode* nodes = \
        <FibonacciNode*> malloc(N * sizeof(FibonacciNode))

    for i from 0 <= i < N:
        initialize_node(&(nodes[i]), i)

    initialize_heap(&heap)

    for i from 0 <= i < N:
        DijkstraDirectedOneRow(i,
                               neighbors,
                               distances,
                               graph,
                               &heap,
                               nodes)


@cython.boundscheck(False)
cdef void DijkstraOneRow(
    unsigned int i_node,
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
    Calculate distances from a single point to all targets using a
    directed graph.

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
    cdef int UNLABELED = 0
    cdef int LABELED = 1
    cdef int SCANNED = 2

    cdef unsigned int N = graph.shape[0]
    cdef unsigned int i
    cdef FibonacciNode *v, *current_neighbor
    cdef DTYPE_t dist

    # initialize nodes
    for i from 0 <= i < N:
        nodes[i].state = UNLABELED
        nodes[i].val = 0

    insert_node(heap, &nodes[i_node])

    while True:
        v = delete_min(heap)
        v.state = SCANNED

        for i from indptr1[v.index] <= i < indptr1[v.index + 1]:
            current_neighbor = &nodes[neighbors1[i]]
            dist = distances1[i]
            if current_neighbor.state != SCANNED:
                if current_neighbor.state == UNLABELED:
                    current_neighbor.state = LABELED
                    current_neighbor.val = v.val + dist
                    insert_node(heap, current_neighbor)
                elif current_neighbor.val > v.val + dist:
                    decrease_val(heap, current_neighbor,
                                 v.val + dist)

        for i from indptr2[v.index] <= i < indptr2[v.index + 1]:
            current_neighbor = &nodes[neighbors2[i]]
            dist = distances2[i]
            if current_neighbor.state != SCANNED:
                if current_neighbor.state == UNLABELED:
                    current_neighbor.state = LABELED
                    current_neighbor.val = v.val + dist
                    insert_node(heap, current_neighbor)
                elif current_neighbor.val > v.val + dist:
                    decrease_val(heap, current_neighbor,
                                 v.val + dist)
        
        if heap.min_node == NULL:
            break

    for i from 0 <= i < N:
        graph[i_node, i] = nodes[i].val
    

@cython.boundscheck(False)
cdef void Dijkstra(np.ndarray[ITYPE_t, ndim=2] neighbors,
                   np.ndarray[DTYPE_t, ndim=2] distances,
                   np.ndarray[DTYPE_t, ndim=2] graph):
    cdef unsigned int N = neighbors.shape[0]
    cdef unsigned int k = neighbors.shape[1]
    cdef unsigned int i

    cdef np.ndarray distances1 = distances.reshape(distances.size)
    cdef np.ndarray neighbors1 = neighbors.reshape(neighbors.size)
    cdef np.ndarray indptr1 = np.arange(0, N * k + 1, k, dtype=ITYPE)

    #use the csr -> csc sparse matrix conversion to quickly get
    # both directions of neigbors
    dist_matrix = csr_matrix((distances1, neighbors1, indptr1))
    dist_matrix_T = dist_matrix.T.tocsr()

    cdef np.ndarray distances2 = np.asarray(dist_matrix_T.data,
                                            dtype=DTYPE)
    cdef np.ndarray neighbors2 = np.asarray(dist_matrix_T.indices,
                                            dtype=ITYPE)
    cdef np.ndarray indptr2 = np.asarray(dist_matrix_T.indptr,
                                         dtype=ITYPE)

    cdef FibonacciHeap heap

    cdef FibonacciNode* nodes = \
        <FibonacciNode*> malloc(N * sizeof(FibonacciNode))

    for i from 0 <= i < N:
        initialize_node(&(nodes[i]), i)

    initialize_heap(&heap)

    for i from 0 <= i < N:
        DijkstraOneRow(i,
                       neighbors1, distances1, indptr1,
                       neighbors2, distances2, indptr2,
                       graph,
                       &heap,
                       nodes)
