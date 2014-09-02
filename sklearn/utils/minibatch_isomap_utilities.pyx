#cython: boundscheck=False
#cython: wraparound=False
#cython: cdivision=True
#cython: profile=False

"""
Author:  Peter Fischer -- <peter.fischer@fau.de>
License: BSD 3 clause, (C) 2014

Helper functions for minibatch isomap
Law and Jain. Incremental Nonlinear Dimensionality Reduction
by Manifold Learning. TPAMI 2006
"""

import numpy as np
cimport numpy as np
np.import_array()
cdef extern from "numpy/npy_math.h":
    double NPY_INFINITY

cimport cython

from scipy.sparse import csr_matrix, isspmatrix_csr

from libc.stdlib cimport malloc, free
from libcpp.vector cimport vector
from libc.math cimport sqrt

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

ITYPE = np.int_
ctypedef np.int_t ITYPE_t

######################################################################
# FibonacciNode structure
#  This structure and the operations on it are the nodes of the
#  Fibonacci heap.
#
cdef struct FibonacciNode:
    Py_ssize_t index
    Py_ssize_t rank
    Py_ssize_t state
    DTYPE_t val
    FibonacciNode* parent
    FibonacciNode* left_sibling
    FibonacciNode* right_sibling
    FibonacciNode* children


cdef void initialize_node(FibonacciNode* node,
                          Py_ssize_t index,
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
    cdef Py_ssize_t i

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


cdef void _construct_fab(ITYPE_t a, ITYPE_t b,
                         Py_ssize_t N,
                         np.ndarray[ITYPE_t, ndim=2, mode='c'] F,
                         np.ndarray[ITYPE_t, ndim=1, mode='c'] neighbors,
                         np.ndarray[DTYPE_t, ndim=1, mode='c'] distances,
                         np.ndarray[ITYPE_t, ndim=1, mode='c'] indptr,
                         np.ndarray[ITYPE_t, ndim=2, mode='c'] pred_,
                         FibonacciHeap* heap,
                         FibonacciNode* nodes):
    '''
    Implementation of ConstructFab from Law-Paper
    similar naming conventions.
    Changed to store F directly in auxiliary matrix 
    pi = self.predecessor_matrix_
    '''
    cdef Py_ssize_t i, j, t
    cdef FibonacciNode *v, *current_neighbor
    cdef ITYPE_t s, u, k

    for i from 0 <= i < N:
        initialize_node(&nodes[i], i)
    heap.min_node = NULL

    if pred_[a, b] == a:
        insert_node(heap, &nodes[a])
    cdef vector[ITYPE_t] Rab
    while heap.min_node:
        v = remove_min(heap)
        v.state = 2  # 2 -> SCANNED
        t = v.index
        Rab.push_back(t)
        for k from indptr[t] <= k < indptr[t + 1]:
            u = neighbors[k]
            current_neighbor = &nodes[u]
            if current_neighbor.state != 2:             # 2 -> SCANNED
                if pred_[u, b] == a and u != a:
                    if current_neighbor.state == 0:     # 0 -> NOT_IN_HEAP
                        current_neighbor.state = 1      # 1 -> IN_HEAP
                        insert_node(heap, current_neighbor)

    for j from 0 <= j < Rab.size():
        u = Rab[j]
        # empty the heap
        for i from 0 <= i < N:
            initialize_node(&nodes[i], i)
            heap.min_node = NULL
        insert_node(heap, &nodes[b])
        while heap.min_node:
            v = remove_min(heap)
            v.state = 2
            t = v.index
            if pred_[u, t] == pred_[a, t]:
                F[u, t] = 1
                F[t, u] = 1
                for k from indptr[t] <= k < indptr[t + 1]:
                    s = neighbors[k]
                    current_neighbor = &nodes[s]
                    if current_neighbor.state != 2:  # 2 -> SCANNED
                        if pred_[a, s] == t and s != b:
                            if current_neighbor.state == 0:  # 0 -> NOT_IN_HEAP
                                current_neighbor.state = 1   # 1 -> IN_HEAP
                                insert_node(heap, current_neighbor)
    return


def _construct_f(deleted_edges,
                 Py_ssize_t N,
                 kng not None,
                 np.ndarray pred_ not None):
    cdef ITYPE_t a, b

    # Initialize heap
    cdef FibonacciHeap heap
    cdef FibonacciNode* nodes = <FibonacciNode*> malloc(N *
                                                        sizeof(FibonacciNode))

    # Convert sparse matrix to np.array which can be handled efficiently
    cdef np.ndarray distances, neighbors, indptr
    if not isspmatrix_csr(kng):
        kng = csr_matrix(kng)
    distances = kng.data  # np.asarray(kng.indices, dtype=ITYPE, order='C')
    neighbors = kng.indices
    indptr = kng.indptr

    F = np.zeros((N, N), dtype=ITYPE)
    for (a, b) in deleted_edges:
        _construct_fab(a, b, N, F, neighbors, distances, indptr, pred_,
                       &heap, nodes)

    free(nodes)
    return F


cdef void _update_insert_ab(ITYPE_t a, ITYPE_t b, Py_ssize_t n,
                            Py_ssize_t N,
                     np.ndarray[DTYPE_t, ndim=2, mode='c'] dist_matrix_,
                     np.ndarray[ITYPE_t, ndim=2, mode='c'] pred_,
                     FibonacciHeap* heap,
                     FibonacciNode* nodes):
    cdef Py_ssize_t i, j, t
    cdef ITYPE_t s, u
    cdef FibonacciNode *v, *current_neighbor
    cdef DTYPE_t dist

    for i from 0 <= i < N:
        initialize_node(&nodes[i], i)
    heap.min_node = NULL
    insert_node(heap, &nodes[a])

    cdef vector[ITYPE_t] S
    while heap.min_node:
        v = remove_min(heap)
        v.state = 2  # -> scanned
        t = v.index
        S.push_back(t)
        for u in xrange(N):
            if pred_[u, n] == t:
                current_neighbor = &nodes[u]
                if current_neighbor.state != 2:             # 2 -> SCANNED
                    if (dist_matrix_[u, n] + dist_matrix_[n, b] <
                        dist_matrix_[u, b]) and u != t:
                        if current_neighbor.state == 0:     # 0 -> NOT_IN_HEAP
                            current_neighbor.state = 1      # 1 -> IN_HEAP
                            insert_node(heap, current_neighbor)

    for j from 0 <= j < S.size():
        u = S[j]
        # empty the heap
        for i from 0 <= i < N:
            initialize_node(&nodes[i], i)
        heap.min_node = NULL
        insert_node(heap, &nodes[b])
        while heap.min_node:
            v = remove_min(heap)
            v.state = 2
            t = v.index
            dist = dist_matrix_[u, n] + dist_matrix_[n, t]
            dist_matrix_[u, t] = dist
            dist_matrix_[t, u] = dist
            pred_[u, t] = pred_[n, t]
            pred_[t, u] = pred_[n, u]

            for s in xrange(N):
                if pred_[s, n] == t:
                    current_neighbor = &nodes[s]
                    if current_neighbor.state != 2:             # 2 -> SCANNED
                    # Removed bug in the pseudocode of the paper
                        if ((dist_matrix_[s, n] + dist_matrix_[n, a] +
                             dist_matrix_[a, u] < dist_matrix_[s, u]) and
                             s != t):
                                if current_neighbor.state == 0:  # NOT_IN_HEAP
                                    current_neighbor.state = 1   # IN_HEAP
                                    insert_node(heap, current_neighbor)
    return


def _update_insert(Py_ssize_t n, Py_ssize_t N,
                   np.ndarray[ITYPE_t, ndim=1, mode='c'] neighbors not None,
                   np.ndarray[DTYPE_t, ndim=2, mode='c'] dist_matrix_ not None,
                   np.ndarray[ITYPE_t, ndim=2, mode='c'] pred_ not None):
    '''
    Implementation of UpdateInsert from Law-Paper
    similar naming conventions.
    pi = self.predecessor_matrix_
    '''
    # Search over all vertex pairs that are neighbors of new point
    # insert these pairs into _update_insert_ab(), where path via new point is
    # better than old path
    cdef ITYPE_t a, b
    cdef Py_ssize_t i, j
    cdef Py_ssize_t neighbors_len = neighbors.shape[0]

    cdef FibonacciHeap heap
    cdef FibonacciNode* nodes = <FibonacciNode*> malloc((N) *
                                                    sizeof(FibonacciNode))

    for i in xrange(neighbors_len):
        for j in xrange(i, neighbors_len):
            a = neighbors[i]
            b = neighbors[j]
            if (dist_matrix_[a, n] + dist_matrix_[n, b] <
                dist_matrix_[a, b]):
                _update_insert_ab(a, b, n, N, dist_matrix_, pred_,
                                  &heap, nodes)
    free(nodes)
    return


def _update_edge(ITYPE_t a,
                 np.ndarray[ITYPE_t, ndim=1, mode='c'] neighbors_a not None,
                 Py_ssize_t n,
                 np.ndarray[ITYPE_t, ndim=1, mode='c'] neighbors_n not None,
                 Py_ssize_t N,
                 np.ndarray[DTYPE_t, ndim=2, mode='c'] dist_matrix_ not None,
                 np.ndarray[ITYPE_t, ndim=2, mode='c'] pred_ not None):
    '''
    Updating the geodesic distances if an edge (a,n) has been added.

    Search over all neighbors of the edge endpoint (n)
    insert these pairs into _update_insert_ab(), where path via new endpoint
    is better than old path. The same is done when switching roles of start (n)
    and endpoint (a)
    '''
    cdef ITYPE_t b
    cdef Py_ssize_t i
    cdef FibonacciHeap heap
    cdef FibonacciNode* nodes = <FibonacciNode*> malloc((N) *
                                                    sizeof(FibonacciNode))

    cdef Py_ssize_t neighbors_len = neighbors_n.shape[0]
    for i in xrange(neighbors_len):
        b = neighbors_n[i]
        if (dist_matrix_[a, n] + dist_matrix_[n, b] <
            dist_matrix_[a, b]):
            _update_insert_ab(a, b, n, N, dist_matrix_, pred_,
                              &heap, nodes)

    neighbors_len = neighbors_a.shape[0]
    for i in xrange(neighbors_len):
        b = neighbors_a[i]
        if (dist_matrix_[n, a] + dist_matrix_[a, b] <
            dist_matrix_[n, b]):
            _update_insert_ab(n, b, a, N, dist_matrix_, pred_,
                              &heap, nodes)
    free(nodes)
    return


def _determine_edge_changes(
                    ITYPE_t i_node,
                    Py_ssize_t N,
                    Py_ssize_t n_neighbors,
                    DTYPE_t[:, :] knn_dist,
                    ITYPE_t[:, :] knn_point,
                    kng not None,
                    DTYPE_t[:, :] training_data_):
    '''
    '''
    cdef ITYPE_t idx, max_idx, max_idx_full, max_j
    cdef DTYPE_t dist, x_dist, max_dist
    cdef Py_ssize_t i, j
    cdef ITYPE_t k

    if not isspmatrix_csr(kng):
        dist_matrix = csr_matrix(kng)

    cdef DTYPE_t[:] distances
    cdef ITYPE_t[:] neighbors, indptr
    distances = dist_matrix.data
    neighbors = dist_matrix.indices
    indptr = dist_matrix.indptr

    added_edges = []
    deleted_edges = []
    # Find added and deleted edges caused by the new point
    for i in xrange(n_neighbors):
        dist = knn_dist[0, i]
        idx = knn_point[0, i]
        added_edges.append(((i_node, idx), dist))
    # Compare if new point is in kNN of old points
    for i in xrange(N - 1):
        idx = knn_point[0, i]
        x_dist = knn_dist[0, i]
        # max for each row of the sparse matrix
        max_idx = -1
        max_dist = -1
        for k from indptr[idx] <= k < indptr[idx + 1]:
            if distances[k] > max_dist:
                max_dist = distances[k]
                max_idx = neighbors[k]
        if (indptr[idx + 1] - indptr[idx]) < n_neighbors:
            # useful for forgetting
            # add edge if idx does not have enough neighbors
            added_edges.append(((idx, i_node), x_dist))
        elif max_dist > x_dist:
            # new point is in kNN --> add edge from old point to new point
            added_edges.append(((idx, i_node), x_dist))
            # --> remove edge from old point to old kNN
            deleted_edges.append((idx, max_idx))
    return added_edges, deleted_edges
