#cython: boundscheck=False, nonecheck=False
# Minimum spanning tree single linkage implementation for hdbscan
# Authors: Leland McInnes, Steve Astels
# License: 3-clause BSD

cimport cython

import numpy as np
cimport numpy as np

cdef np.ndarray[np.double_t, ndim=2] mst_linkage_core(
                               np.ndarray[np.double_t, ndim=2] distance_matrix):

    cdef np.ndarray[np.int64_t, ndim=1] node_labels
    cdef np.ndarray[np.int64_t, ndim=1] current_labels
    cdef np.ndarray[np.double_t, ndim=1] current_distances
    cdef np.ndarray[np.double_t, ndim=1] left
    cdef np.ndarray[np.double_t, ndim=1] right
    cdef np.ndarray[np.double_t, ndim=2] result
    
    cdef np.ndarray label_filter
    
    cdef int current_node
    cdef int new_node_index
    cdef int new_node
    cdef int i
    
    result = np.zeros((distance_matrix.shape[0] - 1, 3))
    node_labels = np.arange(distance_matrix.shape[0], dtype=np.int64)
    current_node = 0
    current_distances = np.infty * np.ones(distance_matrix.shape[0])
    current_labels = node_labels
    for i in range(1,node_labels.shape[0]):
        label_filter = current_labels != current_node
        current_labels = current_labels[label_filter]
        left = current_distances[label_filter]
        right = distance_matrix[current_node][current_labels]
        current_distances = np.where(left < right, left, right)
        
        new_node_index = np.argmin(current_distances)
        new_node = current_labels[new_node_index]
        result[i - 1, 0] = <double> current_node
        result[i - 1, 1] = <double> new_node
        result[i - 1, 2] = current_distances[new_node_index]
        current_node = new_node
        
    return result
    
cdef class UnionFind (object):

    cdef np.ndarray parent
    cdef np.ndarray size
    cdef int next_label
    
    def __init__(self, N):
        self.parent = -1 * np.ones(2 * N - 1, dtype=np.int64)
        self.next_label = N
        self.size = np.hstack((np.ones(N, dtype=np.int64),
                               np.zeros(N-1, dtype=np.int64)))
                               
    cdef void union(self, int m, int n):
        self.size[self.next_label] = self.size[m] + self.size[n]
        self.parent[m] = self.next_label
        self.parent[n] = self.next_label
        self.size[self.next_label] = self.size[m] + self.size[n]
        self.next_label += 1
        
        return
        
    cdef int find(self, int n):
        while self.parent[n] != -1:
            n = self.parent[n]
        return n
        
    cdef int fast_find(self, int n):
        cdef int p
        p = n
        while self.parent[n] != -1:
            n = self.parent[n]
        # label up to the root
        while self.parent[p] != n:
            p, self.parent[p] = self.parent[p], n
        return n
        
cdef np.ndarray[np.double_t, ndim=2] label(np.ndarray[np.double_t, ndim=2] L, 
                                           do_fast_find=True):

    cdef np.ndarray[np.double_t, ndim=2] result

    cdef int N, a, aa, b, bb, idx
    cdef float delta
    
    result = np.zeros((L.shape[0], L.shape[1] + 1))
    N = L.shape[0] + 1
    U = UnionFind(N)
    
    for index, (a, b, delta) in enumerate(L):
        if do_fast_find:
            aa, bb = U.fast_find(a), U.fast_find(b)
        else:
            aa, bb, = U.find(a), U.find(b)
            
        result[index, 0] = aa
        result[index, 1] = bb
        result[index, 2] = delta
        result[index, 3] = U.size[aa] + U.size[bb]
        
        U.union(aa, bb)
       
    return result

cpdef np.ndarray[np.double_t, ndim=2] single_linkage(distance_matrix):
    
    cdef np.ndarray[np.double_t, ndim=2] hierarchy
    cdef np.ndarray[np.double_t, ndim=2] for_labelling
    
    hierarchy = mst_linkage_core(distance_matrix)
    for_labelling = hierarchy[np.argsort(hierarchy.T[2]), :]
    return label(for_labelling)
    
    
    
    
    
        
    
