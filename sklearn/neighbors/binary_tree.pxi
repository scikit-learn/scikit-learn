#!python

# KD Tree and Ball Tree
# =====================
#
#    Author: Jake Vanderplas <jakevdp@cs.washington.edu>, 2012-2013
#    License: BSD
#
# This file is meant to be a literal include in a pyx file.
# See ball_tree.pyx and kd_tree.pyx
#
# The routines here are the core algorithms of the KDTree and BallTree
# structures.  If Cython supported polymorphism, we would be able to
# create a subclass and derive KDTree and BallTree from it.  Because
# polymorphism is not an option, we use this single BinaryTree class
# as a literal include to avoid duplicating the entire file.
#
# A series of functions are implemented in kd_tree.pyx and ball_tree.pyx
# which use the information here to calculate the lower and upper bounds
# between a node and a point, and between two nodes.  These functions are
# used here, and are all that are needed to differentiate between the two
# tree types.
#
# Description of Binary Tree Algorithms
# -------------------------------------
# A binary tree can be thought of as a collection of nodes.  The top node
# contains all the points.  The next level consists of two nodes with half
# the points in each, and this continues recursively.  Each node contains
# metadata which allow fast computation of distance bounds: in the case of
# a ball tree, the metadata is a center and a radius.  In the case of a
# KD tree, the metadata is the minimum and maximum bound along each dimension.
#
# In a typical KD Tree or Ball Tree implementation, the nodes are implemented
# as dynamically allocated structures with pointers linking them.  Here we
# take a different approach, storing all relevant data in a set of arrays
# so that the entire tree object can be saved in a pickle file. For efficiency,
# the data can be stored in such a way that explicit pointers are not
# necessary: for node data stored at index i, the two child nodes are at
# index (2 * i + 1) and (2 * i + 2); the parent node is (i - 1) // 2
# (where // indicates integer division).
#
# The data arrays used here are as follows:
#   data : the [n_samples x n_features] array of data from which the tree
#          is built
#   idx_array : the length n_samples array used to keep track of the indices
#          of data within each node.  Each node has values idx_start and
#          idx_end: the points within the node are given by (using numpy
#          syntax) data[idx_array[idx_start:idx_end]].
#   node_data : the length n_nodes array of structures which store the node
#          indices, node radii, and leaf information for each node.
#   node_bounds : the [* x n_nodes x n_features] array containing the node
#          bound information.  For ball tree, the first dimension is 1, and
#          each row contains the centroid of the node.  For kd tree, the first
#          dimension is 2 and the rows for each point contain the arrays of
#          lower bounds and upper bounds in each direction.
#
# The lack of dynamic allocation means the number of nodes must be computed
# before the building of the tree. This can be done assuming the points are
# divided equally between child nodes at each step; although this removes
# some flexibility in tree creation, it ensures a balanced tree and ensures
# that the number of nodes required can be computed beforehand.  Given a
# specified leaf_size (the minimum number of points in any node), it is
# possible to show that a balanced tree will have
#
#     n_levels = 1 + max(0, floor(log2((n_samples - 1) / leaf_size)))
#
# in order to satisfy
#
#     leaf_size <= min(n_points) <= 2 * leaf_size
#
# with the exception of the special case where n_samples < leaf_size.
# for a given number of levels, the number of nodes in the tree is given by
#
#     n_nodes = 2 ** n_levels - 1
#
# both these results can be straightforwardly shown by induction.  The
# following code uses these values in the construction of the tree.
#
# Distance Metrics
# ----------------
# For flexibility, the trees can be built using a variety of distance metrics.
# The metrics are described in the DistanceMetric class: the standard
# Euclidean distance is the default, and is inlined to be faster than other
# metrics.  In addition, each metric defines both a distance and a
# "reduced distance", which is often faster to compute, and is therefore
# used in the query architecture whenever possible. (For example, in the
# case of the standard Euclidean distance, the reduced distance is the
# squared-distance).

cimport cython
cimport numpy as np
from libc.math cimport fmax, fmin, fabs, sqrt, exp, cos, pow

import numpy as np
import warnings
from ..utils import array2d

from typedefs cimport DTYPE_t, ITYPE_t, DITYPE_t
from typedefs import DTYPE, ITYPE

from dist_metrics cimport (DistanceMetric, euclidean_dist, euclidean_rdist,
                           euclidean_dist_to_rdist, euclidean_rdist_to_dist)

# some handy constants
cdef DTYPE_t INF = np.inf
cdef DTYPE_t PI = np.pi
cdef DTYPE_t ROOT_2PI = sqrt(2 * PI)

######################################################################
# Define doc strings, substituting the appropriate class name using
# the DOC_DICT variable defined in the pyx files.
CLASS_DOC = \
"""{BinaryTree} for fast generalized N-point problems

{BinaryTree}(X, leaf_size=40, metric='minkowski', **kwargs)

Parameters
----------
X : array-like, shape = [n_samples, n_features]
    n_samples is the number of points in the data set, and
    n_features is the dimension of the parameter space.
    Note: if X is a C-contiguous array of doubles then data will
    not be copied. Otherwise, an internal copy will be made.

leaf_size : positive integer (default = 20)
    Number of points at which to switch to brute-force. Changing
    leaf_size will not affect the results of a query, but can
    significantly impact the speed of a query and the memory required
    to store the constructed tree.  The amount of memory needed to
    store the tree scales as approximately n_samples / leaf_size.
    For a specified ``leaf_size``, a leaf node is guaranteed to
    satisfy ``leaf_size <= n_points <= 2 * leaf_size``, except in
    the case that ``n_samples < leaf_size``.

metric : string or DistanceMetric object
    the distance metric to use for the tree.  Default='minkowski'
    with p=2 (that is, a euclidean metric). See the documentation
    of the DistanceMetric class for a list of available metrics.
    {binary_tree}.valid_metrics gives a list of the metrics which
    are valid for {BinaryTree}.

Additional keywords are passed to the distance metric class.

Attributes
----------
data : np.ndarray
    The training data

Examples
--------
Query for k-nearest neighbors

    >>> import numpy as np

    >>> np.random.seed(0)
    >>> X = np.random.random((10,3))  # 10 points in 3 dimensions
    >>> tree = {BinaryTree}(X, leaf_size=2)              # doctest: +SKIP
    >>> dist, ind = tree.query(X[0], k=3)                # doctest: +SKIP
    >>> print ind  # indices of 3 closest neighbors
    [0 3 1]
    >>> print dist  # distances to 3 closest neighbors
    [ 0.          0.19662693  0.29473397]

Pickle and Unpickle a tree.  Note that the state of the tree is saved in the
pickle operation: the tree needs not be rebuilt upon unpickling.

    >>> import numpy as np
    >>> import pickle
    >>> np.random.seed(0)
    >>> X = np.random.random((10,3))  # 10 points in 3 dimensions
    >>> tree = {BinaryTree}(X, leaf_size=2)        # doctest: +SKIP
    >>> s = pickle.dumps(tree)                     # doctest: +SKIP
    >>> tree_copy = pickle.loads(s)                # doctest: +SKIP
    >>> dist, ind = tree_copy.query(X[0], k=3)     # doctest: +SKIP
    >>> print ind  # indices of 3 closest neighbors   
    [0 3 1]
    >>> print dist  # distances to 3 closest neighbors
    [ 0.          0.19662693  0.29473397]

Query for neighbors within a given radius

    >>> import numpy as np
    >>> np.random.seed(0)
    >>> X = np.random.random((10,3))  # 10 points in 3 dimensions
    >>> tree = BinaryTree(X, leaf_size=2)     # doctest: +SKIP
    >>> print tree.query_radius(X[0], r=0.3, count_only=True)
    3
    >>> ind = tree.query_radius(X[0], r=0.3)  # doctest: +SKIP
    >>> print ind  # indices of neighbors within distance 0.3
    [3 0 1]


Compute a gaussian kernel density estimate:

    >>> import numpy as np
    >>> np.random.seed(1)
    >>> X = np.random.random((100, 3))
    >>> tree = BinaryTree(X)                  # doctest: +SKIP
    >>> tree.kernel_density(X[:3], h=0.1, kernel='gaussian')
    array([ 6.94114649,  7.83281226,  7.2071716 ])

Compute a two-point auto-correlation function

    >>> import numpy as np
    >>> np.random.seed(0)
    >>> X = np.random.random((30, 3))
    >>> r = np.linspace(0, 1, 5)
    >>> tree = BinaryTree(X)                  # doctest: +SKIP
    >>> tree.two_point_correlation(X, r)
    array([ 30,  62, 278, 580, 820])

""".format(**DOC_DICT)


######################################################################
# Kernel functions
#
# Note: Kernels assume dist is non-negative and and h is positive
#       All kernel functions are normalized such that K(0, h) = 1.
#       The fully normalized kernel is:
#         kernel_norm(h, kernel) * compute_kernel(dist, h, kernel)
#       The code only works with non-negative kernels: i.e. K(d, h) >= 0
#       for all valid d and h.
cdef enum KernelType:
    GAUSSIAN_KERNEL = 1
    TOPHAT_KERNEL = 2
    EPANECHNIKOV_KERNEL = 3
    EXPONENTIAL_KERNEL = 4
    LINEAR_KERNEL = 5
    COSINE_KERNEL = 6

cdef inline DTYPE_t gaussian_kernel(DTYPE_t dist, DTYPE_t h):
    return exp(-0.5 * (dist * dist) / (h * h))

cdef inline DTYPE_t tophat_kernel(DTYPE_t dist, DTYPE_t h):
    if dist < h:
        return 1.0
    else:
        return 0.0

cdef inline DTYPE_t epanechnikov_kernel(DTYPE_t dist, DTYPE_t h):
    if dist < h:
        return 1.0 - (dist * dist) / (h * h)
    else:
        return 0.0

cdef inline DTYPE_t exponential_kernel(DTYPE_t dist, DTYPE_t h):
    if dist < h:
        return exp(-dist / h)
    else:
        return 0.0

cdef inline DTYPE_t linear_kernel(DTYPE_t dist, DTYPE_t h):
    if dist < h:
        return 1 - dist / h
    else:
        return 0.0

cdef inline DTYPE_t cosine_kernel(DTYPE_t dist, DTYPE_t h):
    if dist < h:
        return cos(0.5 * PI * dist / h)
    else:
        return 0.0


cdef inline DTYPE_t compute_kernel(DTYPE_t dist, DTYPE_t h,
                                   KernelType kernel):
    if kernel == GAUSSIAN_KERNEL:
        return gaussian_kernel(dist, h)
    elif kernel == TOPHAT_KERNEL:
        return tophat_kernel(dist, h)
    if kernel == EPANECHNIKOV_KERNEL:
        return epanechnikov_kernel(dist, h)
    elif kernel == EXPONENTIAL_KERNEL:
        return exponential_kernel(dist, h)
    if kernel == LINEAR_KERNEL:
        return linear_kernel(dist, h)
    elif kernel == COSINE_KERNEL:
        return cosine_kernel(dist, h)


cdef inline DTYPE_t kernel_norm(DTYPE_t h, KernelType kernel):
    if kernel == GAUSSIAN_KERNEL:
        return 1. / (h * ROOT_2PI)
    elif kernel == TOPHAT_KERNEL:
        return 0.5 / h
    if kernel == EPANECHNIKOV_KERNEL:
        return 0.75 / h
    elif kernel == EXPONENTIAL_KERNEL:
        return 0.5 / h
    if kernel == LINEAR_KERNEL:
        return 1.0 / h
    elif kernel == COSINE_KERNEL:
        return 0.25 * PI / h


######################################################################
# Tree Utility Routines
cdef inline void swap(DITYPE_t* arr, ITYPE_t i1, ITYPE_t i2):
    cdef DITYPE_t tmp = arr[i1]
    arr[i1] = arr[i2]
    arr[i2] = tmp


cdef inline void dual_swap(DTYPE_t* darr, ITYPE_t* iarr,
                           ITYPE_t i1, ITYPE_t i2):
    cdef DTYPE_t dtmp = darr[i1]
    darr[i1] = darr[i2]
    darr[i2] = dtmp

    cdef ITYPE_t itmp = iarr[i1]
    iarr[i1] = iarr[i2]
    iarr[i2] = itmp
    

#------------------------------------------------------------
# NeighborsHeap
#  max-heap structure to keep track of distances and indices of neighbors
cdef class NeighborsHeap:
    cdef DTYPE_t[:, ::1] distances
    cdef ITYPE_t[:, ::1] indices

    def __cinit__(self):
        self.distances = np.zeros((1, 1), dtype=DTYPE)
        self.indices = np.zeros((1, 1), dtype=ITYPE)

    def __init__(self, n_pts, n_nbrs):
        self.distances = np.zeros((n_pts, n_nbrs), dtype=DTYPE) + np.inf
        self.indices = np.zeros((n_pts, n_nbrs), dtype=ITYPE)

    def get_arrays(self, sort=True):
        if sort:
            self._sort()
        return map(np.asarray, (self.distances, self.indices))

    cdef inline DTYPE_t largest(self, ITYPE_t row):
        return self.distances[row, 0]

    cpdef push(self, ITYPE_t row, DTYPE_t val, ITYPE_t i_val):
        cdef ITYPE_t i, ic1, ic2, i_swap
        cdef ITYPE_t size = self.distances.shape[1]
        cdef DTYPE_t* dist_arr = &self.distances[row, 0]
        cdef ITYPE_t* ind_arr = &self.indices[row, 0]

        # check if val should be in heap
        if val > dist_arr[0]:
            return

        # insert val at position zero
        dist_arr[0] = val
        ind_arr[0] = i_val

        #descend the heap, swapping values until the max heap criterion is met
        i = 0
        while True:
            ic1 = 2 * i + 1
            ic2 = ic1 + 1

            if ic1 >= size:
                break
            elif ic2 >= size:
                if dist_arr[ic1] > val:
                    i_swap = ic1
                else:
                    break
            elif dist_arr[ic1] >= dist_arr[ic2]:
                if val < dist_arr[ic1]:
                    i_swap = ic1
                else:
                    break
            else:
                if val < dist_arr[ic2]:
                    i_swap = ic2
                else:
                    break

            dist_arr[i] = dist_arr[i_swap]
            ind_arr[i] = ind_arr[i_swap]

            i = i_swap

        dist_arr[i] = val
        ind_arr[i] = i_val

    cdef _sort(self):
        cdef DTYPE_t[:, ::1] distances = self.distances
        cdef ITYPE_t[:, ::1] indices = self.indices
        cdef ITYPE_t row
        for row in range(distances.shape[0]):
            _simultaneous_sort(&distances[row, 0],
                               &indices[row, 0],
                               distances.shape[1])


#------------------------------------------------------------
# simultaneous_sort :
#  this is a simple recursive quicksort implementation which sorts
#  `distances` and simultaneously performs the same swaps on `indices`.
cdef void _simultaneous_sort(DTYPE_t* dist, ITYPE_t* idx, ITYPE_t size):
    cdef ITYPE_t pivot_idx, i, store_idx
    cdef DTYPE_t pivot_val

    # in the small-array case, do things efficiently
    if size <= 1:
        pass
    elif size == 2:
        if dist[0] > dist[1]:
            dual_swap(dist, idx, 0, 1)
    elif size == 3:
        if dist[0] > dist[1]:
            dual_swap(dist, idx, 0, 1)
        if dist[1] > dist[2]:
            dual_swap(dist, idx, 1, 2)
            if dist[0] > dist[1]:
                dual_swap(dist, idx, 0, 1)
    else:
        # Determine the pivot using the median-of-three rule.
        # The smallest of the three is moved to the beginning of the array,
        # the middle (the pivot value) is moved to the end, and the largest
        # is moved to the pivot index.
        pivot_idx = size / 2
        if dist[0] > dist[size - 1]:
            dual_swap(dist, idx, 0, size - 1)
        if dist[size - 1] > dist[pivot_idx]:
            dual_swap(dist, idx, size - 1, pivot_idx)
            if dist[0] > dist[size - 1]:
                dual_swap(dist, idx, 0, size - 1)
        pivot_val = dist[size - 1]

        # partition indices about pivot.  At the end of this operation,
        # pivot_idx will contain the pivot value, everything to the left
        # will be smaller, and everything to the right will be larger.
        store_idx = 0
        for i in range(size - 1):
            if dist[i] < pivot_val:
                dual_swap(dist, idx, i, store_idx)
                store_idx += 1
        dual_swap(dist, idx, store_idx, size - 1)
        pivot_idx = store_idx

        # recursively sort each side of the pivot
        if pivot_idx > 1:
            _simultaneous_sort(dist, idx, pivot_idx)
        if pivot_idx + 2 < size:
            _simultaneous_sort(dist + pivot_idx + 1,
                               idx + pivot_idx + 1,
                               size - pivot_idx - 1)


#------------------------------------------------------------
# find_node_split_dim:
#  this computes the equivalent of the 
#  j_max = np.argmax(np.max(data, 0) - np.min(data, 0))
cdef ITYPE_t find_node_split_dim(DTYPE_t* data,
                                 ITYPE_t* node_indices,
                                 ITYPE_t n_features,
                                 ITYPE_t n_points):
    cdef DTYPE_t min_val, max_val, val, spread, max_spread
    cdef ITYPE_t i, j, j_max

    j_max = 0
    max_spread = 0

    for j in range(n_features):
        max_val = data[node_indices[0] * n_features + j]
        min_val = max_val
        for i in range(1, n_points):
            val = data[node_indices[i] * n_features + j]
            max_val = fmax(max_val, val)
            min_val = fmin(min_val, val)
        spread = max_val - min_val
        if spread > max_spread:
            max_spread = spread
            j_max = j
    return j_max


#------------------------------------------------------------
# partition_node_indices will modify the array node_indices between
# indices 0 and n_points.  Upon return (assuming numpy-style slicing)
#   data[node_indices[0:split_index], split_dim]
#     <= data[node_indices[split_index], split_dim]
# and
#   data[node_indices[split_index], split_dim]
#     <= data[node_indices[split_index:n_points], split_dim]
# will hold.  The algorithm amounts to a partial quicksort
cdef void partition_node_indices(DTYPE_t* data,
                                 ITYPE_t* node_indices,
                                 ITYPE_t split_dim,
                                 ITYPE_t split_index,
                                 ITYPE_t n_features,
                                 ITYPE_t n_points):
    cdef ITYPE_t left, right, midindex, i
    cdef DTYPE_t d1, d2
    left = 0
    right = n_points - 1

    while True:
        midindex = left
        for i in range(left, right):
            d1 = data[node_indices[i] * n_features + split_dim]
            d2 = data[node_indices[right] * n_features + split_dim]
            if d1 < d2:
                swap(node_indices, i, midindex)
                midindex += 1
        swap(node_indices, midindex, right)
        if midindex == split_index:
            break
        elif midindex < split_index:
            left = midindex + 1
        else:
            right = midindex - 1

######################################################################
# NodeHeap : min-heap used to keep track of nodes during
#            breadth-first query

cdef struct NodeHeapData_t:
    DTYPE_t val
    ITYPE_t i1
    ITYPE_t i2

# use a dummy variable to determine the python data type
cdef NodeHeapData_t ndummy
cdef NodeHeapData_t[:] ndummy_view = <NodeHeapData_t[:1]> &ndummy
NodeHeapData = np.asarray(ndummy_view).dtype


cdef inline void swap_nodes(NodeHeapData_t* arr, ITYPE_t i1, ITYPE_t i2):
    cdef NodeHeapData_t tmp = arr[i1]
    arr[i1] = arr[i2]
    arr[i2] = tmp


cdef class NodeHeap:
    """NodeHeap

    This is a min-heap implementation for keeping track of nodes
    during a breadth-first search.

    The heap is an array which is maintained so that the min heap condition
    is met:  heap[i].val < min(heap[2 * i + 1].val, heap[2 * i + 2].val)
    """
    cdef NodeHeapData_t[::1] data
    cdef ITYPE_t n
    def __cinit__(self):
        self.data = np.zeros(0, dtype=NodeHeapData)

    def __init__(self, size_guess=100):
        self.data = np.zeros(size_guess, dtype=NodeHeapData)
        self.n = size_guess
        self.clear()

    cdef void resize(self, ITYPE_t new_size):
        cdef NodeHeapData_t *data_arr, *new_data_arr
        cdef ITYPE_t i
        cdef ITYPE_t size = self.data.shape[0]
        cdef NodeHeapData_t[::1] new_data = np.zeros(new_size,
                                                     dtype=NodeHeapData)

        if size > 0 and new_size > 0:
            data_arr = &self.data[0]
            new_data_arr = &new_data[0]
            for i in range(min(size, new_size)):
                new_data_arr[i] = data_arr[i]
            
        if new_size < size:
            self.n = new_size

        self.data = new_data

    cdef void push(self, NodeHeapData_t data):
        cdef ITYPE_t i, i_parent
        cdef NodeHeapData_t* data_arr
        self.n += 1
        if self.n > self.data.shape[0]:
            self.resize(2 * self.n)

        # put the new element at the end,
        # and then perform swaps until the heap is in order
        data_arr = &self.data[0]
        i = self.n - 1
        data_arr[i] = data
        
        while i > 0:
            i_parent = (i - 1) // 2
            if data_arr[i_parent].val <= data_arr[i].val:
                break
            else:
                swap_nodes(data_arr, i, i_parent)
                i = i_parent

    cdef NodeHeapData_t peek(self):
        return self.data[0]

    cdef NodeHeapData_t pop(self):
        if self.n == 0:
            raise ValueError('cannot pop on empty heap')

        cdef ITYPE_t i, i_child1, i_child2, i_swap
        cdef NodeHeapData_t* data_arr = &self.data[0]
        cdef NodeHeapData_t popped_element = data_arr[0]

        # pop off the first element, move the last element to the front,
        # and then perform swaps until the heap is back in order
        data_arr[0] = data_arr[self.n - 1]
        self.n -= 1

        i = 0

        while (i < self.n):
            i_child1 = 2 * i + 1
            i_child2 = 2 * i + 2
            i_swap = 0

            if i_child2 < self.n:
                if data_arr[i_child1].val <= data_arr[i_child2].val:
                    i_swap = i_child1
                else:
                    i_swap = i_child2
            elif i_child1 < self.n:
                i_swap = i_child1
            else:
                break

            if (i_swap > 0) and (data_arr[i_swap].val <= data_arr[i].val):
                swap_nodes(data_arr, i, i_swap)
                i = i_swap
            else:
                break
        
        return popped_element

    cdef clear(self):
        """Clear the stack"""
        self.n = 0


######################################################################
# newObj function
#  this is a helper function for pickling
def newObj(obj):
    return obj.__new__(obj)


######################################################################
# Binary Tree class

cdef struct NodeData_t:
    ITYPE_t idx_start
    ITYPE_t idx_end
    int is_leaf
    DTYPE_t radius

# use a dummy variable to determine the python data type
cdef NodeData_t dummy
cdef NodeData_t[:] dummy_view = <NodeData_t[:1]> &dummy
NodeData = np.asarray(dummy_view).dtype


cdef class BinaryTree:
    __doc__ = CLASS_DOC

    cdef public object valid_metrics

    cdef readonly DTYPE_t[:, ::1] data
    cdef public ITYPE_t[::1] idx_array
    cdef public NodeData_t[::1] node_data
    cdef public DTYPE_t[:, :, ::1] node_bounds

    cdef ITYPE_t leaf_size
    cdef ITYPE_t n_levels
    cdef ITYPE_t n_nodes

    cdef DistanceMetric dm
    cdef int euclidean

    # variables to keep track of building & querying stats
    cdef int n_trims
    cdef int n_leaves
    cdef int n_splits
    cdef int n_calls

    # Use cinit to initialize all arrays to empty: this will prevent memory
    # errors and seg-faults in rare cases where __init__ is not called
    def __cinit__(self):
        self.data = np.empty((0, 1), dtype=DTYPE, order='C')
        self.idx_array = np.empty(0, dtype=ITYPE, order='C')
        self.node_data = np.empty(0, dtype=NodeData, order='C')
        self.node_bounds = np.empty((0, 1, 1), dtype=DTYPE)

        self.leaf_size = 0
        self.n_levels = 0
        self.n_nodes = 0

        self.euclidean = False

        self.n_trims = 0
        self.n_leaves = 0
        self.n_splits = 0
        self.n_calls = 0

    def __init__(self, data,
                 leaf_size=40, metric='minkowski', **kwargs):
        self.valid_metrics = VALID_METRICS
        self.data = np.asarray(data, dtype=DTYPE, order='C')
        self.leaf_size = leaf_size
        self.dm = DistanceMetric.get_metric(metric, **kwargs)
        self.euclidean = (self.dm.__class__.__name__ == 'EuclideanDistance')

        metric = self.dm.__class__.__name__
        if metric not in self.valid_metrics:
            raise ValueError('metric {metric} is not valid for '
                             '{BinaryTree}'.format(metric=metric,
                                                   **DOC_DICT))

        # validate data
        if self.data.size == 0:
            raise ValueError("X is an empty array")

        if leaf_size < 1:
            raise ValueError("leaf_size must be greater than or equal to 1")
        
        n_samples = self.data.shape[0]
        n_features = self.data.shape[1]

        # determine number of levels in the tree, and from this
        # the number of nodes in the tree.  This results in leaf nodes
        # with numbers of points betweeen leaf_size and 2 * leaf_size
        self.n_levels = np.log2(fmax(1, (n_samples - 1) / self.leaf_size)) + 1
        self.n_nodes = (2 ** self.n_levels) - 1

        # allocate arrays for storage
        self.idx_array = np.arange(n_samples, dtype=ITYPE)
        self.node_data = np.zeros(self.n_nodes, dtype=NodeData)

        # Allocate tree-specific data
        allocate_data(self, self.n_nodes, n_features)        
        self._recursive_build(0, 0, n_samples)

    def __reduce__(self):
        """
        reduce method used for pickling
        """
        return (newObj, (BinaryTree,), self.__getstate__())

    def __getstate__(self):
        """
        get state for pickling
        """
        return (np.asarray(self.data),
                np.asarray(self.idx_array),
                np.asarray(self.node_data),
                np.asarray(self.node_bounds),
                int(self.leaf_size),
                int(self.n_levels),
                int(self.n_nodes),
                int(self.n_trims),
                int(self.n_leaves),
                int(self.n_splits),
                int(self.n_calls),
                self.dm)

    def __setstate__(self, state):
        """
        set state for pickling
        """
        self.data = state[0]
        self.idx_array = state[1]
        self.node_data = state[2]
        self.node_bounds = state[3]
        self.leaf_size = state[4]
        self.n_levels = state[5]
        self.n_nodes = state[6]
        self.n_trims = state[7]
        self.n_leaves = state[8]
        self.n_splits = state[9]
        self.n_calls = state[10]
        self.dm = state[11]
        self.euclidean = (self.dm.__class__.__name__ == 'EuclideanDistance')

    def get_tree_stats(self):
        return (self.n_trims, self.n_leaves, self.n_splits)

    def reset_n_calls(self):
        self.n_calls = 0

    def get_n_calls(self):
        return self.n_calls

    def get_arrays(self):
        return map(np.asarray, (self.data, self.idx_array,
                                self.node_data, self.node_bounds))        

    cdef inline DTYPE_t dist(self, DTYPE_t* x1, DTYPE_t* x2,
                             ITYPE_t size):
        self.n_calls += 1
        if self.euclidean:
            return euclidean_dist(x1, x2, size)
        else:
            return self.dm.dist(x1, x2, size)

    cdef inline DTYPE_t rdist(self, DTYPE_t* x1, DTYPE_t* x2,
                              ITYPE_t size):
        self.n_calls += 1
        if self.euclidean:
            return euclidean_rdist(x1, x2, size)
        else:
            return self.dm.rdist(x1, x2, size)

    cdef void _recursive_build(self, ITYPE_t i_node,
                               ITYPE_t idx_start, ITYPE_t idx_end):
        cdef ITYPE_t imax
        cdef ITYPE_t n_features = self.data.shape[1]
        cdef ITYPE_t n_points = idx_end - idx_start
        cdef ITYPE_t n_mid = n_points / 2
        cdef ITYPE_t* idx_array = &self.idx_array[idx_start]
        cdef DTYPE_t* data = &self.data[0, 0]

        # initialize node data
        init_node(self, i_node, idx_start, idx_end)

        if 2 * i_node + 1 >= self.n_nodes:
            self.node_data[i_node].is_leaf = True
            if idx_end - idx_start > 2 * self.leaf_size:
                # this shouldn't happen if our memory allocation is correct
                # we'll proactively prevent memory errors, but raise a
                # warning saying we're doing so.
                import warnings
                warnings.warn("Internal: memory layout is flawed: "
                              "not enough nodes allocated")

        elif idx_end - idx_start < 2:
            # again, this shouldn't happen if our memory allocation
            # is correct.  Raise a warning.
            import warnings
            warnings.warn("Internal: memory layout is flawed: "
                          "too many nodes allocated")
            self.node_data[i_node].is_leaf = True

        else:
            # split node and recursively construct child nodes.
            self.node_data[i_node].is_leaf = False
            i_max = find_node_split_dim(data, idx_array,
                                        n_features, n_points)
            partition_node_indices(data, idx_array, i_max, n_mid,
                                   n_features, n_points)
            self._recursive_build(2 * i_node + 1,
                                  idx_start, idx_start + n_mid)
            self._recursive_build(2 * i_node + 2,
                                  idx_start + n_mid, idx_end)

    def query(self, X, k=1, return_distance=True,
              dualtree=False, breadth_first=False,
              sort_results=True):
        """
        query(X, k=1, return_distance=True,
              dualtree=False, breadth_first=False)

        query the treeree for the k nearest neighbors

        Parameters
        ----------
        X : array-like, last dimension self.dim
            An array of points to query
        k : integer  (default = 1)
            The number of nearest neighbors to return
        return_distance : boolean (default = True)
            if True, return a tuple (d, i) of distances and indices
            if False, return array i
        dualtree : boolean (default = False)
            if True, use the dual tree formalism for the query: a tree is
            built for the query points, and the pair of trees is used to
            efficiently search this space.  This can lead to better
            performance as the number of points grows large.
        breadth_first : boolean (default = False)
            if True, then query the nodes in a breadth-first manner.
            Otherwise, query the nodes in a depth-first manner.
        sort_results : boolean (default = True)
            if True, then distances and indices of each point are sorted
            on return, so that the first column contains the closest points.
            Otherwise, neighbors are returned in an arbitrary order.

        Returns
        -------
        i    : if return_distance == False
        (d,i) : if return_distance == True

        d : array of doubles - shape: x.shape[:-1] + (k,)
            each entry gives the list of distances to the
            neighbors of the corresponding point

        i : array of integers - shape: x.shape[:-1] + (k,)
            each entry gives the list of indices of
            neighbors of the corresponding point

        Examples
        --------
        Query for k-nearest neighbors

            >>> import numpy as np
            >>> np.random.seed(0)
            >>> X = np.random.random((10,3))  # 10 points in 3 dimensions
            >>> tree = BinaryTree(X, leaf_size=2)    # doctest: +SKIP
            >>> dist, ind = tree.query(X[0], k=3)    # doctest: +SKIP
            >>> print ind  # indices of 3 closest neighbors
            [0 3 1]
            >>> print dist  # distances to 3 closest neighbors
            [ 0.          0.19662693  0.29473397]
        """
        # XXX: we should allow X to be a pre-built tree.
        X = array2d(X, dtype=DTYPE, order='C')

        if X.shape[-1] != self.data.shape[1]:
            raise ValueError("query data dimension must "
                             "match training data dimension")

        if self.data.shape[0] < k:
            raise ValueError("k must be less than or equal "
                             "to the number of training points")

        # flatten X, and save original shape information
        cdef DTYPE_t[:, ::1] Xarr = X.reshape((-1, self.data.shape[1]))
        cdef DTYPE_t reduced_dist_LB
        cdef ITYPE_t i
        cdef DTYPE_t* pt

        # initialize heap for neighbors
        cdef NeighborsHeap heap = NeighborsHeap(Xarr.shape[0], k)

        # node heap for breadth-first queries
        cdef NodeHeap nodeheap
        if breadth_first:
            nodeheap = NodeHeap(self.data.shape[0] // self.leaf_size)

        # bounds is needed for the dual tree algorithm
        cdef DTYPE_t[::1] bounds

        self.n_trims = 0
        self.n_leaves = 0
        self.n_splits = 0

        if dualtree:
            other = self.__class__(Xarr, metric=self.dm,
                                   leaf_size=self.leaf_size)
            if breadth_first:
                self._query_dual_breadthfirst(other, heap, nodeheap)
            else:
                reduced_dist_LB = min_rdist_dual(self, 0, other, 0)
                bounds = np.inf + np.zeros(other.node_data.shape[0])
                self._query_dual_depthfirst(0, other, 0, bounds,
                                            heap, reduced_dist_LB)

        else:
            pt = &Xarr[0, 0]
            if breadth_first:
                for i in range(Xarr.shape[0]):
                    self._query_single_breadthfirst(pt, i, heap, nodeheap)
                    pt += Xarr.shape[1]
            else:
                for i in range(Xarr.shape[0]):
                    reduced_dist_LB = min_rdist(self, 0, pt)
                    self._query_single_depthfirst(0, pt, i, heap,
                                                  reduced_dist_LB)
                    pt += Xarr.shape[1]

        distances, indices = heap.get_arrays(sort=sort_results)
        distances = self.dm.rdist_to_dist_arr(distances)

        # deflatten results
        if return_distance:
            return (distances.reshape(X.shape[:-1] + (k,)),
                    indices.reshape(X.shape[:-1] + (k,)))
        else:
            return indices.reshape(X.shape[:-1] + (k,))

    def query_radius(self, X, r, return_distance=False,
                     int count_only=False, int sort_results=False):
        """
        query_radius(self, X, r, count_only = False):

        query the tree for neighbors within a radius r

        Parameters
        ----------
        X : array-like, last dimension self.dim
            An array of points to query
        r : distance within which neighbors are returned
            r can be a single value, or an array of values of shape
            x.shape[:-1] if different radii are desired for each point.
        return_distance : boolean (default = False)
            if True,  return distances to neighbors of each point
            if False, return only neighbors
            Note that unlike the query() method, setting return_distance=True
            here adds to the computation time.  Not all distances need to be
            calculated explicitly for return_distance=False.  Results are
            not sorted by default: see ``sort_results`` keyword.
        count_only : boolean (default = False)
            if True,  return only the count of points within distance r
            if False, return the indices of all points within distance r
            If return_distance==True, setting count_only=True will
            result in an error.
        sort_results : boolean (default = False)
            if True, the distances and indices will be sorted before being
            returned.  If False, the results will not be sorted.  If
            return_distance == False, setting sort_results = True will
            result in an error.

        Returns
        -------
        count       : if count_only == True
        ind         : if count_only == False and return_distance == False
        (ind, dist) : if count_only == False and return_distance == True

        count : array of integers, shape = X.shape[:-1]
            each entry gives the number of neighbors within
            a distance r of the corresponding point.

        ind : array of objects, shape = X.shape[:-1]
            each element is a numpy integer array listing the indices of
            neighbors of the corresponding point.  Note that unlike
            the results of a k-neighbors query, the returned neighbors
            are not sorted by distance by default.

        dist : array of objects, shape = X.shape[:-1]
            each element is a numpy double array
            listing the distances corresponding to indices in i.

        Examples
        --------
        Query for neighbors in a given radius

        >>> import numpy as np
        >>> np.random.seed(0)
        >>> X = np.random.random((10,3))  # 10 points in 3 dimensions
        >>> tree = BinaryTree(X, leaf_size=2)     # doctest: +SKIP
        >>> print tree.query_radius(X[0], r=0.3, count_only=True)
        3
        >>> ind = tree.query_radius(X[0], r=0.3)  # doctest: +SKIP
        >>> print ind  # indices of neighbors within distance 0.3
        [3 0 1]
        """
        if count_only and return_distance:
            raise ValueError("count_only and return_distance "
                             "cannot both be true")

        if sort_results and not return_distance:
            raise ValueError("return_distance must be True "
                             "if sort_results is True")

        cdef ITYPE_t i, count_i = 0
        cdef ITYPE_t n_features = self.data.shape[1]
        cdef DTYPE_t[::1] dist_arr_i
        cdef ITYPE_t[::1] idx_arr_i, count_arr
        cdef DTYPE_t* pt

        # validate X and prepare for query
        X = array2d(X, dtype=DTYPE, order='C')

        if X.shape[-1] != self.data.shape[1]:
            raise ValueError("query data dimension must "
                             "match training data dimension")

        cdef DTYPE_t[:, ::1] Xarr = X.reshape((-1, self.data.shape[1]))

        # prepare r for query
        r = np.asarray(r, dtype=DTYPE, order='C')
        r = np.atleast_1d(r)
        if r.shape == (1,):
            r = r[0] + np.zeros(X.shape[:-1], dtype=DTYPE)
        else:
            if r.shape != X.shape[:-1]:
                raise ValueError("r must be broadcastable to X.shape")

        cdef DTYPE_t[::1] rarr = r.reshape(-1)

        # prepare variables for iteration
        if not count_only:
            indices = np.zeros(Xarr.shape[0], dtype='object')
            if return_distance:
                distances = np.zeros(Xarr.shape[0], dtype='object')

        np_idx_arr = np.zeros(self.data.shape[0], dtype=ITYPE)
        idx_arr_i = np_idx_arr

        np_dist_arr = np.zeros(self.data.shape[0], dtype=DTYPE)
        dist_arr_i = np_dist_arr

        count_arr = np.zeros(Xarr.shape[0], dtype=ITYPE)

        pt = &Xarr[0, 0]
        for i in range(Xarr.shape[0]):
            count_arr[i] = self._query_radius_single(0, pt, rarr[i],
                                                     &idx_arr_i[0],
                                                     &dist_arr_i[0],
                                                     0, count_only,
                                                     return_distance)
            pt += n_features

            if count_only:
                pass
            else:
                if sort_results:
                    _simultaneous_sort(&dist_arr_i[0], &idx_arr_i[0],
                                       count_arr[i])
 
                indices[i] = np_idx_arr[:count_arr[i]].copy()
                if return_distance:
                    distances[i] = np_dist_arr[:count_arr[i]].copy()

        # deflatten results
        if count_only:
            return np.asarray(count_arr).reshape(X.shape[:-1])
        elif return_distance:
            return (indices.reshape(X.shape[:-1]),
                    distances.reshape(X.shape[:-1]))
        else:
            return indices.reshape(X.shape[:-1])

    def kernel_density(BinaryTree self, X, h, kernel='gaussian',
                       atol=0, rtol=0, dualtree=False, breadth_first=False):
        """
        kernel_density(self, X, h, kernel='gaussian', atol=0, rtol=0,
                       dualtree=False)

        Compute the kernel density estimate at points X with the given kernel,
        using the distance metric specified at tree creation.

        Parameters
        ----------
        X : array_like
            An array of points to query.  Last dimension should match dimension
            of training data.
        h : float
            the bandwidth of the kernel
        kernel : string
            specify the kernel to use.  Options are
            - 'gaussian'
            - 'tophat'
            - 'epanechnikov'
            - 'exponential'
            - 'linear'
            - 'cosine'
            Default is kernel = 'gaussian'
        atol, rtol : float (default = 0)
            Specify the desired relative and absolute tolerance of the result.
            If the true result is K_true, then the returned result K_ret
            satisfies ``abs(K_true - K_ret) < atol + rtol * K_ret``
            The default is zero (i.e. machine precision) for both.
            Note that for dualtree=True, rtol must be set to zero.
        dualtree : boolean  (default = False)
            if True, use the dual tree formalism.  This can be faster for
            large N, but can only be used for rtol = 0
        breadth_first : boolean (default = False)
            if True, use a breadth-first search.  If False (default) use a
            depth-first search.  Breadth-first is generally faster for
            compact kernels and/or high tolerances.

        Returns
        -------
        density : ndarray
            The array of density evaluations.  This has shape X.shape[:-1]

        Examples
        --------
        Compute a gaussian kernel density estimate:

        >>> import numpy as np
        >>> np.random.seed(1)
        >>> X = np.random.random((100, 3))
        >>> tree = BinaryTree(X)           # doctest: +SKIP
        >>> tree.kernel_density(X[:3], h=0.1, kernel='gaussian')
        array([ 6.94114649,  7.83281226,  7.2071716 ])
        """
        cdef DTYPE_t h_c = h
        cdef DTYPE_t atol_c = atol
        cdef DTYPE_t rtol_c = rtol
        cdef DTYPE_t min_bound, max_bound, dist_LB = 0, dist_UB = 0

        cdef ITYPE_t n_samples = self.data.shape[0]
        cdef ITYPE_t n_features = self.data.shape[1]
        cdef KernelType kernel_c

        # validate rtol
        if dualtree and rtol > 0:
            warnings.warn("rtol > 0 not compatible with dualtree")

        # validate kernel
        if kernel == 'gaussian':
            kernel_c = GAUSSIAN_KERNEL
        elif kernel == 'tophat':
            kernel_c = TOPHAT_KERNEL
        elif kernel == 'epanechnikov':
            kernel_c = EPANECHNIKOV_KERNEL
        elif kernel == 'exponential':
            kernel_c = EXPONENTIAL_KERNEL
        elif kernel == 'linear':
            kernel_c = LINEAR_KERNEL
        elif kernel == 'cosine':
            kernel_c = COSINE_KERNEL
        else:
            raise ValueError("kernel = '%s' not recognized" % kernel)

        # validate X and prepare for query
        X = array2d(X, dtype=DTYPE, order='C')

        if X.shape[-1] != n_features:
            raise ValueError("query data dimension must "
                             "match training data dimension")
        cdef DTYPE_t[:, ::1] Xarr = X.reshape((-1, n_features))
        
        cdef DTYPE_t[::1] density = np.zeros(Xarr.shape[0], dtype=DTYPE)

        cdef DTYPE_t* pt = &Xarr[0, 0]

        cdef NodeHeap nodeheap
        if breadth_first:
            nodeheap = NodeHeap(self.data.shape[0] // self.leaf_size)
        cdef DTYPE_t[::1] node_min_bounds
        cdef DTYPE_t[::1] node_max_bounds

        if dualtree:
            # dualtree algorithms assume atol is the tolerance per point
            other = self.__class__(Xarr, metric=self.dm,
                                   leaf_size=self.leaf_size)
            if breadth_first:
                self._kde_dual_breadthfirst(other, kernel_c, h_c,
                                            atol_c / self.data.shape[0],
                                            &density[0], nodeheap)
            else:
                self._kde_dual_depthfirst(0, other, 0, kernel_c,
                                          h_c, atol_c / self.data.shape[0],
                                          &density[0])
        else:
            if breadth_first:
                node_min_bounds = np.zeros(self.n_nodes)
                node_max_bounds = np.zeros(self.n_nodes)
                for i in range(Xarr.shape[0]):
                    density[i] = self._kde_single_breadthfirst(
                                            pt, kernel_c,
                                            h_c, atol_c,
                                            rtol_c, nodeheap,
                                            &node_min_bounds[0],
                                            &node_max_bounds[0])
                    pt += n_features
            else:
                for i in range(Xarr.shape[0]):
                    min_max_dist(self, 0, pt, &dist_LB, &dist_UB)
                    # compute max & min bounds on density within top node
                    min_bound = n_samples * compute_kernel(dist_UB,
                                                           h_c, kernel_c)
                    max_bound = n_samples * compute_kernel(dist_LB,
                                                           h_c, kernel_c)
                    self._kde_single_depthfirst(0, pt, kernel_c, h_c,
                                                atol_c, rtol_c,
                                                min_bound, max_bound,
                                                &min_bound, &max_bound)
                    density[i] = 0.5 * (min_bound + max_bound)
                    pt += n_features

        # normalize the results
        cdef DTYPE_t knorm = kernel_norm(h_c, kernel_c)
        for i in range(density.shape[0]):
            density[i] *= knorm

        return np.asarray(density).reshape(X.shape[:-1])

    def two_point_correlation(self, X, r, dualtree=False):
        """Compute the two-point correlation function

        Parameters
        ----------
        X : array_like
            An array of points to query.  Last dimension should match dimension
            of training data.
        r : array_like
            A one-dimensional array of distances
        dualtree : boolean (default = False)
            If true, use a dualtree algorithm.  Otherwise, use a single-tree
            algorithm.  Dual tree algorithms can have better scaling for
            large N.

        Returns
        -------
        counts : ndarray
            counts[i] contains the number of pairs of points with distance
            less than or equal to r[i]

        Examples
        --------
        Compute the two-point autocorrelation function of X:

        >>> import numpy as np
        >>> np.random.seed(0)
        >>> X = np.random.random((30, 3))
        >>> r = np.linspace(0, 1, 5)
        >>> tree = BinaryTree(X)     # doctest: +SKIP
        >>> tree.two_point_correlation(X, r)
        array([ 30,  62, 278, 580, 820])
        """
        cdef ITYPE_t n_features = self.data.shape[1]

        # validate X and prepare for query
        X = array2d(X, dtype=DTYPE, order='C')

        if X.shape[-1] != self.data.shape[1]:
            raise ValueError("query data dimension must "
                             "match training data dimension")

        cdef DTYPE_t[:, ::1] Xarr = X.reshape((-1, self.data.shape[1]))

        # prepare r for query
        r = np.asarray(r, dtype=DTYPE, order='C')
        r = np.atleast_1d(r).astype(DTYPE)
        if r.ndim != 1:
            raise ValueError("r must be a 1-dimensional array")
        i_rsort = np.argsort(r)
        cdef DTYPE_t[::1] rarr = r[i_rsort]

        # create array to hold counts
        count = np.zeros(r.shape[0], dtype=ITYPE)
        cdef ITYPE_t[::1] carr = count

        cdef DTYPE_t* pt = &Xarr[0, 0]

        if dualtree:
            other = self.__class__(Xarr, metric=self.dm,
                                   leaf_size=self.leaf_size)
            self._two_point_dual(0, other, 0, &rarr[0], &carr[0],
                                 0, rarr.shape[0])
        else:
            for i in range(Xarr.shape[0]):
                self._two_point_single(0, pt, &rarr[0], &carr[0],
                                       0, rarr.shape[0])
                pt += n_features
    
        return count

    cdef void _query_single_depthfirst(BinaryTree self, ITYPE_t i_node,
                                       DTYPE_t* pt, ITYPE_t i_pt,
                                       NeighborsHeap heap,
                                       DTYPE_t reduced_dist_LB):
        cdef NodeData_t node_info = self.node_data[i_node]

        cdef DTYPE_t dist_pt, reduced_dist_LB_1, reduced_dist_LB_2
        cdef ITYPE_t i, i1, i2
        
        cdef DTYPE_t* data = &self.data[0, 0]

        #------------------------------------------------------------
        # Case 1: query point is outside node radius:
        #         trim it from the query
        if reduced_dist_LB > heap.largest(i_pt):
            self.n_trims += 1

        #------------------------------------------------------------
        # Case 2: this is a leaf node.  Update set of nearby points
        elif node_info.is_leaf:
            self.n_leaves += 1
            for i in range(node_info.idx_start, node_info.idx_end):
                dist_pt = self.rdist(pt,
                                     &self.data[self.idx_array[i], 0],
                                     self.data.shape[1])
                if dist_pt < heap.largest(i_pt):
                    heap.push(i_pt, dist_pt, self.idx_array[i])

        #------------------------------------------------------------
        # Case 3: Node is not a leaf.  Recursively query subnodes
        #         starting with the closest
        else:
            self.n_splits += 1
            i1 = 2 * i_node + 1
            i2 = i1 + 1
            reduced_dist_LB_1 = min_rdist(self, i1, pt)
            reduced_dist_LB_2 = min_rdist(self, i2, pt)

            # recursively query subnodes
            if reduced_dist_LB_1 <= reduced_dist_LB_2:
                self._query_single_depthfirst(i1, pt, i_pt, heap,
                                              reduced_dist_LB_1)
                self._query_single_depthfirst(i2, pt, i_pt, heap,
                                              reduced_dist_LB_2)
            else:
                self._query_single_depthfirst(i2, pt, i_pt, heap,
                                              reduced_dist_LB_2)
                self._query_single_depthfirst(i1, pt, i_pt, heap,
                                              reduced_dist_LB_1)

    cdef void _query_single_breadthfirst(BinaryTree self, DTYPE_t* pt,
                                         ITYPE_t i_pt,
                                         NeighborsHeap heap,
                                         NodeHeap nodeheap):
        cdef ITYPE_t i, i_node
        cdef DTYPE_t dist_pt, reduced_dist_LB
        cdef NodeData_t* node_data = &self.node_data[0]
        cdef DTYPE_t* data = &self.data[0, 0]

        # Set up the node heap and push the head node onto it
        cdef NodeHeapData_t nodeheap_item
        nodeheap_item.val = min_rdist(self, 0, pt)
        nodeheap_item.i1 = 0
        nodeheap.push(nodeheap_item)

        while nodeheap.n > 0:
            nodeheap_item = nodeheap.pop()
            reduced_dist_LB = nodeheap_item.val
            i_node = nodeheap_item.i1
            node_info = node_data[i_node]

            #------------------------------------------------------------
            # Case 1: query point is outside node radius:
            #         trim it from the query
            if reduced_dist_LB > heap.largest(i_pt):
                self.n_trims += 1

            #------------------------------------------------------------
            # Case 2: this is a leaf node.  Update set of nearby points
            elif node_data[i_node].is_leaf:
                self.n_leaves += 1
                for i in range(node_data[i_node].idx_start,
                               node_data[i_node].idx_end):
                    dist_pt = self.rdist(pt,
                                         &self.data[self.idx_array[i], 0],
                                         self.data.shape[1])
                    if dist_pt < heap.largest(i_pt):
                        heap.push(i_pt, dist_pt, self.idx_array[i])

            #------------------------------------------------------------
            # Case 3: Node is not a leaf.  Add subnodes to the node heap
            else:
                self.n_splits += 1
                for i in range(2 * i_node + 1, 2 * i_node + 3):
                    nodeheap_item.i1 = i
                    nodeheap_item.val = min_rdist(self, i, pt)
                    nodeheap.push(nodeheap_item)

    cdef void _query_dual_depthfirst(BinaryTree self, ITYPE_t i_node1,
                                     BinaryTree other, ITYPE_t i_node2,
                                     DTYPE_t[::1] bounds,
                                     NeighborsHeap heap,
                                     DTYPE_t reduced_dist_LB):
        # note that the array `bounds` is maintained such that
        # bounds[i] is the largest distance among any of the
        # current neighbors in node i of the other tree.
        cdef NodeData_t node_info1 = self.node_data[i_node1]
        cdef NodeData_t node_info2 = other.node_data[i_node2]

        cdef DTYPE_t* data1 = &self.data[0, 0]
        cdef DTYPE_t* data2 = &other.data[0, 0]
        cdef ITYPE_t n_features = self.data.shape[1]

        cdef DTYPE_t bound_max, dist_pt, reduced_dist_LB1, reduced_dist_LB2
        cdef ITYPE_t i1, i2, i_pt, i_parent

        #------------------------------------------------------------
        # Case 1: nodes are further apart than the current bound:
        #         trim both from the query
        if reduced_dist_LB > bounds[i_node2]:
            pass

        #------------------------------------------------------------
        # Case 2: both nodes are leaves:
        #         do a brute-force search comparing all pairs
        elif node_info1.is_leaf and node_info2.is_leaf:
            bounds[i_node2] = 0

            for i2 in range(node_info2.idx_start, node_info2.idx_end):
                i_pt = other.idx_array[i2]
                
                if heap.largest(i_pt) <= reduced_dist_LB:
                    continue

                for i1 in range(node_info1.idx_start, node_info1.idx_end):
                    dist_pt = self.rdist(
                        data1 + n_features * self.idx_array[i1],
                        data2 + n_features * i_pt,
                        n_features)
                    if dist_pt < heap.largest(i_pt):
                        heap.push(i_pt, dist_pt, self.idx_array[i1])
                
                # keep track of node bound
                bounds[i_node2] = fmax(bounds[i_node2],
                                       heap.largest(i_pt))

            # update bounds up the tree
            while i_node2 > 0:
                i_parent = (i_node2 - 1) // 2
                bound_max = fmax(bounds[2 * i_parent + 1],
                                 bounds[2 * i_parent + 2])
                if bound_max < bounds[i_parent]:
                    bounds[i_parent] = bound_max
                    i_node2 = i_parent
                else:
                    break
            
        #------------------------------------------------------------
        # Case 3a: node 1 is a leaf: split node 2 and recursively
        #          query, starting with the nearest node
        elif node_info1.is_leaf:
            reduced_dist_LB1 = min_rdist_dual(self, i_node1,
                                              other, 2 * i_node2 + 1)
            reduced_dist_LB2 = min_rdist_dual(self, i_node1,
                                              other, 2 * i_node2 + 2)

            if reduced_dist_LB1 < reduced_dist_LB2:
                self._query_dual_depthfirst(i_node1, other, 2 * i_node2 + 1,
                                            bounds, heap, reduced_dist_LB1)
                self._query_dual_depthfirst(i_node1, other, 2 * i_node2 + 2,
                                            bounds, heap, reduced_dist_LB2)
            else:
                self._query_dual_depthfirst(i_node1, other, 2 * i_node2 + 2,
                                            bounds, heap, reduced_dist_LB2)
                self._query_dual_depthfirst(i_node1, other, 2 * i_node2 + 1,
                                            bounds, heap, reduced_dist_LB1)
            
        #------------------------------------------------------------
        # Case 3b: node 2 is a leaf: split node 1 and recursively
        #          query, starting with the nearest node
        elif node_info2.is_leaf:
            reduced_dist_LB1 = min_rdist_dual(self, 2 * i_node1 + 1,
                                              other, i_node2)
            reduced_dist_LB2 = min_rdist_dual(self, 2 * i_node1 + 2,
                                              other, i_node2)

            if reduced_dist_LB1 < reduced_dist_LB2:
                self._query_dual_depthfirst(2 * i_node1 + 1, other, i_node2,
                                            bounds, heap, reduced_dist_LB1)
                self._query_dual_depthfirst(2 * i_node1 + 2, other, i_node2,
                                            bounds, heap, reduced_dist_LB2)
            else:
                self._query_dual_depthfirst(2 * i_node1 + 2, other, i_node2,
                                            bounds, heap, reduced_dist_LB2)
                self._query_dual_depthfirst(2 * i_node1 + 1, other, i_node2,
                                            bounds, heap, reduced_dist_LB1)
        
        #------------------------------------------------------------
        # Case 4: neither node is a leaf:
        #         split both and recursively query all four pairs
        else:
            reduced_dist_LB1 = min_rdist_dual(self, 2 * i_node1 + 1,
                                              other, 2 * i_node2 + 1)
            reduced_dist_LB2 = min_rdist_dual(self, 2 * i_node1 + 2,
                                              other, 2 * i_node2 + 1)

            if reduced_dist_LB1 < reduced_dist_LB2:
                self._query_dual_depthfirst(2 * i_node1 + 1,
                                            other, 2 * i_node2 + 1,
                                            bounds, heap, reduced_dist_LB1)
                self._query_dual_depthfirst(2 * i_node1 + 2,
                                            other, 2 * i_node2 + 1,
                                            bounds, heap, reduced_dist_LB2)
            else:
                self._query_dual_depthfirst(2 * i_node1 + 2,
                                            other, 2 * i_node2 + 1,
                                            bounds, heap, reduced_dist_LB2)
                self._query_dual_depthfirst(2 * i_node1 + 1,
                                            other, 2 * i_node2 + 1,
                                            bounds, heap, reduced_dist_LB1)

            reduced_dist_LB1 = min_rdist_dual(self, 2 * i_node1 + 1,
                                              other, 2 * i_node2 + 2)
            reduced_dist_LB2 = min_rdist_dual(self, 2 * i_node1 + 2,
                                              other, 2 * i_node2 + 2)
            if reduced_dist_LB1 < reduced_dist_LB2:
                self._query_dual_depthfirst(2 * i_node1 + 1,
                                            other, 2 * i_node2 + 2,
                                            bounds, heap, reduced_dist_LB1)
                self._query_dual_depthfirst(2 * i_node1 + 2,
                                            other, 2 * i_node2 + 2,
                                            bounds, heap, reduced_dist_LB2)
            else:
                self._query_dual_depthfirst(2 * i_node1 + 2,
                                            other, 2 * i_node2 + 2,
                                            bounds, heap, reduced_dist_LB2)
                self._query_dual_depthfirst(2 * i_node1 + 1,
                                            other, 2 * i_node2 + 2,
                                            bounds, heap, reduced_dist_LB1)

    cdef void _query_dual_breadthfirst(BinaryTree self, BinaryTree other,
                                       NeighborsHeap heap, NodeHeap nodeheap):
        cdef ITYPE_t i, i1, i2, i_node1, i_node2, i_pt
        cdef DTYPE_t dist_pt, reduced_dist_LB
        cdef DTYPE_t[::1] bounds = np.inf + np.zeros(other.node_data.shape[0])
        cdef NodeData_t* node_data1 = &self.node_data[0]
        cdef NodeData_t* node_data2 = &other.node_data[0]
        cdef NodeData_t node_info1, node_info2
        cdef DTYPE_t* data1 = &self.data[0, 0]
        cdef DTYPE_t* data2 = &other.data[0, 0]
        cdef ITYPE_t n_features = self.data.shape[1]

        # Set up the node heap and push the head nodes onto it
        cdef NodeHeapData_t nodeheap_item
        nodeheap_item.val = min_rdist_dual(self, 0, other, 0)
        nodeheap_item.i1 = 0
        nodeheap_item.i2 = 0
        nodeheap.push(nodeheap_item)

        while nodeheap.n > 0:
            nodeheap_item = nodeheap.pop()
            reduced_dist_LB = nodeheap_item.val
            i_node1 = nodeheap_item.i1
            i_node2 = nodeheap_item.i2
            
            node_info1 = node_data1[i_node1]
            node_info2 = node_data2[i_node2]

            #------------------------------------------------------------
            # Case 1: nodes are further apart than the current bound:
            #         trim both from the query
            if reduced_dist_LB > bounds[i_node2]:
                pass

            #------------------------------------------------------------
            # Case 2: both nodes are leaves:
            #         do a brute-force search comparing all pairs
            elif node_info1.is_leaf and node_info2.is_leaf:
                bounds[i_node2] = -1

                for i2 in range(node_info2.idx_start, node_info2.idx_end):
                    i_pt = other.idx_array[i2]
                
                    if heap.largest(i_pt) <= reduced_dist_LB:
                        continue

                    for i1 in range(node_info1.idx_start, node_info1.idx_end):
                        dist_pt = self.rdist(
                            data1 + n_features * self.idx_array[i1],
                            data2 + n_features * i_pt,
                            n_features)
                        if dist_pt < heap.largest(i_pt):
                            heap.push(i_pt, dist_pt, self.idx_array[i1])
                
                    # keep track of node bound
                    bounds[i_node2] = fmax(bounds[i_node2],
                                           heap.largest(i_pt))
            
            #------------------------------------------------------------
            # Case 3a: node 1 is a leaf: split node 2 and recursively
            #          query, starting with the nearest node
            elif node_info1.is_leaf:
                nodeheap_item.i1 = i_node1
                for i2 in range(2 * i_node2 + 1, 2 * i_node2 + 3):
                    nodeheap_item.i2 = i2
                    nodeheap_item.val = min_rdist_dual(self, i_node1, other, i2)
                    nodeheap.push(nodeheap_item)
            
            #------------------------------------------------------------
            # Case 3b: node 2 is a leaf: split node 1 and recursively
            #          query, starting with the nearest node
            elif node_info2.is_leaf:
                nodeheap_item.i2 = i_node2
                for i1 in range(2 * i_node1 + 1, 2 * i_node1 + 3):
                    nodeheap_item.i1 = i1
                    nodeheap_item.val = min_rdist_dual(self, i1, other, i_node2)
                    nodeheap.push(nodeheap_item)
        
            #------------------------------------------------------------
            # Case 4: neither node is a leaf:
            #         split both and recursively query all four pairs
            else:
                for i1 in range(2 * i_node1 + 1, 2 * i_node1 + 3):
                    for i2 in range(2 * i_node2 + 1, 2 * i_node2 + 3):
                        nodeheap_item.i1 = i1
                        nodeheap_item.i2 = i2
                        nodeheap_item.val = min_rdist_dual(self, i1, other, i2)
                        nodeheap.push(nodeheap_item)

    cdef ITYPE_t _query_radius_single(BinaryTree self,
                                      ITYPE_t i_node,
                                      DTYPE_t* pt, DTYPE_t r,
                                      ITYPE_t* indices,
                                      DTYPE_t* distances,
                                      ITYPE_t count,
                                      int count_only,
                                      int return_distance):
        cdef DTYPE_t* data = &self.data[0, 0]
        cdef ITYPE_t* idx_array = &self.idx_array[0]
        cdef ITYPE_t n_features = self.data.shape[1]
        cdef NodeData_t node_info = self.node_data[i_node]

        cdef ITYPE_t i
        cdef DTYPE_t reduced_r

        cdef DTYPE_t dist_pt, dist_LB = 0, dist_UB = 0
        min_max_dist(self, i_node, pt, &dist_LB, &dist_UB)

        #------------------------------------------------------------
        # Case 1: all node points are outside distance r.
        #         prune this branch.
        if dist_LB > r:
            pass

        #------------------------------------------------------------
        # Case 2: all node points are within distance r
        #         add all points to neighbors
        elif dist_UB <= r:
            if count_only:
                count += (node_info.idx_end - node_info.idx_start)
            else:
                for i in range(node_info.idx_start, node_info.idx_end):
                    if (count < 0) or (count >= self.data.shape[0]):
                        raise ValueError("Fatal: count too big: "
                                         "this should never happen")
                    indices[count] = idx_array[i]
                    if return_distance:
                        distances[count] = self.dist(pt, (data + n_features
                                                          * idx_array[i]),
                                                     n_features)
                    count += 1

        #------------------------------------------------------------
        # Case 3: this is a leaf node.  Go through all points to
        #         determine if they fall within radius
        elif node_info.is_leaf:
            reduced_r = self.dm.dist_to_rdist(r)

            for i in range(node_info.idx_start, node_info.idx_end):
                dist_pt = self.rdist(pt, (data + n_features * idx_array[i]),
                                     n_features)
                if dist_pt <= reduced_r:
                    if (count < 0) or (count >= self.data.shape[0]):
                        raise ValueError("Fatal: count out of range. "
                                         "This should never happen.")
                    if count_only:
                        pass
                    else:
                        indices[count] = idx_array[i]
                        if return_distance:
                            distances[count] = self.dm.rdist_to_dist(dist_pt)
                    count += 1

        #------------------------------------------------------------
        # Case 4: Node is not a leaf.  Recursively query subnodes
        else:
            count = self._query_radius_single(2 * i_node + 1, pt, r,
                                              indices, distances, count,
                                              count_only, return_distance)
            count = self._query_radius_single(2 * i_node + 2, pt, r,
                                              indices, distances, count,
                                              count_only, return_distance)

        return count

    cdef DTYPE_t _kde_single_breadthfirst(self, DTYPE_t* pt,
                                          KernelType kernel, DTYPE_t h,
                                          DTYPE_t atol, DTYPE_t rtol,
                                          NodeHeap nodeheap,
                                          DTYPE_t* node_min_bounds,
                                          DTYPE_t* node_max_bounds):
        cdef ITYPE_t i, i1, i2, N1, N2, i_node
        cdef DTYPE_t global_min_bound, global_max_bound

        cdef DTYPE_t* data = &self.data[0, 0]
        cdef ITYPE_t* idx_array = &self.idx_array[0]
        cdef NodeData_t* node_data = &self.node_data[0]
        cdef ITYPE_t N = self.data.shape[0]
        cdef ITYPE_t n_features = self.data.shape[1]
        cdef DTYPE_t knorm = kernel_norm(h, kernel)

        cdef NodeData_t node_info
        cdef DTYPE_t dist_pt, dens_contribution
        cdef DTYPE_t dist_LB_1 = 0, dist_LB_2 = 0
        cdef DTYPE_t dist_UB_1 = 0, dist_UB_2 = 0

        cdef DTYPE_t dist_UB, dist_LB

        # push the top node to the heap
        cdef NodeHeapData_t nodeheap_item
        nodeheap_item.val = min_dist(self, 0, pt)
        nodeheap_item.i1 = 0
        nodeheap.push(nodeheap_item)

        global_min_bound = N * compute_kernel(max_dist(self, 0, pt),
                                              h, kernel)
        global_max_bound = N * compute_kernel(nodeheap_item.val,
                                              h, kernel)

        node_min_bounds[0] = global_min_bound
        node_max_bounds[0] = global_max_bound

        while nodeheap.n > 0:
            nodeheap_item = nodeheap.pop()
            i_node = nodeheap_item.i1
            node_info = node_data[i_node]
            N1 = node_info.idx_end - node_info.idx_start

            #------------------------------------------------------------
            # Case 1: local bounds are equal to within per-point tolerance.
            if (knorm * (node_max_bounds[i_node] - node_min_bounds[i_node])
                <= (atol + rtol * knorm * node_min_bounds[i_node]) * N1 / N):
                pass

            #------------------------------------------------------------
            # Case 2: global bounds are within rtol & atol.
            elif (knorm * (global_max_bound - global_min_bound)
                  <= (atol + rtol * knorm * global_min_bound)):
                break

            #------------------------------------------------------------
            # Case 3: node is a leaf. Count contributions from all points
            elif node_info.is_leaf:
                global_min_bound -= node_min_bounds[i_node]
                global_max_bound -= node_max_bounds[i_node]
                for i in range(node_info.idx_start, node_info.idx_end):
                    dist_pt = self.dist(pt, data + n_features * idx_array[i],
                                        n_features)
                    dens_contribution = compute_kernel(dist_pt, h, kernel)
                    global_min_bound += dens_contribution
                    global_max_bound += dens_contribution

            #------------------------------------------------------------
            # Case 4: split node and query subnodes
            else:
                i1 = 2 * i_node + 1
                i2 = 2 * i_node + 2

                N1 = node_data[i1].idx_end - node_data[i1].idx_start
                N2 = node_data[i2].idx_end - node_data[i2].idx_start

                min_max_dist(self, i1, pt, &dist_LB_1, &dist_UB_1)
                min_max_dist(self, i2, pt, &dist_LB_2, &dist_UB_2)

                node_max_bounds[i1] = N1 * compute_kernel(dist_LB_1, h, kernel)
                node_min_bounds[i1] = N1 * compute_kernel(dist_UB_1, h, kernel)

                node_max_bounds[i2] = N2 * compute_kernel(dist_LB_2, h, kernel)
                node_min_bounds[i2] = N2 * compute_kernel(dist_UB_2, h, kernel)
            
                global_min_bound += (node_min_bounds[i1] + node_min_bounds[i2]
                                     - node_min_bounds[i_node])
                global_max_bound += (node_max_bounds[i1] + node_max_bounds[i2]
                                     - node_max_bounds[i_node])

                nodeheap_item.val = dist_LB_1
                nodeheap_item.i1 = i1
                nodeheap.push(nodeheap_item)

                nodeheap_item.val = dist_LB_2
                nodeheap_item.i1 = i2
                nodeheap.push(nodeheap_item)

        nodeheap.clear()
        return 0.5 * (global_max_bound + global_min_bound)


    cdef void _kde_single_depthfirst(self, ITYPE_t i_node, DTYPE_t* pt,
                                     KernelType kernel, DTYPE_t h,
                                     DTYPE_t atol, DTYPE_t rtol,
                                     DTYPE_t local_min_bound,
                                     DTYPE_t local_max_bound,
                                     DTYPE_t* global_min_bound,
                                     DTYPE_t* global_max_bound):
        cdef ITYPE_t i, i1, i2, N1, N2

        cdef DTYPE_t* data = &self.data[0, 0]
        cdef ITYPE_t* idx_array = &self.idx_array[0]
        cdef ITYPE_t n_features = self.data.shape[1]
        cdef DTYPE_t knorm = kernel_norm(h, kernel)

        cdef NodeData_t node_info = self.node_data[i_node]
        cdef DTYPE_t dist_pt, dens_contribution

        cdef DTYPE_t child1_min_bound, child2_min_bound, dist_UB = 0
        cdef DTYPE_t child1_max_bound, child2_max_bound, dist_LB = 0

        N1 = node_info.idx_end - node_info.idx_start
        N2 = self.data.shape[0]

        #------------------------------------------------------------
        # Case 1: local bounds are equal to within errors.  Return
        if ((knorm * local_min_bound) >=
            (knorm * local_max_bound
             - (atol + rtol * knorm * local_min_bound) * N1 / N2)):
            pass

        #------------------------------------------------------------
        # Case 2: global bounds are within rtol & atol. Return
        elif (knorm * (global_max_bound[0] - global_min_bound[0])
            <= (atol + rtol * knorm * global_min_bound[0])):
            pass

        #------------------------------------------------------------
        # Case 3: node is a leaf. Count contributions from all points
        elif node_info.is_leaf:
            global_min_bound[0] -= local_min_bound
            global_max_bound[0] -= local_max_bound
            for i in range(node_info.idx_start, node_info.idx_end):
                dist_pt = self.dist(pt, (data + n_features * idx_array[i]),
                                    n_features)
                dens_contribution = compute_kernel(dist_pt, h, kernel)
                global_min_bound[0] += dens_contribution
                global_max_bound[0] += dens_contribution

        #------------------------------------------------------------
        # Case 4: split node and query subnodes
        else:
            i1 = 2 * i_node + 1
            i2 = 2 * i_node + 2

            N1 = self.node_data[i1].idx_end - self.node_data[i1].idx_start
            N2 = self.node_data[i2].idx_end - self.node_data[i2].idx_start

            min_max_dist(self, i1, pt, &dist_LB, &dist_UB)
            child1_min_bound = N1 * compute_kernel(dist_UB, h, kernel)
            child1_max_bound = N1 * compute_kernel(dist_LB, h, kernel)

            min_max_dist(self, i2, pt, &dist_LB, &dist_UB)
            child2_min_bound = N2 * compute_kernel(dist_UB, h, kernel)
            child2_max_bound = N2 * compute_kernel(dist_LB, h, kernel)
            
            global_min_bound[0] += (child1_min_bound + child2_min_bound
                                    - local_min_bound)
            global_max_bound[0] += (child1_max_bound + child2_max_bound
                                    - local_max_bound)

            self._kde_single_depthfirst(i1, pt, kernel, h,
                                        atol, rtol,
                                        child1_min_bound, child1_max_bound,
                                        global_min_bound, global_max_bound)
            self._kde_single_depthfirst(i2, pt, kernel, h,
                                        atol, rtol,
                                        child2_min_bound, child2_max_bound,
                                        global_min_bound, global_max_bound)


    cdef void _kde_dual_breadthfirst(BinaryTree self, BinaryTree other,
                                     KernelType kernel,
                                     DTYPE_t h, DTYPE_t atol,
                                     DTYPE_t* density, NodeHeap nodeheap):
        # note that atol here is absolute tolerance *per point*
        cdef ITYPE_t i1, i2, i_node1, i_node2

        cdef DTYPE_t* data1 = &self.data[0, 0]
        cdef DTYPE_t* data2 = &other.data[0, 0]

        cdef ITYPE_t* idx_array1 = &self.idx_array[0]
        cdef ITYPE_t* idx_array2 = &other.idx_array[0]

        cdef ITYPE_t n_features = self.data.shape[1]

        cdef NodeData_t node_info1, node_info2

        cdef NodeData_t* node_data1 = &self.node_data[0]
        cdef NodeData_t* node_data2 = &other.node_data[0]

        cdef NodeHeapData_t nodeheap_item

        cdef DTYPE_t dist_LB, dist_UB, dens_LB, dens_UB
        
        nodeheap_item.val = min_dist_dual(self, 0, other, 0)
        nodeheap_item.i1 = 0
        nodeheap_item.i2 = 0
        nodeheap.push(nodeheap_item)

        cdef DTYPE_t knorm = kernel_norm(h, kernel)

        while nodeheap.n > 0:
            nodeheap_item = nodeheap.pop()
            dist_LB = nodeheap_item.val
            i_node1 = nodeheap_item.i1
            i_node2 = nodeheap_item.i2
            node_info1 = node_data1[i_node1]
            node_info2 = node_data1[i_node2]

            dist_UB = max_dist_dual(self, i_node1, other, i_node2)
            dens_LB = compute_kernel(dist_UB, h, kernel)
            dens_UB = compute_kernel(dist_LB, h, kernel)

            #------------------------------------------------------------
            # Case 1: nodes are within desired tolerance
            if dens_LB * knorm + atol >= dens_UB * knorm:
                dens_LB *= (node_info1.idx_end - node_info1.idx_start)
                for i2 in range(node_info2.idx_start, node_info2.idx_end):
                    density[idx_array2[i2]] += dens_LB

            #------------------------------------------------------------
            # Case 2: both nodes are leaves: go through all pairs
            elif node_info1.is_leaf and node_info2.is_leaf:
                for i2 in range(node_info2.idx_start, node_info2.idx_end):
                    for i1 in range(node_info1.idx_start, node_info1.idx_end):
                        dist_pt = self.dist(data1 + n_features * idx_array1[i1],
                                            data2 + n_features * idx_array2[i2],
                                            n_features)
                        density[idx_array2[i2]] += compute_kernel(dist_pt, h,
                                                                  kernel)

            #------------------------------------------------------------
            # Case 3a: only one node is a leaf: split the other
            elif node_info1.is_leaf:
                for i2 in range(2 * i_node2 + 1, 2 * i_node2 + 3):
                    nodeheap_item.val = min_dist_dual(self, i_node1, other, i2)
                    nodeheap_item.i1 = i_node1
                    nodeheap_item.i2 = i2
                    nodeheap.push(nodeheap_item)

            elif node_info2.is_leaf:
                for i1 in range(2 * i_node1 + 1, 2 * i_node1 + 3):
                    nodeheap_item.val = min_dist_dual(self, i1, other, i_node2)
                    nodeheap_item.i1 = i1
                    nodeheap_item.i2 = i_node2
                    nodeheap.push(nodeheap_item)

            #------------------------------------------------------------
            # Case 3b: both nodes need to be split
            else:
                for i1 in range(2 * i_node1 + 1, 2 * i_node1 + 3):
                    for i2 in range(2 * i_node2 + 1, 2 * i_node2 + 3):
                        nodeheap_item.val = min_dist_dual(self, i1, other, i2)
                        nodeheap_item.i1 = i1
                        nodeheap_item.i2 = i2
                        nodeheap.push(nodeheap_item)


    cdef void _kde_dual_depthfirst(BinaryTree self, ITYPE_t i_node1,
                                   BinaryTree other, ITYPE_t i_node2,
                                   KernelType kernel,
                                   DTYPE_t h, DTYPE_t atol,
                                   DTYPE_t* density):
        # note that atol here is absolute tolerance *per point*
        cdef ITYPE_t i1, i2

        cdef DTYPE_t* data1 = &self.data[0, 0]
        cdef DTYPE_t* data2 = &other.data[0, 0]

        cdef ITYPE_t* idx_array1 = &self.idx_array[0]
        cdef ITYPE_t* idx_array2 = &other.idx_array[0]

        cdef ITYPE_t n_features = self.data.shape[1]

        cdef NodeData_t node_info1 = self.node_data[i_node1]
        cdef NodeData_t node_info2 = other.node_data[i_node2]
        cdef DTYPE_t dist_pt

        # XXX: for efficiency, calculate dist_LB and dist_UB at the same time
        cdef DTYPE_t dist_LB = min_dist_dual(self, i_node1, other, i_node2)
        cdef DTYPE_t dist_UB = max_dist_dual(self, i_node1, other, i_node2)

        cdef DTYPE_t dens_UB = compute_kernel(dist_LB, h, kernel)
        cdef DTYPE_t dens_LB = compute_kernel(dist_UB, h, kernel)
        cdef DTYPE_t knorm = kernel_norm(h, kernel)

        #------------------------------------------------------------
        # Case 1: points are all far enough that contribution is zero
        #         to within the desired tolerance
        if dens_LB * knorm + atol >= dens_UB * knorm:
            dens_LB *= (node_info1.idx_end - node_info1.idx_start)
            for i2 in range(node_info2.idx_start, node_info2.idx_end):
                density[idx_array2[i2]] += dens_LB

        #------------------------------------------------------------
        # Case 2: both nodes are leaves: go through all pairs
        elif node_info1.is_leaf and node_info2.is_leaf:
            for i2 in range(node_info2.idx_start, node_info2.idx_end):
                for i1 in range(node_info1.idx_start, node_info1.idx_end):
                    dist_pt = self.dist(data1 + n_features * idx_array1[i1],
                                        data2 + n_features * idx_array2[i2],
                                        n_features)
                    density[idx_array2[i2]] += compute_kernel(dist_pt, h,
                                                              kernel)

        #------------------------------------------------------------
        # Case 3a: only one node is a leaf: split the other
        elif node_info1.is_leaf:
            for i2 in range(2 * i_node2 + 1, 2 * i_node2 + 3):
                self._kde_dual_depthfirst(i_node1, other, i2, kernel,
                                          h, atol, density)

        elif node_info2.is_leaf:
            for i1 in range(2 * i_node1 + 1, 2 * i_node1 + 3):
                self._kde_dual_depthfirst(i1, other, i_node2, kernel,
                                          h, atol, density)

        #------------------------------------------------------------
        # Case 3b: both nodes need to be split
        else:
            for i1 in range(2 * i_node1 + 1, 2 * i_node1 + 3):
                for i2 in range(2 * i_node2 + 1, 2 * i_node2 + 3):
                    self._kde_dual_depthfirst(i1, other, i2, kernel,
                                              h, atol, density)

    cdef void _two_point_single(self, ITYPE_t i_node, DTYPE_t* pt, DTYPE_t* r,
                                ITYPE_t* count, ITYPE_t i_min, ITYPE_t i_max):
        cdef DTYPE_t* data = &self.data[0, 0]
        cdef ITYPE_t* idx_array = &self.idx_array[0]
        cdef ITYPE_t n_features = self.data.shape[1]
        cdef NodeData_t node_info = self.node_data[i_node]

        cdef ITYPE_t i, j, Npts
        cdef DTYPE_t reduced_r

        cdef DTYPE_t dist_pt, dist_LB = 0, dist_UB = 0
        min_max_dist(self, i_node, pt, &dist_LB, &dist_UB)

        #------------------------------------------------------------
        # Go through bounds and check for cuts
        while i_min < i_max:
            if dist_LB > r[i_min]:
                i_min += 1
            else:
                break

        while i_max > i_min:
            Npts = (node_info.idx_end - node_info.idx_start)
            if dist_UB <= r[i_max - 1]:
                count[i_max - 1] += Npts
                i_max -= 1
            else:
                break

        if i_min < i_max:
            # If node is a leaf, go through all points
            if node_info.is_leaf:
                for i in range(node_info.idx_start, node_info.idx_end):
                    dist_pt = self.dist(pt, (data + n_features * idx_array[i]),
                                        n_features)
                    j = i_max - 1
                    while (j >= i_min) and (dist_pt <= r[j]):
                        count[j] += 1
                        j -= 1

            else:
                self._two_point_single(2 * i_node + 1, pt, r,
                                       count, i_min, i_max)
                self._two_point_single(2 * i_node + 2, pt, r,
                                       count, i_min, i_max)

    cdef void _two_point_dual(self, ITYPE_t i_node1,
                              BinaryTree other, ITYPE_t i_node2,
                              DTYPE_t* r, ITYPE_t* count,
                              ITYPE_t i_min, ITYPE_t i_max):
        cdef DTYPE_t* data1 = &self.data[0, 0]
        cdef DTYPE_t* data2 = &other.data[0, 0]
        cdef ITYPE_t* idx_array1 = &self.idx_array[0]
        cdef ITYPE_t* idx_array2 = &other.idx_array[0]
        cdef NodeData_t node_info1 = self.node_data[i_node1]
        cdef NodeData_t node_info2 = other.node_data[i_node2]

        cdef ITYPE_t n_features = self.data.shape[1]

        cdef ITYPE_t i1, i2, j, Npts
        cdef DTYPE_t reduced_r

        cdef DTYPE_t dist_pt, dist_LB = 0, dist_UB = 0
        dist_LB = min_dist_dual(self, i_node1, other, i_node2)
        dist_UB = max_dist_dual(self, i_node1, other, i_node2)

        #------------------------------------------------------------
        # Go through bounds and check for cuts
        while i_min < i_max:
            if dist_LB > r[i_min]:
                i_min += 1
            else:
                break

        while i_max > i_min:
            Npts = ((node_info1.idx_end - node_info1.idx_start)
                    * (node_info2.idx_end - node_info2.idx_start))
            if dist_UB <= r[i_max - 1]:
                count[i_max - 1] += Npts
                i_max -= 1
            else:
                break

        if i_min < i_max:
            if node_info1.is_leaf and node_info2.is_leaf:
                # If both nodes are leaves, go through all points
                for i1 in range(node_info1.idx_start, node_info1.idx_end):
                    for i2 in range(node_info2.idx_start, node_info2.idx_end):
                        dist_pt = self.dist((data1 + n_features
                                             * idx_array1[i1]),
                                            (data2 + n_features
                                             * idx_array2[i2]),
                                            n_features)
                        j = i_max - 1
                        while (j >= i_min) and (dist_pt <= r[j]):
                            count[j] += 1
                            j -= 1

            elif node_info1.is_leaf:
                # If only one is a leaf, split the other
                for i2 in range(2 * i_node2 + 1, 2 * i_node2 + 3):
                    self._two_point_dual(i_node1, other, i2,
                                         r, count, i_min, i_max)

            elif node_info2.is_leaf:
                for i1 in range(2 * i_node1 + 1, 2 * i_node1 + 3):
                    self._two_point_dual(i1, other, i_node2,
                                         r, count, i_min, i_max)

            else:
                 # neither is a leaf: split & query both
                for i1 in range(2 * i_node1 + 1, 2 * i_node1 + 3):
                    for i2 in range(2 * i_node2 + 1, 2 * i_node2 + 3):
                        self._two_point_dual(i1, other, i2,
                                             r, count, i_min, i_max)

######################################################################
# Python functions for benchmarking and testing

def load_heap(DTYPE_t[:, ::1] X, ITYPE_t k):
    """test fully loading the heap"""
    assert k <= X.shape[1]
    cdef NeighborsHeap heap = NeighborsHeap(X.shape[0], k)
    cdef ITYPE_t i, j
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            heap.push(i, X[i, j], j)
    return heap.get_arrays()


def simultaneous_sort(DTYPE_t[:, ::1] distances, ITYPE_t[:, ::1] indices):
    """In-place simultaneous sort the given row of the arrays
    
    This python wrapper exists primarily to enable unit testing
    of the _simultaneous_sort C routine.
    """
    assert distances.shape[0] == indices.shape[0]
    assert distances.shape[1] == indices.shape[1]
    cdef ITYPE_t row
    for row in range(distances.shape[0]):
        _simultaneous_sort(&distances[row, 0],
                           &indices[row, 0],
                           distances.shape[1])

def nodeheap_sort(DTYPE_t[::1] vals):
    """In-place reverse sort of vals using NodeHeap"""
    cdef ITYPE_t[::1] indices = np.zeros(vals.shape[0], dtype=ITYPE)
    cdef DTYPE_t[::1] vals_sorted = np.zeros_like(vals)

    # use initial size 0 to check corner case
    cdef NodeHeap heap = NodeHeap(0)
    cdef NodeHeapData_t data
    cdef ITYPE_t i
    for i in range(vals.shape[0]):
        data.val = vals[i]
        data.i1 = i
        data.i2 = i + 1
        heap.push(data)

    for i in range(vals.shape[0]):
        data = heap.pop()
        vals_sorted[i] = data.val
        indices[i] = data.i1

    return np.asarray(vals_sorted), np.asarray(indices)
