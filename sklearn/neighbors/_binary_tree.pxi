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
#
# Implementation Notes
# --------------------
# This implementation uses the common object-oriented approach of having an
# abstract base class which is extended by the KDTree and BallTree
# specializations.
#
# The BinaryTree "base class" is defined here and then subclassed in the BallTree
# and KDTree pyx files. These files include implementations of the
# "abstract" methods.

# Necessary Helper Functions
# --------------------------
# These are the names and descriptions of the "abstract" functions which are
# defined in kd_tree.pyx and ball_tree.pyx:

# cdef int allocate_data(BinaryTree tree, ITYPE_t n_nodes, ITYPE_t n_features):
#     """Allocate arrays needed for the KD Tree"""

# cdef int init_node(BinaryTree tree, ITYPE_t i_node,
#                    ITYPE_t idx_start, ITYPE_t idx_end):
#    """Initialize the node for the dataset stored in tree.data"""

# cdef DTYPE_t min_rdist(BinaryTree tree, ITYPE_t i_node, DTYPE_t* pt):
#     """Compute the minimum reduced-distance between a point and a node"""

# cdef DTYPE_t min_dist(BinaryTree tree, ITYPE_t i_node, DTYPE_t* pt):
#     """Compute the minimum distance between a point and a node"""

# cdef DTYPE_t max_rdist(BinaryTree tree, ITYPE_t i_node, DTYPE_t* pt):
#     """Compute the maximum reduced-distance between a point and a node"""

# cdef DTYPE_t max_dist(BinaryTree tree, ITYPE_t i_node, DTYPE_t* pt):
#     """Compute the maximum distance between a point and a node"""

# cdef inline int min_max_dist(BinaryTree tree, ITYPE_t i_node, DTYPE_t* pt,
#                              DTYPE_t* min_dist, DTYPE_t* max_dist):
#     """Compute the minimum and maximum distance between a point and a node"""

# cdef inline DTYPE_t min_rdist_dual(BinaryTree tree1, ITYPE_t i_node1,
#                                    BinaryTree tree2, ITYPE_t i_node2):
#     """Compute the minimum reduced distance between two nodes"""

# cdef inline DTYPE_t min_dist_dual(BinaryTree tree1, ITYPE_t i_node1,
#                                   BinaryTree tree2, ITYPE_t i_node2):
#     """Compute the minimum distance between two nodes"""

# cdef inline DTYPE_t max_rdist_dual(BinaryTree tree1, ITYPE_t i_node1,
#                                    BinaryTree tree2, ITYPE_t i_node2):
#     """Compute the maximum reduced distance between two nodes"""

# cdef inline DTYPE_t max_dist_dual(BinaryTree tree1, ITYPE_t i_node1,
#                                   BinaryTree tree2, ITYPE_t i_node2):
#     """Compute the maximum distance between two nodes"""

cimport cython
cimport numpy as np
from libc.math cimport fabs, sqrt, exp, cos, pow, log, lgamma
from libc.math cimport fmin, fmax
from libc.stdlib cimport calloc, malloc, free
from libc.string cimport memcpy

import numpy as np
import warnings
from ..utils import check_array

from ._typedefs cimport DTYPE_t, ITYPE_t, DITYPE_t
from ._typedefs import DTYPE, ITYPE

from ._dist_metrics cimport (DistanceMetric, euclidean_dist, euclidean_rdist,
                             euclidean_dist_to_rdist, euclidean_rdist_to_dist)

cdef extern from "numpy/arrayobject.h":
    void PyArray_ENABLEFLAGS(np.ndarray arr, int flags)

np.import_array()

# some handy constants
cdef DTYPE_t INF = np.inf
cdef DTYPE_t NEG_INF = -np.inf
cdef DTYPE_t PI = np.pi
cdef DTYPE_t ROOT_2PI = sqrt(2 * PI)
cdef DTYPE_t LOG_PI = log(PI)
cdef DTYPE_t LOG_2PI = log(2 * PI)


# Some compound datatypes used below:
cdef struct NodeHeapData_t:
    DTYPE_t val
    ITYPE_t i1
    ITYPE_t i2

# build the corresponding numpy dtype for NodeHeapData
cdef NodeHeapData_t nhd_tmp
NodeHeapData = np.asarray(<NodeHeapData_t[:1]>(&nhd_tmp)).dtype

cdef struct NodeData_t:
    ITYPE_t idx_start
    ITYPE_t idx_end
    ITYPE_t is_leaf
    DTYPE_t radius

# build the corresponding numpy dtype for NodeData
cdef NodeData_t nd_tmp
NodeData = np.asarray(<NodeData_t[:1]>(&nd_tmp)).dtype


######################################################################
# Numpy 1.3-1.4 compatibility utilities
cdef DTYPE_t[::1] get_memview_DTYPE_1D(
                               np.ndarray[DTYPE_t, ndim=1, mode='c'] X):
    return <DTYPE_t[:X.shape[0]:1]> (<DTYPE_t*> X.data)


cdef DTYPE_t[:, ::1] get_memview_DTYPE_2D(
                               np.ndarray[DTYPE_t, ndim=2, mode='c'] X):
    return <DTYPE_t[:X.shape[0], :X.shape[1]:1]> (<DTYPE_t*> X.data)


cdef DTYPE_t[:, :, ::1] get_memview_DTYPE_3D(
                               np.ndarray[DTYPE_t, ndim=3, mode='c'] X):
    return <DTYPE_t[:X.shape[0], :X.shape[1], :X.shape[2]:1]>\
                                                       (<DTYPE_t*> X.data)


cdef ITYPE_t[::1] get_memview_ITYPE_1D(
                               np.ndarray[ITYPE_t, ndim=1, mode='c'] X):
    return <ITYPE_t[:X.shape[0]:1]> (<ITYPE_t*> X.data)


cdef ITYPE_t[:, ::1] get_memview_ITYPE_2D(
                               np.ndarray[ITYPE_t, ndim=2, mode='c'] X):
    return <ITYPE_t[:X.shape[0], :X.shape[1]:1]> (<ITYPE_t*> X.data)


cdef NodeHeapData_t[::1] get_memview_NodeHeapData_1D(
                    np.ndarray[NodeHeapData_t, ndim=1, mode='c'] X):
    return <NodeHeapData_t[:X.shape[0]:1]> (<NodeHeapData_t*> X.data)


cdef NodeData_t[::1] get_memview_NodeData_1D(
                    np.ndarray[NodeData_t, ndim=1, mode='c'] X):
    return <NodeData_t[:X.shape[0]:1]> (<NodeData_t*> X.data)

######################################################################



######################################################################
# Define doc strings, substituting the appropriate class name using
# the DOC_DICT variable defined in the pyx files.
CLASS_DOC = \
"""
{BinaryTree}(X, leaf_size=40, metric='minkowski', **kwargs)

{BinaryTree} for fast generalized N-point problems

Parameters
----------
X : array-like of shape (n_samples, n_features)
    n_samples is the number of points in the data set, and
    n_features is the dimension of the parameter space.
    Note: if X is a C-contiguous array of doubles then data will
    not be copied. Otherwise, an internal copy will be made.

leaf_size : positive int, default=40
    Number of points at which to switch to brute-force. Changing
    leaf_size will not affect the results of a query, but can
    significantly impact the speed of a query and the memory required
    to store the constructed tree.  The amount of memory needed to
    store the tree scales as approximately n_samples / leaf_size.
    For a specified ``leaf_size``, a leaf node is guaranteed to
    satisfy ``leaf_size <= n_points <= 2 * leaf_size``, except in
    the case that ``n_samples < leaf_size``.

metric : str or DistanceMetric object
    the distance metric to use for the tree.  Default='minkowski'
    with p=2 (that is, a euclidean metric). See the documentation
    of the DistanceMetric class for a list of available metrics.
    {binary_tree}.valid_metrics gives a list of the metrics which
    are valid for {BinaryTree}.

Additional keywords are passed to the distance metric class.
Note: Callable functions in the metric parameter are NOT supported for KDTree
and Ball Tree. Function call overhead will result in very poor performance.

Attributes
----------
data : memory view
    The training data

Examples
--------
Query for k-nearest neighbors

    >>> import numpy as np
    >>> rng = np.random.RandomState(0)
    >>> X = rng.random_sample((10, 3))  # 10 points in 3 dimensions
    >>> tree = {BinaryTree}(X, leaf_size=2)              # doctest: +SKIP
    >>> dist, ind = tree.query(X[:1], k=3)                # doctest: +SKIP
    >>> print(ind)  # indices of 3 closest neighbors
    [0 3 1]
    >>> print(dist)  # distances to 3 closest neighbors
    [ 0.          0.19662693  0.29473397]

Pickle and Unpickle a tree.  Note that the state of the tree is saved in the
pickle operation: the tree needs not be rebuilt upon unpickling.

    >>> import numpy as np
    >>> import pickle
    >>> rng = np.random.RandomState(0)
    >>> X = rng.random_sample((10, 3))  # 10 points in 3 dimensions
    >>> tree = {BinaryTree}(X, leaf_size=2)        # doctest: +SKIP
    >>> s = pickle.dumps(tree)                     # doctest: +SKIP
    >>> tree_copy = pickle.loads(s)                # doctest: +SKIP
    >>> dist, ind = tree_copy.query(X[:1], k=3)     # doctest: +SKIP
    >>> print(ind)  # indices of 3 closest neighbors
    [0 3 1]
    >>> print(dist)  # distances to 3 closest neighbors
    [ 0.          0.19662693  0.29473397]

Query for neighbors within a given radius

    >>> import numpy as np
    >>> rng = np.random.RandomState(0)
    >>> X = rng.random_sample((10, 3))  # 10 points in 3 dimensions
    >>> tree = {BinaryTree}(X, leaf_size=2)     # doctest: +SKIP
    >>> print(tree.query_radius(X[:1], r=0.3, count_only=True))
    3
    >>> ind = tree.query_radius(X[:1], r=0.3)  # doctest: +SKIP
    >>> print(ind)  # indices of neighbors within distance 0.3
    [3 0 1]


Compute a gaussian kernel density estimate:

    >>> import numpy as np
    >>> rng = np.random.RandomState(42)
    >>> X = rng.random_sample((100, 3))
    >>> tree = {BinaryTree}(X)                # doctest: +SKIP
    >>> tree.kernel_density(X[:3], h=0.1, kernel='gaussian')
    array([ 6.94114649,  7.83281226,  7.2071716 ])

Compute a two-point auto-correlation function

    >>> import numpy as np
    >>> rng = np.random.RandomState(0)
    >>> X = rng.random_sample((30, 3))
    >>> r = np.linspace(0, 1, 5)
    >>> tree = {BinaryTree}(X)                # doctest: +SKIP
    >>> tree.two_point_correlation(X, r)
    array([ 30,  62, 278, 580, 820])

"""


######################################################################
# Utility functions
cdef DTYPE_t logaddexp(DTYPE_t x1, DTYPE_t x2):
    """logaddexp(x1, x2) -> log(exp(x1) + exp(x2))"""
    cdef DTYPE_t a = fmax(x1, x2)
    if a == NEG_INF:
        return NEG_INF
    else:
        return a + log(exp(x1 - a) + exp(x2 - a))

cdef DTYPE_t logsubexp(DTYPE_t x1, DTYPE_t x2):
    """logsubexp(x1, x2) -> log(exp(x1) - exp(x2))"""
    if x1 <= x2:
        return NEG_INF
    else:
        return x1 + log(1 - exp(x2 - x1))


######################################################################
# Kernel functions
#
# Note: Kernels assume dist is non-negative and h is positive
#       All kernel functions are normalized such that K(0, h) = 1.
#       The fully normalized kernel is:
#         K = exp[kernel_norm(h, d, kernel) + compute_kernel(dist, h, kernel)]
#       The code only works with non-negative kernels: i.e. K(d, h) >= 0
#       for all valid d and h.  Note that for precision, the log of both
#       the kernel and kernel norm is returned.
cdef enum KernelType:
    GAUSSIAN_KERNEL = 1
    TOPHAT_KERNEL = 2
    EPANECHNIKOV_KERNEL = 3
    EXPONENTIAL_KERNEL = 4
    LINEAR_KERNEL = 5
    COSINE_KERNEL = 6


cdef inline DTYPE_t log_gaussian_kernel(DTYPE_t dist, DTYPE_t h):
    """log of the gaussian kernel for bandwidth h (unnormalized)"""
    return -0.5 * (dist * dist) / (h * h)


cdef inline DTYPE_t log_tophat_kernel(DTYPE_t dist, DTYPE_t h):
    """log of the tophat kernel for bandwidth h (unnormalized)"""
    if dist < h:
        return 0.0
    else:
        return NEG_INF


cdef inline DTYPE_t log_epanechnikov_kernel(DTYPE_t dist, DTYPE_t h):
    """log of the epanechnikov kernel for bandwidth h (unnormalized)"""
    if dist < h:
        return log(1.0 - (dist * dist) / (h * h))
    else:
        return NEG_INF


cdef inline DTYPE_t log_exponential_kernel(DTYPE_t dist, DTYPE_t h):
    """log of the exponential kernel for bandwidth h (unnormalized)"""
    return -dist / h


cdef inline DTYPE_t log_linear_kernel(DTYPE_t dist, DTYPE_t h):
    """log of the linear kernel for bandwidth h (unnormalized)"""
    if dist < h:
        return log(1 - dist / h)
    else:
        return NEG_INF


cdef inline DTYPE_t log_cosine_kernel(DTYPE_t dist, DTYPE_t h):
    """log of the cosine kernel for bandwidth h (unnormalized)"""
    if dist < h:
        return log(cos(0.5 * PI * dist / h))
    else:
        return NEG_INF


cdef inline DTYPE_t compute_log_kernel(DTYPE_t dist, DTYPE_t h,
                                       KernelType kernel):
    """Given a KernelType enumeration, compute the appropriate log-kernel"""
    if kernel == GAUSSIAN_KERNEL:
        return log_gaussian_kernel(dist, h)
    elif kernel == TOPHAT_KERNEL:
        return log_tophat_kernel(dist, h)
    elif kernel == EPANECHNIKOV_KERNEL:
        return log_epanechnikov_kernel(dist, h)
    elif kernel == EXPONENTIAL_KERNEL:
        return log_exponential_kernel(dist, h)
    elif kernel == LINEAR_KERNEL:
        return log_linear_kernel(dist, h)
    elif kernel == COSINE_KERNEL:
        return log_cosine_kernel(dist, h)


#------------------------------------------------------------
# Kernel norms are defined via the volume element V_n
# and surface element S_(n-1) of an n-sphere.
cdef DTYPE_t logVn(ITYPE_t n):
    """V_n = pi^(n/2) / gamma(n/2 - 1)"""
    return 0.5 * n * LOG_PI - lgamma(0.5 * n + 1)


cdef DTYPE_t logSn(ITYPE_t n):
    """V_(n+1) = int_0^1 S_n r^n dr"""
    return LOG_2PI + logVn(n - 1)


cdef DTYPE_t _log_kernel_norm(DTYPE_t h, ITYPE_t d,
                              KernelType kernel) except -1:
    """Given a KernelType enumeration, compute the kernel normalization.

    h is the bandwidth, d is the dimension.
    """
    cdef DTYPE_t tmp, factor = 0
    cdef ITYPE_t k
    if kernel == GAUSSIAN_KERNEL:
        factor = 0.5 * d * LOG_2PI
    elif kernel == TOPHAT_KERNEL:
        factor = logVn(d)
    elif kernel == EPANECHNIKOV_KERNEL:
        factor = logVn(d) + log(2. / (d + 2.))
    elif kernel == EXPONENTIAL_KERNEL:
        factor = logSn(d - 1) + lgamma(d)
    elif kernel == LINEAR_KERNEL:
        factor = logVn(d) - log(d + 1.)
    elif kernel == COSINE_KERNEL:
        # this is derived from a chain rule integration
        factor = 0
        tmp = 2. / PI
        for k in range(1, d + 1, 2):
            factor += tmp
            tmp *= -(d - k) * (d - k - 1) * (2. / PI) ** 2
        factor = log(factor) + logSn(d - 1)
    else:
        raise ValueError("Kernel code not recognized")
    return -factor - d * log(h)


def kernel_norm(h, d, kernel, return_log=False):
    """Given a string specification of a kernel, compute the normalization.

    Parameters
    ----------
    h : float
        The bandwidth of the kernel.
    d : int
        The dimension of the space in which the kernel norm is computed.
    kernel : str
        The kernel identifier.  Must be one of
        ['gaussian'|'tophat'|'epanechnikov'|
         'exponential'|'linear'|'cosine']
    return_log : bool, default=False
        If True, return the log of the kernel norm.  Otherwise, return the
        kernel norm.
    Returns
    -------
    knorm or log_knorm : float
        the kernel norm or logarithm of the kernel norm.
    """
    if kernel == 'gaussian':
        result = _log_kernel_norm(h, d, GAUSSIAN_KERNEL)
    elif kernel == 'tophat':
        result = _log_kernel_norm(h, d, TOPHAT_KERNEL)
    elif kernel == 'epanechnikov':
        result = _log_kernel_norm(h, d, EPANECHNIKOV_KERNEL)
    elif kernel == 'exponential':
        result = _log_kernel_norm(h, d, EXPONENTIAL_KERNEL)
    elif kernel == 'linear':
        result = _log_kernel_norm(h, d, LINEAR_KERNEL)
    elif kernel == 'cosine':
        result = _log_kernel_norm(h, d, COSINE_KERNEL)
    else:
        raise ValueError('kernel not recognized')

    if return_log:
        return result
    else:
        return np.exp(result)


######################################################################
# Tree Utility Routines
cdef inline void swap(DITYPE_t* arr, ITYPE_t i1, ITYPE_t i2):
    """swap the values at index i1 and i2 of arr"""
    cdef DITYPE_t tmp = arr[i1]
    arr[i1] = arr[i2]
    arr[i2] = tmp


cdef inline void dual_swap(DTYPE_t* darr, ITYPE_t* iarr,
                           ITYPE_t i1, ITYPE_t i2) nogil:
    """swap the values at inex i1 and i2 of both darr and iarr"""
    cdef DTYPE_t dtmp = darr[i1]
    darr[i1] = darr[i2]
    darr[i2] = dtmp

    cdef ITYPE_t itmp = iarr[i1]
    iarr[i1] = iarr[i2]
    iarr[i2] = itmp


cdef class NeighborsHeap:
    """A max-heap structure to keep track of distances/indices of neighbors

    This implements an efficient pre-allocated set of fixed-size heaps
    for chasing neighbors, holding both an index and a distance.
    When any row of the heap is full, adding an additional point will push
    the furthest point off the heap.

    Parameters
    ----------
    n_pts : int
        the number of heaps to use
    n_nbrs : int
        the size of each heap.
    """
    cdef np.ndarray distances_arr
    cdef np.ndarray indices_arr

    cdef DTYPE_t[:, ::1] distances
    cdef ITYPE_t[:, ::1] indices

    def __cinit__(self):
        self.distances_arr = np.zeros((1, 1), dtype=DTYPE, order='C')
        self.indices_arr = np.zeros((1, 1), dtype=ITYPE, order='C')
        self.distances = get_memview_DTYPE_2D(self.distances_arr)
        self.indices = get_memview_ITYPE_2D(self.indices_arr)

    def __init__(self, n_pts, n_nbrs):
        self.distances_arr = np.full((n_pts, n_nbrs), np.inf, dtype=DTYPE,
                                     order='C')
        self.indices_arr = np.zeros((n_pts, n_nbrs), dtype=ITYPE, order='C')
        self.distances = get_memview_DTYPE_2D(self.distances_arr)
        self.indices = get_memview_ITYPE_2D(self.indices_arr)

    def get_arrays(self, sort=True):
        """Get the arrays of distances and indices within the heap.

        If sort=True, then simultaneously sort the indices and distances,
        so the closer points are listed first.
        """
        if sort:
            self._sort()
        return self.distances_arr, self.indices_arr

    cdef inline DTYPE_t largest(self, ITYPE_t row) nogil except -1:
        """Return the largest distance in the given row"""
        return self.distances[row, 0]

    def push(self, ITYPE_t row, DTYPE_t val, ITYPE_t i_val):
        return self._push(row, val, i_val)

    cdef int _push(self, ITYPE_t row, DTYPE_t val,
                   ITYPE_t i_val) nogil except -1:
        """push (val, i_val) into the given row"""
        cdef ITYPE_t i, ic1, ic2, i_swap
        cdef ITYPE_t size = self.distances.shape[1]
        cdef DTYPE_t* dist_arr = &self.distances[row, 0]
        cdef ITYPE_t* ind_arr = &self.indices[row, 0]

        # check if val should be in heap
        if val > dist_arr[0]:
            return 0

        # insert val at position zero
        dist_arr[0] = val
        ind_arr[0] = i_val

        # descend the heap, swapping values until the max heap criterion is met
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

        return 0

    cdef int _sort(self) except -1:
        """simultaneously sort the distances and indices"""
        cdef DTYPE_t[:, ::1] distances = self.distances
        cdef ITYPE_t[:, ::1] indices = self.indices
        cdef ITYPE_t row
        for row in range(distances.shape[0]):
            _simultaneous_sort(&distances[row, 0],
                               &indices[row, 0],
                               distances.shape[1])
        return 0


cdef int _simultaneous_sort(DTYPE_t* dist, ITYPE_t* idx,
                            ITYPE_t size) nogil except -1:
    """
    Perform a recursive quicksort on the dist array, simultaneously
    performing the same swaps on the idx array.  The equivalent in
    numpy (though quite a bit slower) is

    def simultaneous_sort(dist, idx):
        i = np.argsort(dist)
        return dist[i], idx[i]
    """
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
    return 0

#------------------------------------------------------------
# find_node_split_dim:
#  this computes the equivalent of
#  j_max = np.argmax(np.max(data, 0) - np.min(data, 0))
cdef ITYPE_t find_node_split_dim(DTYPE_t* data,
                                 ITYPE_t* node_indices,
                                 ITYPE_t n_features,
                                 ITYPE_t n_points) except -1:
    """Find the dimension with the largest spread.

    Parameters
    ----------
    data : double pointer
        Pointer to a 2D array of the training data, of shape [N, n_features].
        N must be greater than any of the values in node_indices.
    node_indices : int pointer
        Pointer to a 1D array of length n_points.  This lists the indices of
        each of the points within the current node.

    Returns
    -------
    i_max : int
        The index of the feature (dimension) within the node that has the
        largest spread.

    Notes
    -----
    In numpy, this operation is equivalent to

    def find_node_split_dim(data, node_indices):
        return np.argmax(data[node_indices].max(0) - data[node_indices].min(0))

    The cython version is much more efficient in both computation and memory.
    """
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


cdef int partition_node_indices(DTYPE_t* data,
                                ITYPE_t* node_indices,
                                ITYPE_t split_dim,
                                ITYPE_t split_index,
                                ITYPE_t n_features,
                                ITYPE_t n_points) except -1:
    """Partition points in the node into two equal-sized groups.

    Upon return, the values in node_indices will be rearranged such that
    (assuming numpy-style indexing):

        data[node_indices[0:split_index], split_dim]
          <= data[node_indices[split_index], split_dim]

    and

        data[node_indices[split_index], split_dim]
          <= data[node_indices[split_index:n_points], split_dim]

    The algorithm is essentially a partial in-place quicksort around a
    set pivot.

    Parameters
    ----------
    data : double pointer
        Pointer to a 2D array of the training data, of shape [N, n_features].
        N must be greater than any of the values in node_indices.
    node_indices : int pointer
        Pointer to a 1D array of length n_points.  This lists the indices of
        each of the points within the current node.  This will be modified
        in-place.
    split_dim : int
        the dimension on which to split.  This will usually be computed via
        the routine ``find_node_split_dim``
    split_index : int
        the index within node_indices around which to split the points.

    Returns
    -------
    status : int
        integer exit status.  On return, the contents of node_indices are
        modified as noted above.
    """
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

    return 0


######################################################################
# NodeHeap : min-heap used to keep track of nodes during
#            breadth-first query
cdef inline void swap_nodes(NodeHeapData_t* arr, ITYPE_t i1, ITYPE_t i2):
    cdef NodeHeapData_t tmp = arr[i1]
    arr[i1] = arr[i2]
    arr[i2] = tmp


cdef class NodeHeap:
    """NodeHeap

    This is a min-heap implementation for keeping track of nodes
    during a breadth-first search.  Unlike the NeighborsHeap above,
    the NodeHeap does not have a fixed size and must be able to grow
    as elements are added.

    Internally, the data is stored in a simple binary heap which meets
    the min heap condition:

        heap[i].val < min(heap[2 * i + 1].val, heap[2 * i + 2].val)
    """
    cdef np.ndarray data_arr
    cdef NodeHeapData_t[::1] data
    cdef ITYPE_t n

    def __cinit__(self):
        self.data_arr = np.zeros(1, dtype=NodeHeapData, order='C')
        self.data = get_memview_NodeHeapData_1D(self.data_arr)

    def __init__(self, size_guess=100):
        size_guess = max(size_guess, 1)  # need space for at least one item
        self.data_arr = np.zeros(size_guess, dtype=NodeHeapData, order='C')
        self.data = get_memview_NodeHeapData_1D(self.data_arr)
        self.n = size_guess
        self.clear()

    cdef int resize(self, ITYPE_t new_size) except -1:
        """Resize the heap to be either larger or smaller"""
        cdef NodeHeapData_t *data_ptr
        cdef NodeHeapData_t *new_data_ptr
        cdef ITYPE_t i
        cdef ITYPE_t size = self.data.shape[0]
        cdef np.ndarray new_data_arr = np.zeros(new_size,
                                                dtype=NodeHeapData)
        cdef NodeHeapData_t[::1] new_data =\
                                    get_memview_NodeHeapData_1D(new_data_arr)

        if size > 0 and new_size > 0:
            data_ptr = &self.data[0]
            new_data_ptr = &new_data[0]
            for i in range(min(size, new_size)):
                new_data_ptr[i] = data_ptr[i]

        if new_size < size:
            self.n = new_size

        self.data = new_data
        self.data_arr = new_data_arr
        return 0

    cdef int push(self, NodeHeapData_t data) except -1:
        """Push a new item onto the heap"""
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
        return 0

    cdef NodeHeapData_t peek(self):
        """Peek at the root of the heap, without removing it"""
        return self.data[0]

    cdef NodeHeapData_t pop(self):
        """Remove the root of the heap, and update the remaining nodes"""
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

    cdef void clear(self):
        """Clear the heap"""
        self.n = 0


######################################################################
# newObj function
#  this is a helper function for pickling
def newObj(obj):
    return obj.__new__(obj)


######################################################################
# define the reverse mapping of VALID_METRICS
from ._dist_metrics import get_valid_metric_ids
VALID_METRIC_IDS = get_valid_metric_ids(VALID_METRICS)


######################################################################
# Binary Tree class
cdef class BinaryTree:

    cdef np.ndarray data_arr
    cdef np.ndarray sample_weight_arr
    cdef np.ndarray idx_array_arr
    cdef np.ndarray node_data_arr
    cdef np.ndarray node_bounds_arr

    cdef readonly DTYPE_t[:, ::1] data
    cdef readonly DTYPE_t[::1] sample_weight
    cdef public DTYPE_t sum_weight
    cdef public ITYPE_t[::1] idx_array
    cdef public NodeData_t[::1] node_data
    cdef public DTYPE_t[:, :, ::1] node_bounds

    cdef ITYPE_t leaf_size
    cdef ITYPE_t n_levels
    cdef ITYPE_t n_nodes

    cdef DistanceMetric dist_metric
    cdef int euclidean

    # variables to keep track of building & querying stats
    cdef int n_trims
    cdef int n_leaves
    cdef int n_splits
    cdef int n_calls

    valid_metrics = VALID_METRIC_IDS

    # Use cinit to initialize all arrays to empty: this will prevent memory
    # errors and seg-faults in rare cases where __init__ is not called
    def __cinit__(self):
        self.data_arr = np.empty((1, 1), dtype=DTYPE, order='C')
        self.sample_weight_arr = np.empty(1, dtype=DTYPE, order='C')
        self.idx_array_arr = np.empty(1, dtype=ITYPE, order='C')
        self.node_data_arr = np.empty(1, dtype=NodeData, order='C')
        self.node_bounds_arr = np.empty((1, 1, 1), dtype=DTYPE)

        self.data = get_memview_DTYPE_2D(self.data_arr)
        self.sample_weight = get_memview_DTYPE_1D(self.sample_weight_arr)
        self.idx_array = get_memview_ITYPE_1D(self.idx_array_arr)
        self.node_data = get_memview_NodeData_1D(self.node_data_arr)
        self.node_bounds = get_memview_DTYPE_3D(self.node_bounds_arr)

        self.leaf_size = 0
        self.n_levels = 0
        self.n_nodes = 0

        self.euclidean = False

        self.n_trims = 0
        self.n_leaves = 0
        self.n_splits = 0
        self.n_calls = 0

    def __init__(self, data,
                 leaf_size=40, metric='minkowski', sample_weight=None, **kwargs):
        # validate data
        if data.size == 0:
            raise ValueError("X is an empty array")

        if leaf_size < 1:
            raise ValueError("leaf_size must be greater than or equal to 1")

        n_samples = data.shape[0]
        n_features = data.shape[1]

        self.data_arr = np.asarray(data, dtype=DTYPE, order='C')
        self.leaf_size = leaf_size
        self.dist_metric = DistanceMetric.get_metric(metric, **kwargs)
        self.euclidean = (self.dist_metric.__class__.__name__
                          == 'EuclideanDistance')

        metric = self.dist_metric.__class__.__name__
        if metric not in VALID_METRICS:
            raise ValueError('metric {metric} is not valid for '
                             '{BinaryTree}'.format(metric=metric,
                                                   **DOC_DICT))
        self.dist_metric._validate_data(data)

        # determine number of levels in the tree, and from this
        # the number of nodes in the tree.  This results in leaf nodes
        # with numbers of points between leaf_size and 2 * leaf_size
        self.n_levels = np.log2(fmax(1, (n_samples - 1) / self.leaf_size)) + 1
        self.n_nodes = (2 ** self.n_levels) - 1

        # allocate arrays for storage
        self.idx_array_arr = np.arange(n_samples, dtype=ITYPE)
        self.node_data_arr = np.zeros(self.n_nodes, dtype=NodeData)

        self._update_sample_weight(n_samples, sample_weight)
        self._update_memviews()

        # Allocate tree-specific data
        allocate_data(self, self.n_nodes, n_features)
        self._recursive_build(0, 0, n_samples)

    def _update_sample_weight(self, n_samples, sample_weight):
        if sample_weight is not None:
            self.sample_weight_arr = np.asarray(
                sample_weight, dtype=DTYPE, order='C')
            self.sample_weight = get_memview_DTYPE_1D(
                self.sample_weight_arr)
            self.sum_weight = np.sum(self.sample_weight)
        else:
            self.sample_weight = None
            self.sample_weight_arr = np.empty(1, dtype=DTYPE, order='C')
            self.sum_weight = <DTYPE_t> n_samples

    def _update_memviews(self):
        self.data = get_memview_DTYPE_2D(self.data_arr)
        self.idx_array = get_memview_ITYPE_1D(self.idx_array_arr)
        self.node_data = get_memview_NodeData_1D(self.node_data_arr)
        self.node_bounds = get_memview_DTYPE_3D(self.node_bounds_arr)


    def __reduce__(self):
        """
        reduce method used for pickling
        """
        return (newObj, (type(self),), self.__getstate__())

    def __getstate__(self):
        """
        get state for pickling
        """
        if self.sample_weight is not None:
            # pass the numpy array
            sample_weight_arr = self.sample_weight_arr
        else:
            # pass None to avoid confusion with the empty place holder
            # of size 1 from __cinit__
            sample_weight_arr = None
        return (self.data_arr,
                self.idx_array_arr,
                self.node_data_arr,
                self.node_bounds_arr,
                int(self.leaf_size),
                int(self.n_levels),
                int(self.n_nodes),
                int(self.n_trims),
                int(self.n_leaves),
                int(self.n_splits),
                int(self.n_calls),
                self.dist_metric,
                sample_weight_arr)

    def __setstate__(self, state):
        """
        set state for pickling
        """
        self.data_arr = state[0]
        self.idx_array_arr = state[1]
        self.node_data_arr = state[2]
        self.node_bounds_arr = state[3]
        self.leaf_size = state[4]
        self.n_levels = state[5]
        self.n_nodes = state[6]
        self.n_trims = state[7]
        self.n_leaves = state[8]
        self.n_splits = state[9]
        self.n_calls = state[10]
        self.dist_metric = state[11]
        sample_weight_arr = state[12]

        self.euclidean = (self.dist_metric.__class__.__name__
                          == 'EuclideanDistance')
        n_samples = self.data_arr.shape[0]
        self._update_sample_weight(n_samples, sample_weight_arr)
        self._update_memviews()

    def get_tree_stats(self):
        """
        get_tree_stats(self)

        Get tree status.

        Returns
        -------
        tree_stats: tuple of int
            (number of trims, number of leaves, number of splits)
        """
        return (self.n_trims, self.n_leaves, self.n_splits)

    def reset_n_calls(self):
        """
        reset_n_calls(self)

        Reset number of calls to 0.
        """
        self.n_calls = 0

    def get_n_calls(self):
        """
        get_n_calls(self)

        Get number of calls.

        Returns
        -------
        n_calls: int
            number of distance computation calls
        """
        return self.n_calls

    def get_arrays(self):
        """
        get_arrays(self)

        Get data and node arrays.

        Returns
        -------
        arrays: tuple of array
            Arrays for storing tree data, index, node data and node bounds.
        """
        return (self.data_arr, self.idx_array_arr,
                self.node_data_arr, self.node_bounds_arr)

    cdef inline DTYPE_t dist(self, DTYPE_t* x1, DTYPE_t* x2,
                             ITYPE_t size) nogil except -1:
        """Compute the distance between arrays x1 and x2"""
        self.n_calls += 1
        if self.euclidean:
            return euclidean_dist(x1, x2, size)
        else:
            return self.dist_metric.dist(x1, x2, size)

    cdef inline DTYPE_t rdist(self, DTYPE_t* x1, DTYPE_t* x2,
                              ITYPE_t size) nogil except -1:
        """Compute the reduced distance between arrays x1 and x2.

        The reduced distance, defined for some metrics, is a quantity which
        is more efficient to compute than the distance, but preserves the
        relative rankings of the true distance.  For example, the reduced
        distance for the Euclidean metric is the squared-euclidean distance.
        """
        self.n_calls += 1
        if self.euclidean:
            return euclidean_rdist(x1, x2, size)
        else:
            return self.dist_metric.rdist(x1, x2, size)

    cdef int _recursive_build(self, ITYPE_t i_node, ITYPE_t idx_start,
                              ITYPE_t idx_end) except -1:
        """Recursively build the tree.

        Parameters
        ----------
        i_node : int
            the node for the current step
        idx_start, idx_end : int
            the bounding indices in the idx_array which define the points that
            belong to this node.
        """
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

        query the tree for the k nearest neighbors

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            An array of points to query
        k : int, default=1
            The number of nearest neighbors to return
        return_distance : bool, default=True
            if True, return a tuple (d, i) of distances and indices
            if False, return array i
        dualtree : bool, default=False
            if True, use the dual tree formalism for the query: a tree is
            built for the query points, and the pair of trees is used to
            efficiently search this space.  This can lead to better
            performance as the number of points grows large.
        breadth_first : bool, default=False
            if True, then query the nodes in a breadth-first manner.
            Otherwise, query the nodes in a depth-first manner.
        sort_results : bool, default=True
            if True, then distances and indices of each point are sorted
            on return, so that the first column contains the closest points.
            Otherwise, neighbors are returned in an arbitrary order.

        Returns
        -------
        i    : if return_distance == False
        (d,i) : if return_distance == True

        d : ndarray of shape X.shape[:-1] + k, dtype=double
            Each entry gives the list of distances to the neighbors of the
            corresponding point.

        i : ndarray of shape X.shape[:-1] + k, dtype=int
            Each entry gives the list of indices of neighbors of the
            corresponding point.
        """
        # XXX: we should allow X to be a pre-built tree.
        X = check_array(X, dtype=DTYPE, order='C')

        if X.shape[X.ndim - 1] != self.data.shape[1]:
            raise ValueError("query data dimension must "
                             "match training data dimension")

        if self.data.shape[0] < k:
            raise ValueError("k must be less than or equal "
                             "to the number of training points")

        # flatten X, and save original shape information
        np_Xarr = X.reshape((-1, self.data.shape[1]))
        cdef DTYPE_t[:, ::1] Xarr = get_memview_DTYPE_2D(np_Xarr)
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
            other = self.__class__(np_Xarr, metric=self.dist_metric,
                                   leaf_size=self.leaf_size)
            if breadth_first:
                self._query_dual_breadthfirst(other, heap, nodeheap)
            else:
                reduced_dist_LB = min_rdist_dual(self, 0, other, 0)
                bounds = np.full(other.node_data.shape[0], np.inf)
                self._query_dual_depthfirst(0, other, 0, bounds,
                                            heap, reduced_dist_LB)

        else:
            pt = &Xarr[0, 0]
            if breadth_first:
                for i in range(Xarr.shape[0]):
                    self._query_single_breadthfirst(pt, i, heap, nodeheap)
                    pt += Xarr.shape[1]
            else:
                with nogil:
                    for i in range(Xarr.shape[0]):
                        reduced_dist_LB = min_rdist(self, 0, pt)
                        self._query_single_depthfirst(0, pt, i, heap,
                                                      reduced_dist_LB)
                        pt += Xarr.shape[1]

        distances, indices = heap.get_arrays(sort=sort_results)
        distances = self.dist_metric.rdist_to_dist(distances)

        # deflatten results
        if return_distance:
            return (distances.reshape(X.shape[:X.ndim - 1] + (k,)),
                    indices.reshape(X.shape[:X.ndim - 1] + (k,)))
        else:
            return indices.reshape(X.shape[:X.ndim - 1] + (k,))

    def query_radius(self, X, r, int return_distance=False,
                     int count_only=False, int sort_results=False):
        """
        query_radius(X, r, return_distance=False,
        count_only=False, sort_results=False)

        query the tree for neighbors within a radius r

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            An array of points to query
        r : distance within which neighbors are returned
            r can be a single value, or an array of values of shape
            x.shape[:-1] if different radii are desired for each point.
        return_distance : bool, default=False
            if True,  return distances to neighbors of each point
            if False, return only neighbors
            Note that unlike the query() method, setting return_distance=True
            here adds to the computation time.  Not all distances need to be
            calculated explicitly for return_distance=False.  Results are
            not sorted by default: see ``sort_results`` keyword.
        count_only : bool, default=False
            if True,  return only the count of points within distance r
            if False, return the indices of all points within distance r
            If return_distance==True, setting count_only=True will
            result in an error.
        sort_results : bool, default=False
            if True, the distances and indices will be sorted before being
            returned.  If False, the results will not be sorted.  If
            return_distance == False, setting sort_results = True will
            result in an error.

        Returns
        -------
        count       : if count_only == True
        ind         : if count_only == False and return_distance == False
        (ind, dist) : if count_only == False and return_distance == True

        count : ndarray of shape X.shape[:-1], dtype=int
            Each entry gives the number of neighbors within a distance r of the
            corresponding point.

        ind : ndarray of shape X.shape[:-1], dtype=object
            Each element is a numpy integer array listing the indices of
            neighbors of the corresponding point.  Note that unlike
            the results of a k-neighbors query, the returned neighbors
            are not sorted by distance by default.

        dist : ndarray of shape X.shape[:-1], dtype=object
            Each element is a numpy double array listing the distances
            corresponding to indices in i.
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
        cdef ITYPE_t[::1] idx_arr_i, counts
        cdef DTYPE_t* pt
        cdef ITYPE_t** indices = NULL
        cdef DTYPE_t** distances = NULL

        # validate X and prepare for query
        X = check_array(X, dtype=DTYPE, order='C')

        if X.shape[X.ndim - 1] != self.data.shape[1]:
            raise ValueError("query data dimension must "
                             "match training data dimension")

        cdef DTYPE_t[:, ::1] Xarr =\
                get_memview_DTYPE_2D(X.reshape((-1, self.data.shape[1])))

        # prepare r for query
        r = np.asarray(r, dtype=DTYPE, order='C')
        r = np.atleast_1d(r)
        if r.shape == (1,):
            r = np.full(X.shape[:X.ndim - 1], r[0], dtype=DTYPE)
        else:
            if r.shape != X.shape[:X.ndim - 1]:
                raise ValueError("r must be broadcastable to X.shape")

        rarr_np = r.reshape(-1)  # store explicitly to keep in scope
        cdef DTYPE_t[::1] rarr = get_memview_DTYPE_1D(rarr_np)

        if not count_only:
            indices = <ITYPE_t**>calloc(Xarr.shape[0], sizeof(ITYPE_t*))
            if indices == NULL:
                raise MemoryError()
            if return_distance:
                distances = <DTYPE_t**>calloc(Xarr.shape[0], sizeof(DTYPE_t*))
                if distances == NULL:
                    free(indices)
                    raise MemoryError()

        np_idx_arr = np.zeros(self.data.shape[0], dtype=ITYPE)
        idx_arr_i = get_memview_ITYPE_1D(np_idx_arr)

        np_dist_arr = np.zeros(self.data.shape[0], dtype=DTYPE)
        dist_arr_i = get_memview_DTYPE_1D(np_dist_arr)

        counts_arr = np.zeros(Xarr.shape[0], dtype=ITYPE)
        counts = get_memview_ITYPE_1D(counts_arr)

        pt = &Xarr[0, 0]
        memory_error = False
        with nogil:
            for i in range(Xarr.shape[0]):
                counts[i] = self._query_radius_single(0, pt, rarr[i],
                                                      &idx_arr_i[0],
                                                      &dist_arr_i[0],
                                                      0, count_only,
                                                      return_distance)
                pt += n_features

                if count_only:
                    continue

                if sort_results:
                    _simultaneous_sort(&dist_arr_i[0], &idx_arr_i[0],
                                       counts[i])

                # equivalent to: indices[i] = np_idx_arr[:counts[i]].copy()
                indices[i] = <ITYPE_t*>malloc(counts[i] * sizeof(ITYPE_t))
                if indices[i] == NULL:
                    memory_error = True
                    break
                memcpy(indices[i], &idx_arr_i[0], counts[i] * sizeof(ITYPE_t))

                if return_distance:
                    # equivalent to: distances[i] = np_dist_arr[:counts[i]].copy()
                    distances[i] = <DTYPE_t*>malloc(counts[i] * sizeof(DTYPE_t))
                    if distances[i] == NULL:
                        memory_error = True
                        break
                    memcpy(distances[i], &dist_arr_i[0], counts[i] * sizeof(DTYPE_t))

        try:
            if memory_error:
                raise MemoryError()

            if count_only:
                # deflatten results
                return counts_arr.reshape(X.shape[:X.ndim - 1])
            elif return_distance:
                indices_npy = np.zeros(Xarr.shape[0], dtype='object')
                distances_npy = np.zeros(Xarr.shape[0], dtype='object')
                for i in range(Xarr.shape[0]):
                    # make a new numpy array that wraps the existing data
                    indices_npy[i] = np.PyArray_SimpleNewFromData(1, &counts[i], np.NPY_INTP, indices[i])
                    # make sure the data will be freed when the numpy array is garbage collected
                    PyArray_ENABLEFLAGS(indices_npy[i], np.NPY_OWNDATA)
                    # make sure the data is not freed twice
                    indices[i] = NULL

                    # make a new numpy array that wraps the existing data
                    distances_npy[i] = np.PyArray_SimpleNewFromData(1, &counts[i], np.NPY_DOUBLE, distances[i])
                    # make sure the data will be freed when the numpy array is garbage collected
                    PyArray_ENABLEFLAGS(distances_npy[i], np.NPY_OWNDATA)
                    # make sure the data is not freed twice
                    distances[i] = NULL

                # deflatten results
                return (indices_npy.reshape(X.shape[:X.ndim - 1]),
                        distances_npy.reshape(X.shape[:X.ndim - 1]))
            else:
                indices_npy = np.zeros(Xarr.shape[0], dtype='object')
                for i in range(Xarr.shape[0]):
                    # make a new numpy array that wraps the existing data
                    indices_npy[i] = np.PyArray_SimpleNewFromData(1, &counts[i], np.NPY_INTP, indices[i])
                    # make sure the data will be freed when the numpy array is garbage collected
                    PyArray_ENABLEFLAGS(indices_npy[i], np.NPY_OWNDATA)
                    # make sure the data is not freed twice
                    indices[i] = NULL

                # deflatten results
                return indices_npy.reshape(X.shape[:X.ndim - 1])
        except:
            # free any buffer that is not owned by a numpy array
            for i in range(Xarr.shape[0]):
                free(indices[i])
                if return_distance:
                    free(distances[i])
            raise
        finally:
            free(indices)
            free(distances)


    def kernel_density(self, X, h, kernel='gaussian',
                       atol=0, rtol=1E-8,
                       breadth_first=True, return_log=False):
        """
        kernel_density(self, X, h, kernel='gaussian', atol=0, rtol=1E-8,
                       breadth_first=True, return_log=False)

        Compute the kernel density estimate at points X with the given kernel,
        using the distance metric specified at tree creation.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            An array of points to query.  Last dimension should match dimension
            of training data.
        h : float
            the bandwidth of the kernel
        kernel : str, default="gaussian"
            specify the kernel to use.  Options are
            - 'gaussian'
            - 'tophat'
            - 'epanechnikov'
            - 'exponential'
            - 'linear'
            - 'cosine'
            Default is kernel = 'gaussian'
        atol, rtol : float, default=0, 1e-8
            Specify the desired relative and absolute tolerance of the result.
            If the true result is K_true, then the returned result K_ret
            satisfies ``abs(K_true - K_ret) < atol + rtol * K_ret``
            The default is zero (i.e. machine precision) for both.
        breadth_first : bool, default=False
            If True, use a breadth-first search.  If False (default) use a
            depth-first search.  Breadth-first is generally faster for
            compact kernels and/or high tolerances.
        return_log : bool, default=False
            Return the logarithm of the result.  This can be more accurate
            than returning the result itself for narrow kernels.

        Returns
        -------
        density : ndarray of shape X.shape[:-1]
            The array of (log)-density evaluations
        """
        cdef DTYPE_t h_c = h
        cdef DTYPE_t log_atol = log(atol)
        cdef DTYPE_t log_rtol = log(rtol)
        cdef DTYPE_t log_min_bound, log_max_bound, log_bound_spread
        cdef DTYPE_t dist_LB = 0, dist_UB = 0

        cdef ITYPE_t n_samples = self.data.shape[0]
        cdef ITYPE_t n_features = self.data.shape[1]
        cdef ITYPE_t i
        cdef KernelType kernel_c

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

        cdef DTYPE_t log_knorm = _log_kernel_norm(h_c, n_features, kernel_c)

        # validate X and prepare for query
        X = check_array(X, dtype=DTYPE, order='C')

        if X.shape[X.ndim - 1] != n_features:
            raise ValueError("query data dimension must "
                             "match training data dimension")
        Xarr_np = X.reshape((-1, n_features))
        cdef DTYPE_t[:, ::1] Xarr = get_memview_DTYPE_2D(Xarr_np)

        log_density_arr = np.zeros(Xarr.shape[0], dtype=DTYPE)
        cdef DTYPE_t[::1] log_density = get_memview_DTYPE_1D(log_density_arr)

        cdef DTYPE_t* pt = &Xarr[0, 0]

        cdef NodeHeap nodeheap
        if breadth_first:
            nodeheap = NodeHeap(self.data.shape[0] // self.leaf_size)
        cdef DTYPE_t[::1] node_log_min_bounds
        cdef DTYPE_t[::1] node_bound_widths
        # TODO: implement dual tree approach.
        #       this is difficult because of the need to cache values
        #       computed between node pairs.
        if breadth_first:
            node_log_min_bounds_arr = np.full(self.n_nodes, -np.inf)
            node_log_min_bounds = get_memview_DTYPE_1D(node_log_min_bounds_arr)
            node_bound_widths_arr = np.zeros(self.n_nodes)
            node_bound_widths = get_memview_DTYPE_1D(node_bound_widths_arr)
            for i in range(Xarr.shape[0]):
                log_density[i] = self._kde_single_breadthfirst(
                                            pt, kernel_c, h_c,
                                            log_knorm, log_atol, log_rtol,
                                            nodeheap,
                                            &node_log_min_bounds[0],
                                            &node_bound_widths[0])
                pt += n_features
        else:
            for i in range(Xarr.shape[0]):
                min_max_dist(self, 0, pt, &dist_LB, &dist_UB)
                # compute max & min bounds on density within top node
                log_min_bound = (log(self.sum_weight) +
                                 compute_log_kernel(dist_UB,
                                                    h_c, kernel_c))
                log_max_bound = (log(self.sum_weight) +
                                 compute_log_kernel(dist_LB,
                                                    h_c, kernel_c))
                log_bound_spread = logsubexp(log_max_bound, log_min_bound)
                self._kde_single_depthfirst(0, pt, kernel_c, h_c,
                                            log_knorm, log_atol, log_rtol,
                                            log_min_bound,
                                            log_bound_spread,
                                            &log_min_bound,
                                            &log_bound_spread)
                log_density[i] = logaddexp(log_min_bound,
                                           log_bound_spread - log(2))
                pt += n_features

        # normalize the results
        for i in range(log_density.shape[0]):
            log_density[i] += log_knorm

        log_density_arr = log_density_arr.reshape(X.shape[:X.ndim - 1])

        if return_log:
            return log_density_arr
        else:
            return np.exp(log_density_arr)

    def two_point_correlation(self, X, r, dualtree=False):
        """
        two_point_correlation(X, r, dualtree=False)

        Compute the two-point correlation function

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            An array of points to query.  Last dimension should match dimension
            of training data.
        r : array-like
            A one-dimensional array of distances
        dualtree : bool, default=False
            If True, use a dualtree algorithm.  Otherwise, use a single-tree
            algorithm.  Dual tree algorithms can have better scaling for
            large N.

        Returns
        -------
        counts : ndarray
            counts[i] contains the number of pairs of points with distance
            less than or equal to r[i]
        """
        cdef ITYPE_t n_features = self.data.shape[1]
        cdef ITYPE_t i

        # validate X and prepare for query
        X = check_array(X, dtype=DTYPE, order='C')

        if X.shape[X.ndim - 1] != self.data.shape[1]:
            raise ValueError("query data dimension must "
                             "match training data dimension")

        np_Xarr = X.reshape((-1, self.data.shape[1]))
        cdef DTYPE_t[:, ::1] Xarr = get_memview_DTYPE_2D(np_Xarr)

        # prepare r for query
        r = np.asarray(r, dtype=DTYPE, order='C')
        r = np.atleast_1d(r)
        if r.ndim != 1:
            raise ValueError("r must be a 1-dimensional array")
        i_rsort = np.argsort(r)
        rarr_np = r[i_rsort]  # needed to keep memory in scope
        cdef DTYPE_t[::1] rarr = get_memview_DTYPE_1D(rarr_np)

        # create array to hold counts
        count = np.zeros(r.shape[0], dtype=ITYPE)
        cdef ITYPE_t[::1] carr = get_memview_ITYPE_1D(count)

        cdef DTYPE_t* pt = &Xarr[0, 0]

        if dualtree:
            other = self.__class__(Xarr, metric=self.dist_metric,
                                   leaf_size=self.leaf_size)
            self._two_point_dual(0, other, 0, &rarr[0], &carr[0],
                                 0, rarr.shape[0])
        else:
            for i in range(Xarr.shape[0]):
                self._two_point_single(0, pt, &rarr[0], &carr[0],
                                       0, rarr.shape[0])
                pt += n_features

        return count

    cdef int _query_single_depthfirst(self, ITYPE_t i_node,
                                      DTYPE_t* pt, ITYPE_t i_pt,
                                      NeighborsHeap heap,
                                      DTYPE_t reduced_dist_LB) nogil except -1:
        """Recursive Single-tree k-neighbors query, depth-first approach"""
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
                    heap._push(i_pt, dist_pt, self.idx_array[i])

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
        return 0

    cdef int _query_single_breadthfirst(self, DTYPE_t* pt,
                                        ITYPE_t i_pt,
                                        NeighborsHeap heap,
                                        NodeHeap nodeheap) except -1:
        """Non-recursive single-tree k-neighbors query, breadth-first search"""
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
                        heap._push(i_pt, dist_pt, self.idx_array[i])

            #------------------------------------------------------------
            # Case 3: Node is not a leaf.  Add subnodes to the node heap
            else:
                self.n_splits += 1
                for i in range(2 * i_node + 1, 2 * i_node + 3):
                    nodeheap_item.i1 = i
                    nodeheap_item.val = min_rdist(self, i, pt)
                    nodeheap.push(nodeheap_item)
        return 0

    cdef int _query_dual_depthfirst(self, ITYPE_t i_node1,
                                    BinaryTree other, ITYPE_t i_node2,
                                    DTYPE_t[::1] bounds,
                                    NeighborsHeap heap,
                                    DTYPE_t reduced_dist_LB) except -1:
        """Recursive dual-tree k-neighbors query, depth-first"""
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
                        heap._push(i_pt, dist_pt, self.idx_array[i1])

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
        # Case 3a: node 1 is a leaf or is smaller: split node 2 and
        #          recursively query, starting with the nearest subnode
        elif node_info1.is_leaf or (not node_info2.is_leaf
                                    and node_info2.radius > node_info1.radius):
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
        # Case 3b: node 2 is a leaf or is smaller: split node 1 and
        #          recursively query, starting with the nearest subnode
        else:
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
        return 0

    cdef int _query_dual_breadthfirst(self, BinaryTree other,
                                      NeighborsHeap heap,
                                      NodeHeap nodeheap) except -1:
        """Non-recursive dual-tree k-neighbors query, breadth-first"""
        cdef ITYPE_t i, i1, i2, i_node1, i_node2, i_pt
        cdef DTYPE_t dist_pt, reduced_dist_LB
        cdef DTYPE_t[::1] bounds = np.full(other.node_data.shape[0], np.inf)
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
                            heap._push(i_pt, dist_pt, self.idx_array[i1])

                    # keep track of node bound
                    bounds[i_node2] = fmax(bounds[i_node2],
                                           heap.largest(i_pt))

            #------------------------------------------------------------
            # Case 3a: node 1 is a leaf or is smaller: split node 2 and
            #          recursively query, starting with the nearest subnode
            elif node_info1.is_leaf or (not node_info2.is_leaf
                                        and (node_info2.radius
                                             > node_info1.radius)):
                nodeheap_item.i1 = i_node1
                for i2 in range(2 * i_node2 + 1, 2 * i_node2 + 3):
                    nodeheap_item.i2 = i2
                    nodeheap_item.val = min_rdist_dual(self, i_node1,
                                                       other, i2)
                    nodeheap.push(nodeheap_item)

            #------------------------------------------------------------
            # Case 3b: node 2 is a leaf or is smaller: split node 1 and
            #          recursively query, starting with the nearest subnode
            else:
                nodeheap_item.i2 = i_node2
                for i1 in range(2 * i_node1 + 1, 2 * i_node1 + 3):
                    nodeheap_item.i1 = i1
                    nodeheap_item.val = min_rdist_dual(self, i1,
                                                       other, i_node2)
                    nodeheap.push(nodeheap_item)
        return 0

    cdef ITYPE_t _query_radius_single(self,
                                      ITYPE_t i_node,
                                      DTYPE_t* pt, DTYPE_t r,
                                      ITYPE_t* indices,
                                      DTYPE_t* distances,
                                      ITYPE_t count,
                                      int count_only,
                                      int return_distance) nogil:
        """recursive single-tree radius query, depth-first"""
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
                        return -1
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
            reduced_r = self.dist_metric._dist_to_rdist(r)

            for i in range(node_info.idx_start, node_info.idx_end):
                dist_pt = self.rdist(pt, (data + n_features * idx_array[i]),
                                     n_features)
                if dist_pt <= reduced_r:
                    if (count < 0) or (count >= self.data.shape[0]):
                        return -1
                    if count_only:
                        pass
                    else:
                        indices[count] = idx_array[i]
                        if return_distance:
                            distances[count] =\
                                self.dist_metric._rdist_to_dist(dist_pt)
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
                                          DTYPE_t log_knorm,
                                          DTYPE_t log_atol, DTYPE_t log_rtol,
                                          NodeHeap nodeheap,
                                          DTYPE_t* node_log_min_bounds,
                                          DTYPE_t* node_log_bound_spreads):
        """non-recursive single-tree kernel density estimation"""
        # For the given point, node_log_min_bounds and node_log_bound_spreads
        # will encode the current bounds on the density between the point
        # and the associated node.
        # The variables global_log_min_bound and global_log_bound_spread
        # keep track of the global bounds on density.  The procedure here is
        # to split nodes, updating these bounds, until the bounds are within
        # atol & rtol.
        cdef ITYPE_t i, i1, i2, i_node
        cdef DTYPE_t N1, N2
        cdef DTYPE_t global_log_min_bound, global_log_bound_spread
        cdef DTYPE_t global_log_max_bound

        cdef DTYPE_t* data = &self.data[0, 0]
        cdef bint with_sample_weight = self.sample_weight is not None
        cdef DTYPE_t* sample_weight
        if with_sample_weight:
            sample_weight = &self.sample_weight[0]
        cdef ITYPE_t* idx_array = &self.idx_array[0]
        cdef NodeData_t* node_data = &self.node_data[0]
        cdef DTYPE_t N
        cdef DTYPE_t log_weight
        if with_sample_weight:
            N = self.sum_weight
        else:
            N = <DTYPE_t> self.data.shape[0]
        cdef ITYPE_t n_features = self.data.shape[1]

        cdef NodeData_t node_info
        cdef DTYPE_t dist_pt, log_density
        cdef DTYPE_t dist_LB_1 = 0, dist_LB_2 = 0
        cdef DTYPE_t dist_UB_1 = 0, dist_UB_2 = 0

        cdef DTYPE_t dist_UB, dist_LB

        # push the top node to the heap
        cdef NodeHeapData_t nodeheap_item
        nodeheap_item.val = min_dist(self, 0, pt)
        nodeheap_item.i1 = 0
        nodeheap.push(nodeheap_item)

        global_log_min_bound = log(N) + compute_log_kernel(max_dist(self,
                                                                    0, pt),
                                                           h, kernel)
        global_log_max_bound = log(N) + compute_log_kernel(nodeheap_item.val,
                                                           h, kernel)
        global_log_bound_spread = logsubexp(global_log_max_bound,
                                            global_log_min_bound)

        node_log_min_bounds[0] = global_log_min_bound
        node_log_bound_spreads[0] = global_log_bound_spread

        while nodeheap.n > 0:
            nodeheap_item = nodeheap.pop()
            i_node = nodeheap_item.i1

            node_info = node_data[i_node]
            if with_sample_weight:
                N1 = _total_node_weight(node_data, sample_weight,
                                        idx_array, i_node)
            else:
                N1 = node_info.idx_end - node_info.idx_start

            #------------------------------------------------------------
            # Case 1: local bounds are equal to within per-point tolerance.
            if (log_knorm + node_log_bound_spreads[i_node] - log(N1) + log(N)
                <= logaddexp(log_atol, (log_rtol + log_knorm
                                        + node_log_min_bounds[i_node]))):
                pass

            #------------------------------------------------------------
            # Case 2: global bounds are within rtol & atol.
            elif (log_knorm + global_log_bound_spread
                  <= logaddexp(log_atol,
                               log_rtol + log_knorm + global_log_min_bound)):
                break

            #------------------------------------------------------------
            # Case 3: node is a leaf. Count contributions from all points
            elif node_info.is_leaf:
                global_log_min_bound =\
                    logsubexp(global_log_min_bound,
                              node_log_min_bounds[i_node])
                global_log_bound_spread =\
                    logsubexp(global_log_bound_spread,
                              node_log_bound_spreads[i_node])
                for i in range(node_info.idx_start, node_info.idx_end):
                    dist_pt = self.dist(pt, data + n_features * idx_array[i],
                                        n_features)
                    log_density = compute_log_kernel(dist_pt, h, kernel)
                    if with_sample_weight:
                        log_weight = np.log(sample_weight[idx_array[i]])
                    else:
                        log_weight = 0.
                    global_log_min_bound = logaddexp(global_log_min_bound,
                                                     log_density + log_weight)

            #------------------------------------------------------------
            # Case 4: split node and query subnodes
            else:
                i1 = 2 * i_node + 1
                i2 = 2 * i_node + 2

                if with_sample_weight:
                    N1 = _total_node_weight(node_data, sample_weight,
                                            idx_array, i1)
                    N2 = _total_node_weight(node_data, sample_weight,
                                            idx_array, i2)
                else:
                    N1 = node_data[i1].idx_end - node_data[i1].idx_start
                    N2 = node_data[i2].idx_end - node_data[i2].idx_start

                min_max_dist(self, i1, pt, &dist_LB_1, &dist_UB_1)
                min_max_dist(self, i2, pt, &dist_LB_2, &dist_UB_2)

                node_log_min_bounds[i1] = (log(N1) +
                                           compute_log_kernel(dist_UB_1,
                                                              h, kernel))
                node_log_bound_spreads[i1] = (log(N1) +
                                              compute_log_kernel(dist_LB_1,
                                                                 h, kernel))

                node_log_min_bounds[i2] = (log(N2) +
                                           compute_log_kernel(dist_UB_2,
                                                              h, kernel))
                node_log_bound_spreads[i2] = (log(N2) +
                                              compute_log_kernel(dist_LB_2,
                                                                 h, kernel))

                global_log_min_bound = logsubexp(global_log_min_bound,
                                                 node_log_min_bounds[i_node])
                global_log_min_bound = logaddexp(global_log_min_bound,
                                                 node_log_min_bounds[i1])
                global_log_min_bound = logaddexp(global_log_min_bound,
                                                 node_log_min_bounds[i2])

                global_log_bound_spread =\
                    logsubexp(global_log_bound_spread,
                              node_log_bound_spreads[i_node])
                global_log_bound_spread = logaddexp(global_log_bound_spread,
                                                    node_log_bound_spreads[i1])
                global_log_bound_spread = logaddexp(global_log_bound_spread,
                                                    node_log_bound_spreads[i2])

                # TODO: rank by the spread rather than the distance?
                nodeheap_item.val = dist_LB_1
                nodeheap_item.i1 = i1
                nodeheap.push(nodeheap_item)

                nodeheap_item.val = dist_LB_2
                nodeheap_item.i1 = i2
                nodeheap.push(nodeheap_item)

        nodeheap.clear()
        return logaddexp(global_log_min_bound,
                         global_log_bound_spread - log(2))

    cdef int _kde_single_depthfirst(
                   self, ITYPE_t i_node, DTYPE_t* pt,
                   KernelType kernel, DTYPE_t h,
                   DTYPE_t log_knorm,
                   DTYPE_t log_atol, DTYPE_t log_rtol,
                   DTYPE_t local_log_min_bound,
                   DTYPE_t local_log_bound_spread,
                   DTYPE_t* global_log_min_bound,
                   DTYPE_t* global_log_bound_spread) except -1:
        """recursive single-tree kernel density estimate, depth-first"""
        # For the given point, local_min_bound and local_max_bound give the
        # minimum and maximum density for the current node, while
        # global_min_bound and global_max_bound give the minimum and maximum
        # density over the entire tree.  We recurse down until global_min_bound
        # and global_max_bound are within rtol and atol.
        cdef ITYPE_t i, i1, i2, iw, start, end
        cdef DTYPE_t N1, N2

        cdef DTYPE_t* data = &self.data[0, 0]
        cdef NodeData_t* node_data = &self.node_data[0]
        cdef bint with_sample_weight = self.sample_weight is not None
        cdef DTYPE_t* sample_weight
        cdef DTYPE_t log_weight
        if with_sample_weight:
            sample_weight = &self.sample_weight[0]
        cdef ITYPE_t* idx_array = &self.idx_array[0]
        cdef ITYPE_t n_features = self.data.shape[1]

        cdef NodeData_t node_info = self.node_data[i_node]
        cdef DTYPE_t dist_pt, log_dens_contribution

        cdef DTYPE_t child1_log_min_bound, child2_log_min_bound
        cdef DTYPE_t child1_log_bound_spread, child2_log_bound_spread
        cdef DTYPE_t dist_UB = 0, dist_LB = 0

        if with_sample_weight:
            N1  = _total_node_weight(node_data, sample_weight,
                                     idx_array, i_node)
            N2 = self.sum_weight
        else:
            N1 = <DTYPE_t>(node_info.idx_end - node_info.idx_start)
            N2 = <DTYPE_t>self.data.shape[0]

        #------------------------------------------------------------
        # Case 1: local bounds are equal to within errors.  Return
        if (log_knorm + local_log_bound_spread - log(N1) + log(N2)
            <= logaddexp(log_atol, (log_rtol + log_knorm
                                    + local_log_min_bound))):
            pass

        #------------------------------------------------------------
        # Case 2: global bounds are within rtol & atol. Return
        elif (log_knorm + global_log_bound_spread[0]
            <= logaddexp(log_atol, (log_rtol + log_knorm
                                    + global_log_min_bound[0]))):
            pass

        #------------------------------------------------------------
        # Case 3: node is a leaf. Count contributions from all points
        elif node_info.is_leaf:
            global_log_min_bound[0] = logsubexp(global_log_min_bound[0],
                                                local_log_min_bound)
            global_log_bound_spread[0] = logsubexp(global_log_bound_spread[0],
                                                   local_log_bound_spread)
            for i in range(node_info.idx_start, node_info.idx_end):
                dist_pt = self.dist(pt, (data + n_features * idx_array[i]),
                                    n_features)
                log_dens_contribution = compute_log_kernel(dist_pt, h, kernel)
                if with_sample_weight:
                    log_weight = np.log(sample_weight[idx_array[i]])
                else:
                    log_weight = 0.
                global_log_min_bound[0] = logaddexp(global_log_min_bound[0],
                                                    (log_dens_contribution +
                                                     log_weight))

        #------------------------------------------------------------
        # Case 4: split node and query subnodes
        else:
            i1 = 2 * i_node + 1
            i2 = 2 * i_node + 2

            if with_sample_weight:
                N1 = _total_node_weight(node_data, sample_weight,
                                        idx_array, i1)
                N2 = _total_node_weight(node_data, sample_weight,
                                        idx_array, i2)
            else:
                N1 = <DTYPE_t>(self.node_data[i1].idx_end - self.node_data[i1].idx_start)
                N2 = <DTYPE_t>(self.node_data[i2].idx_end - self.node_data[i2].idx_start)

            min_max_dist(self, i1, pt, &dist_LB, &dist_UB)
            child1_log_min_bound = log(N1) + compute_log_kernel(dist_UB, h,
                                                                kernel)
            child1_log_bound_spread = logsubexp(log(N1) +
                                                compute_log_kernel(dist_LB, h,
                                                                   kernel),
                                                child1_log_min_bound)

            min_max_dist(self, i2, pt, &dist_LB, &dist_UB)
            child2_log_min_bound = log(N2) + compute_log_kernel(dist_UB, h,
                                                                kernel)
            child2_log_bound_spread = logsubexp(log(N2) +
                                                compute_log_kernel(dist_LB, h,
                                                                   kernel),
                                                child2_log_min_bound)

            global_log_min_bound[0] = logsubexp(global_log_min_bound[0],
                                                local_log_min_bound)
            global_log_min_bound[0] = logaddexp(global_log_min_bound[0],
                                                child1_log_min_bound)
            global_log_min_bound[0] = logaddexp(global_log_min_bound[0],
                                                child2_log_min_bound)

            global_log_bound_spread[0] = logsubexp(global_log_bound_spread[0],
                                                   local_log_bound_spread)
            global_log_bound_spread[0] = logaddexp(global_log_bound_spread[0],
                                                   child1_log_bound_spread)
            global_log_bound_spread[0] = logaddexp(global_log_bound_spread[0],
                                                   child2_log_bound_spread)

            self._kde_single_depthfirst(i1, pt, kernel, h, log_knorm,
                                        log_atol, log_rtol,
                                        child1_log_min_bound,
                                        child1_log_bound_spread,
                                        global_log_min_bound,
                                        global_log_bound_spread)
            self._kde_single_depthfirst(i2, pt, kernel, h, log_knorm,
                                        log_atol, log_rtol,
                                        child2_log_min_bound,
                                        child2_log_bound_spread,
                                        global_log_min_bound,
                                        global_log_bound_spread)
        return 0

    cdef int _two_point_single(self, ITYPE_t i_node, DTYPE_t* pt, DTYPE_t* r,
                               ITYPE_t* count, ITYPE_t i_min,
                               ITYPE_t i_max) except -1:
        """recursive single-tree two-point correlation function query"""
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
        return 0

    cdef int _two_point_dual(self, ITYPE_t i_node1,
                             BinaryTree other, ITYPE_t i_node2,
                             DTYPE_t* r, ITYPE_t* count,
                             ITYPE_t i_min, ITYPE_t i_max) except -1:
        """recursive dual-tree two-point correlation function query"""
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
        return 0


######################################################################
# Python functions for benchmarking and testing C implementations

def load_heap(DTYPE_t[:, ::1] X, ITYPE_t k):
    """test fully loading the heap"""
    assert k <= X.shape[1]
    cdef NeighborsHeap heap = NeighborsHeap(X.shape[0], k)
    cdef ITYPE_t i, j
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            heap._push(i, X[i, j], j)
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


cdef inline DTYPE_t _total_node_weight(NodeData_t* node_data,
                                       DTYPE_t* sample_weight,
                                       ITYPE_t* idx_array,
                                       ITYPE_t i_node):
    cdef ITYPE_t i
    cdef DTYPE_t N = 0.0
    for i in range(node_data[i_node].idx_start, node_data[i_node].idx_end):
        N += sample_weight[idx_array[i]]
    return N
