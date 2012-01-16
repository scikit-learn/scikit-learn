
# Author: Jake Vanderplas <vanderplas@astro.washington.edu>
# License: BSD

"""
=========
Ball Tree
=========
A ball tree is a data object which speeds up nearest neighbor
searches in high dimensions (see scikit-learn neighbors module
documentation for an overview of neighbor trees). There are many
types of ball trees.  This package provides a basic implementation
in cython.

Implementation Notes
--------------------

A ball tree can be thought of as a collection of nodes.  Each node
stores a centroid, a radius, and the pointers to two child nodes.

* centroid : the centroid of a node is the mean of all the locations
    of points within the node
* radius : the radius of a node is the distance from the centroid
    to the furthest point in the node.
* subnodes : each node has a maximum of 2 child nodes.  The data within
    the parent node is divided between the two child nodes.

In a typical tree implementation, nodes may be classes or structures which
are dynamically allocated as needed.  This offers flexibility in the number
of nodes, and leads to very straightforward and readable code.  It also means
that the tree can be dynamically augmented or pruned with new data, in an
in-line fashion.  This approach generally leads to recursive code: upon
construction, the head node constructs its child nodes, the child nodes
construct their child nodes, and so-on.

The current package uses a different approach: all node data is stored in
a set of numpy arrays which are pre-allocated.  The main advantage of this
approach is that the whole object can be quickly and easily saved to disk
and reconstructed from disk.  This also allows for an iterative interface
which gives more control over the heap, and leads to speed.  There are a
few disadvantages, however: once the tree is built, augmenting or pruning it
is not as straightforward.  Also, the size of the tree must be known from the
start, so there is not as much flexibility in building it.

BallTree Pseudo-code
~~~~~~~~~~~~~~~~~~~~
Because understanding a ball tree is simpler with recursive code, here is some
pseudo-code to show the structure of the main functionality

    # Ball Tree pseudo code

    class Node:
        #class data:
        centroid
        radius
        child1, child2

        #class methods:
        def construct(data):
            centroid = compute_centroid(data)
            radius = compute_radius(centroid, data)

            # Divide the data into two approximately equal sets.
            # This is often done by splitting along a single dimension.
            data1, data2 = divide(data)

            if number_of_points(data1) > 0:
                child1.construct(data1)

            if number_of_points(data2) > 0:
                child2.construct(data2)

        def query(pt, neighbors_heap):
            # compute the minimum distance from pt to any point in this node
            d = distance(point, centroid)
            if d < radius:
                min_distance = 0
            else:
                min_distance = d - radius

            if min_distance > max_distance_in(neighbors_heap):
                # all these points are too far away.  cut off the search here
                return
            elif node_size > 1:
                child1.query(pt, neighbors_heap)
                child2.query(pt, neighbors_heap)


    object BallTree:
        #class data:
        data
        root_node

        #class methods
        def construct(data, num_leaves):
            root_node.construct(data)

        def query(point, num_neighbors):
            neighbors_heap = empty_heap_of_size(num_neighbors)
            root_node.query(point, neighbors_heap)

This certainly is not a complete description, but should give the basic idea
of the form of the algorithm.  The implementation below is much faster than
anything mirroring the pseudo-code above, but for that reason is much more
opaque.  Here's the basic idea:

BallTree Storage
~~~~~~~~~~~~~~~~
The BallTree information is stored using a combination of
"Array of Structures" and "Structure of Arrays" to maximize speed.
Given input data of size ``(n_samples, n_features)``, BallTree computes the
expected number of nodes ``n_nodes`` (see below), and allocates the
following arrays:

* ``data`` : a float array of shape ``(n_samples, n_features)``
    This is simply the input data.  If the input matrix is well-formed
    (contiguous, c-ordered, correct data type) then no copy is needed
* ``idx_array`` : an integer array of size ``n_samples``
    This can be thought of as an array of pointers to the data in ``data``.
    Rather than shuffling around the data itself, we shuffle around pointers
    to the rows in data.
* ``node_centroid_arr`` : a float array of shape ``(n_nodes, n_features)``
    This stores the centroid of the data in each node.
* ``node_info_arr`` : a size-``n_nodes`` array of ``NodeInfo`` structures.
    This stores information associated with each node.  Each ``NodeInfo``
    instance has the following attributes:
    - ``idx_start``
    - ``idx_end`` : ``idx_start`` and ``idx_end`` reference the part of
      ``idx_array`` which point to the data associated with the node.
      The data in node with index ``i_node`` is given by
      ``data[idx_array[idx_start:idx_end]]``
    - ``is_leaf`` : a boolean value which tells whether this node is a leaf:
      that is, whether or not it has children.
    - ``radius`` : a floating-point value which gives the distance from
      the node centroid to the furthest point in the node.

One feature here is that there are no stored pointers from parent nodes to
child nodes and vice-versa.  These pointers are implemented implicitly:
For a node with index ``i``, the two children are found at indices
``2 * i + 1`` and ``2 * i + 2``, while the parent is found at index
``floor((i - 1) / 2)``.  The root node has no parent.

With this data structure in place, the functionality of the above BallTree
pseudo-code can be implemented in a much more efficient manner.
Most of the data passing done in this code uses raw data pointers.
Using numpy arrays would be preferable for safety, but the
overhead of array slicing and sub-array construction leads to execution
time which is several orders of magnitude slower than the current
implementation.

Priority Queue vs Max-heap
~~~~~~~~~~~~~~~~~~~~~~~~~~
When querying for more than one neighbor, the code must maintain a list of
the current k nearest points.  The BallTree code implements this in two ways.

- A priority queue: this is just a sorted list.  When an item is added,
  it is inserted in the appropriate location.  The cost of the search plus
  insert averages O[k].
- A max-heap: this is a binary tree structure arranged such that each node is
  greater than its children.  The cost of adding an item is O[log(k)].
  At the end of the iterations, the results must be sorted: a quicksort is
  used, which averages O[k log(k)].  Quicksort has worst-case O[k^2]
  performance, but because the input is already structured in a max-heap,
  the worst case will not be realized.  Thus the sort is a one-time operation
  with cost O[k log(k)].

Each insert is performed an average of log(N) times per query, where N is
the number of training points.  Because of this, for a single query, the
priority-queue approach costs O[k log(N)], and the max-heap approach costs
O[log(k)log(N)] + O[k log(k)].  Tests show that for sufficiently large k,
the max-heap approach out-performs the priority queue approach by a factor
of a few.  In light of these tests, the code uses a priority queue for
k < 5, and a max-heap otherwise.

Memory Allocation
~~~~~~~~~~~~~~~~~
It is desirable to construct a tree in as balanced a way as possible.
Given a training set with n_samples and a user-supplied leaf_size, if
the points in each node are divided as evenly as possible between the
two children, the maximum depth needed so that leaf nodes satisfy
``leaf_size <= n_points <= 2 * leaf_size`` is given by
``n_levels = 1 + max(0, floor(log2((n_samples - 1) / leaf_size)))``
(with the exception of the special case where ``n_samples < leaf_size``)
For a given number of levels, the number of points in a tree is given by
``n_nodes = 2 ** n_levels - 1``.  Both of these results can be shown
by induction.  Using them, the correct amount of memory can be pre-allocated
for a given ``n_samples`` and ``leaf_size``.
"""

import numpy as np

cimport numpy as np
cimport cython
from libc cimport stdlib

from ..utils import array2d

######################################################################
# global definitions
#
# type used for data
DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

# type used for indices & counts
# warning: there will be problems if this is switched to an unsigned type!
ITYPE = np.int32
ctypedef np.int32_t ITYPE_t

# infinity
cdef DTYPE_t infinity = np.inf


######################################################################
# utility functions: fast max, min, and absolute value
#
@cython.profile(False)
cdef inline DTYPE_t dmax(DTYPE_t x, DTYPE_t y):
    if x >= y:
        return x
    else:
        return y


@cython.profile(False)
cdef inline DTYPE_t dmin(DTYPE_t x, DTYPE_t y):
    if x <= y:
        return x
    else:
        return y


@cython.profile(False)
cdef inline DTYPE_t dabs(DTYPE_t x):
    if x >= 0:
        return x
    else:
        return -x


######################################################################
# distance functions
#  These implement the Minkowski p-distance given by
#    dist = sum((x - y) ** p) ** (1 / p)
#  To compare distances, the raising to the (1 / p) is not necessary
#  therefore, for speed, we also define a function dist_p() given by
#    dist_p = sum((x - y) ** p)
#  there are also functions dist_from_dist_p() and dist_p_from_dist()
#  which convert between these forms.
@cython.cdivision(True)
cdef DTYPE_t dist(DTYPE_t *x1, DTYPE_t *x2, ITYPE_t n, DTYPE_t p):
    cdef ITYPE_t i
    cdef DTYPE_t r, d
    r = 0
    if p == 2:
        for i from 0 <= i < n:
            d = x1[i] - x2[i]
            r += d * d
        r = r ** 0.5
    elif p == infinity:
        for i from 0 <= i < n:
            r = dmax(r, dabs(x1[i] - x2[i]))
    elif p == 1:
        for i from 0 <= i < n:
            r += dabs(x1[i] - x2[i])
    else:
        for i from 0 <= i < n:
            d = dabs(x1[i] - x2[i])
            r += d ** p
        r = r ** (1. / p)
    return r


@cython.cdivision(True)
cdef DTYPE_t dist_p(DTYPE_t *x1, DTYPE_t *x2, ITYPE_t n, DTYPE_t p):
    cdef ITYPE_t i
    cdef DTYPE_t r, d
    r = 0
    if p == 2:
        for i from 0 <= i < n:
            d = x1[i] - x2[i]
            r += d * d
    elif p == infinity:
        for i from 0 <= i < n:
            r = dmax(r, dabs(x1[i] - x2[i]))
    elif p == 1:
        for i from 0 <= i < n:
            r += dabs(x1[i] - x2[i])
    else:
        for i from 0 <= i < n:
            d = dabs(x1[i] - x2[i])
            r += d ** p
    return r


@cython.cdivision(True)
cdef DTYPE_t dist_from_dist_p(DTYPE_t r, DTYPE_t p):
    if p == 2:
        return r ** 0.5
    elif p == infinity:
        return r
    elif p == 1:
        return r
    else:
        return r ** (1. / p)


@cython.cdivision(True)
cdef DTYPE_t dist_p_from_dist(DTYPE_t r, DTYPE_t p):
    if p == 2:
        return r ** 2
    elif p == infinity:
        return r
    elif p == 1:
        return r
    else:
        return r ** p


######################################################################
# NodeInfo struct
#  used to keep track of node information.
#  there is also a centroid for each node: this is kept in a separate
#  array for efficiency.  This is a hybrid of the "Array of Structures"
#  and "Structure of Arrays" styles.
cdef struct NodeInfo:
    ITYPE_t idx_start
    ITYPE_t idx_end
    ITYPE_t is_leaf
    DTYPE_t radius


######################################################################
# stack struct
#  This is used to keep track of the recursion stack in Node_query
cdef struct stack_item:
    ITYPE_t i_node
    DTYPE_t dist_p_LB


cdef struct stack:
    int n
    stack_item* heap
    int size


@cython.profile(False)
cdef inline void stack_create(stack* self, int size):
    self.size = size
    self.heap = <stack_item*> stdlib.malloc(sizeof(stack_item) * size)
    self.n = 0


@cython.profile(False)
cdef inline void stack_destroy(stack* self):
    stdlib.free(self.heap)


@cython.profile(False)
cdef inline void stack_resize(stack* self, int new_size):
    #print "resize", self.n, new_size
    if new_size < self.n:
        raise ValueError("new_size smaller than current")

    self.size = new_size
    self.heap = <stack_item*>stdlib.realloc(<void*> self.heap,
                                            new_size * sizeof(stack_item))


@cython.profile(False)
cdef inline void stack_push(stack* self, stack_item item):
    if self.n >= self.size:
        stack_resize(self, 2 * self.size + 1)

    self.heap[self.n] = item
    self.n += 1


@cython.profile(False)
cdef inline stack_item stack_pop(stack* self):
    if self.n == 0:
        raise ValueError("popping empty stack")

    self.n -= 1
    return self.heap[self.n]


######################################################################
# newObj function
#  this is a helper function for pickling
def newObj(obj):
    return obj.__new__(obj)


######################################################################
# BallTree class
#
cdef class BallTree(object):
    """
    Ball Tree for fast nearest-neighbor searches :

    BallTree(X, leaf_size=20, p=2.0)

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
        to store the built ball tree.  The amount of memory needed to
        store the tree scales as
        2 ** (1 + floor(log2((n_samples - 1) / leaf_size))) - 1
        For a specified ``leaf_size``, a leaf node is guaranteed to
        satisfy ``leaf_size <= n_points <= 2 * leaf_size``, except in
        the case that ``n_samples < leaf_size``.

    p : distance metric for the BallTree.  ``p`` encodes the Minkowski
        p-distance::

            D = sum((X[i] - X[j]) ** p) ** (1. / p)

        p must be greater than or equal to 1, so that the triangle
        inequality will hold.  If ``p == np.inf``, then the distance is
        equivalent to::

            D = max(X[i] - X[j])

    Attributes
    ----------
    data : np.ndarray
        The training data

    warning_flag : bool
        Warning flag is set to true during query(...) if results are
        dependent on the order of the training cases.
        For classification or regression based on k-neighbors, if
        neighbor k and neighbor k+1 have identical distances but different
        labels, then the result will be dependent on the ordering of the
        training data.  In this case, ``warning_flag`` will be set to True.

    Examples
    --------
    Query for k-nearest neighbors

        >>> import numpy as np
        >>> np.random.seed(0)
        >>> X = np.random.random((10,3))  # 10 points in 3 dimensions
        >>> ball_tree = BallTree(X, leaf_size=2)              # doctest: +SKIP
        >>> dist, ind = ball_tree.query(X[0], n_neighbors=3)  # doctest: +SKIP
        >>> print ind  # indices of 3 closest neighbors
        [0 3 1]
        >>> print dist  # distances to 3 closest neighbors
        [ 0.          0.19662693  0.29473397]

    Pickle and Unpickle a ball tree (using protocol = 2).  Note that the
    state of the tree is saved in the pickle operation: the tree is not
    rebuilt on un-pickling

        >>> import numpy as np
        >>> import pickle
        >>> np.random.seed(0)
        >>> X = np.random.random((10,3))  # 10 points in 3 dimensions
        >>> ball_tree = BallTree(X, leaf_size=2)          # doctest: +SKIP
        >>> s = pickle.dumps(ball_tree, protocol=2)       # doctest: +SKIP
        >>> ball_tree_copy = pickle.loads(s)              # doctest: +SKIP
        >>> dist, ind = ball_tree_copy.query(X[0], k=3)   # doctest: +SKIP
        >>> print ind  # indices of 3 closest neighbors   
        [0 3 1]
        >>> print dist  # distances to 3 closest neighbors
        [ 0.          0.19662693  0.29473397]
    """
    cdef readonly np.ndarray data
    cdef np.ndarray idx_array
    cdef np.ndarray node_centroid_arr
    cdef np.ndarray node_info_arr
    cdef DTYPE_t p
    cdef ITYPE_t leaf_size
    cdef ITYPE_t n_levels
    cdef ITYPE_t n_nodes
    cdef readonly int warning_flag

    def __cinit__(self):
        """
        initialize all arrays to empty.  This will prevent memory errors
        in rare cases where __init__ is not called
        """
        self.data = np.empty((0,0), dtype=DTYPE)
        self.idx_array = np.empty(0, dtype=ITYPE)
        self.node_centroid_arr = np.empty((0,0), dtype=DTYPE)
        self.node_info_arr = np.empty(0, dtype='c')

    def __init__(self, X, ITYPE_t leaf_size=20, DTYPE_t p=2):
        self.data = np.asarray(X, dtype=DTYPE, order='C')
        self.warning_flag = True

        if X.size == 0:
            raise ValueError("X is an empty array")

        if self.data.ndim != 2:
            raise ValueError("X should have two dimensions")

        if p < 1:
            raise ValueError("p must be greater than or equal to 1")
        self.p = p

        if leaf_size < 1:
            raise ValueError("leaf_size must be greater than or equal to 1")
        self.leaf_size = leaf_size

        cdef ITYPE_t n_samples = self.data.shape[0]
        cdef ITYPE_t n_features = self.data.shape[1]

        # determine number of levels in the ball tree, and from this
        # the number of nodes in the ball tree
        self.n_levels = np.log2(max(1, (n_samples - 1)/self.leaf_size)) + 1
        self.n_nodes = (2 ** self.n_levels) - 1

        self.idx_array = np.arange(n_samples, dtype=ITYPE)

        self.node_centroid_arr = np.empty((self.n_nodes, n_features),
                                          dtype=DTYPE, order='C')

        self.node_info_arr = np.empty(self.n_nodes * sizeof(NodeInfo),
                                      dtype='c', order='C')
        self.build_tree_()

    def __reduce__(self):
        """
        reduce method used for pickling
        """
        return (newObj, (BallTree,), self.__getstate__())

    def __getstate__(self):
        """
        get state for pickling
        """
        return (self.data,
                self.idx_array,
                self.node_centroid_arr,
                self.node_info_arr,
                self.p,
                self.leaf_size,
                self.n_levels,
                self.n_nodes)

    def __setstate__(self, state):
        """
        set state for pickling
        """
        self.data = state[0]
        self.idx_array = state[1]
        self.node_centroid_arr = state[2]
        self.node_info_arr = state[3]
        self.p = state[4]
        self.leaf_size = state[5]
        self.n_levels = state[6]
        self.n_nodes = state[7]

    def query(self, X, k=1, return_distance=True):
        """
        query(X, k=1, return_distance=True)

        query the Ball Tree for the k nearest neighbors

        Parameters
        ----------
        X : array-like, last dimension self.dim
            An array of points to query
        k : integer  (default = 1)
            The number of nearest neighbors to return
        return_distance : boolean (default = True)
            if True, return a tuple (d,i)
            if False, return array i

        Returns
        -------
        i    : if return_distance == False
        (d,i) : if return_distance == True

        d : array of doubles - shape: x.shape[:-1] + (k,)
            each entry gives the list of distances to the
            neighbors of the corresponding point
            (note that distances are not sorted)

        i : array of integers - shape: x.shape[:-1] + (k,)
            each entry gives the list of indices of
            neighbors of the corresponding point
            (note that neighbors are not sorted)

        Examples
        --------
        Query for k-nearest neighbors

            >>> import numpy as np
            >>> np.random.seed(0)
            >>> X = np.random.random((10,3))  # 10 points in 3 dimensions
            >>> ball_tree = BallTree(X, leaf_size=2)    # doctest: +SKIP
            >>> dist, ind = ball_tree.query(X[0], k=3)  # doctest: +SKIP
            >>> print ind  # indices of 3 closest neighbors
            [0 3 1]
            >>> print dist  # distances to 3 closest neighbors
            [ 0.          0.19662693  0.29473397]
        """
        self.warning_flag = False

        X = array2d(X, dtype=DTYPE, order='C')

        if X.shape[-1] != self.data.shape[1]:
            raise ValueError("query data dimension must match BallTree "
                             "data dimension")

        if k > self.data.shape[0]:
            raise ValueError("k must be less than or equal "
                             "to the number of training points")

        # flatten X for iteration
        orig_shape = X.shape
        X = X.reshape((-1, X.shape[-1]))

        # for k less than 5, a priority queue is slightly faster
        # for more neighbors, a max-heap implementation is faster
        cdef ITYPE_t use_max_heap = (k >= 5)

        cdef ITYPE_t i
        cdef ITYPE_t n_neighbors = k
        cdef np.ndarray distances = np.empty((X.shape[0], n_neighbors),
                                             dtype=DTYPE)
        cdef np.ndarray idx_array = np.empty((X.shape[0], n_neighbors),
                                             dtype=ITYPE)
        cdef np.ndarray Xi

        distances[:] = np.inf

        cdef DTYPE_t* dist_ptr = <DTYPE_t*> distances.data
        cdef ITYPE_t* idx_ptr = <ITYPE_t*> idx_array.data

        cdef stack node_stack
        stack_create(&node_stack, self.n_levels + 1)

        for i, Xi in enumerate(X):
            self.query_one_(<DTYPE_t*>Xi.data, n_neighbors,
                            dist_ptr, idx_ptr, &node_stack, use_max_heap)

            # if max-heap is used, results must be sorted
            if use_max_heap:
                sort_dist_idx(dist_ptr, idx_ptr, n_neighbors)

            dist_ptr += n_neighbors
            idx_ptr += n_neighbors

        stack_destroy(&node_stack)

        # deflatten results
        if return_distance:
            return (distances.reshape((orig_shape[:-1]) + (k,)),
                    idx_array.reshape((orig_shape[:-1]) + (k,)))
        else:
            return idx_array.reshape((orig_shape[:-1]) + (k,))

    def query_radius(self, X, r, return_distance=False,
                     count_only=False, sort_results=False):
        """
        query_radius(self, X, r, count_only = False):

        query the Ball Tree for neighbors within a ball of size r

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
            Note that unlike BallTree.query(), setting return_distance=True
            adds to the computation time.  Not all distances need to be
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
            the results of BallTree.query(), the returned neighbors
            are not sorted by distance

        dist : array of objects, shape = X.shape[:-1]
            each element is a numpy double array
            listing the distances corresponding to indices in i.

        Examples
        --------
        Query for neighbors in a given radius

        >>> import numpy as np
        >>> np.random.seed(0)
        >>> X = np.random.random((10,3))  # 10 points in 3 dimensions
        >>> ball_tree = BallTree(X, leaf_size=2)        # doctest: +SKIP
        >>> print ball_tree.query_radius(X[0], r=0.3, count_only=True)
        3
        >>> ind = ball_tree.query_radius(X[0], r=0.3)  # doctest: +SKIP
        >>> print ind  # indices of neighbors within distance 0.3
        [3 0 1]
        """
        if count_only and return_distance:
            raise ValueError("count_only and return_distance "
                             "cannot both be true")

        if sort_results and not return_distance:
            raise ValueError("return_distance must be True if sort_distances "
                             "is True")

        cdef np.ndarray idx_array, idx_array_i, distances, distances_i
        cdef np.ndarray pt, count
        cdef ITYPE_t count_i

        # prepare X for query
        X = array2d(X, dtype=DTYPE, order='C')
        if X.shape[-1] != self.data.shape[1]:
            raise ValueError("query data dimension must match BallTree "
                             "data dimension")

        # prepare r for query
        r = np.asarray(r, dtype=DTYPE, order='C')
        r = np.atleast_1d(r)
        if r.shape == (1,):
            r = r[0] * np.ones(X.shape[:-1], dtype=np.double)
        else:
            if r.shape != X.shape[:-1]:
                raise ValueError("r must be broadcastable to X.shape")

        # flatten X and r for iteration
        orig_shape = X.shape
        X = X.reshape((-1, X.shape[-1]))
        r = r.reshape(-1)

        cdef stack node_stack
        stack_create(&node_stack, self.n_levels + 1)

        if count_only:
            count = np.zeros(X.shape[0], ITYPE)
            for pt_idx, pt in enumerate(X):
                count[pt_idx] = self.query_radius_count_(<DTYPE_t*>pt.data,
                                                         r[pt_idx],
                                                         &node_stack)
        elif not return_distance:
            idx_array = np.empty(X.shape[0], dtype='object')
            idx_array_i = np.empty(self.data.shape[0], dtype=ITYPE)
            for pt_idx, pt in enumerate(X):
                count_i = self.query_radius_idx_only_(
                    <DTYPE_t*>pt.data,
                    r[pt_idx],
                    <ITYPE_t*>idx_array_i.data,
                    &node_stack)
                idx_array[pt_idx] = idx_array_i[:count_i].copy()

        else:
            idx_array = np.empty(X.shape[0], dtype='object')
            distances = np.empty(X.shape[0], dtype='object')
            idx_array_i = np.empty(self.data.shape[0], dtype=ITYPE)
            distances_i = np.empty(self.data.shape[0], dtype=DTYPE)
            for pt_idx, pt in enumerate(X):
                count_i = self.query_radius_distances_(
                    <DTYPE_t*>pt.data,
                    r[pt_idx],
                    <ITYPE_t*>idx_array_i.data,
                    <DTYPE_t*>distances_i.data,
                    &node_stack)
                if sort_results:
                    sort_dist_idx(<DTYPE_t*>distances_i.data,
                                  <ITYPE_t*>idx_array_i.data,
                                  count_i)

                idx_array[pt_idx] = idx_array_i[:count_i].copy()
                distances[pt_idx] = distances_i[:count_i].copy()

        stack_destroy(&node_stack)

        # deflatten results
        if count_only:
            return count.reshape(orig_shape[:-1])
        elif return_distance:
            return (idx_array.reshape(orig_shape[:-1]),
                    distances.reshape(orig_shape[:-1]))
        else:
            return idx_array.reshape(orig_shape[:-1])

    @cython.cdivision(True)
    cdef void build_tree_(BallTree self):
        cdef DTYPE_t* data = <DTYPE_t*> self.data.data
        cdef ITYPE_t* idx_array = <ITYPE_t*> self.idx_array.data
        cdef DTYPE_t* node_centroid_arr = <DTYPE_t*>self.node_centroid_arr.data
        cdef NodeInfo* node_info_arr = <NodeInfo*> self.node_info_arr.data

        cdef DTYPE_t p = self.p
        cdef ITYPE_t n_samples = self.data.shape[0]
        cdef ITYPE_t n_features = self.data.shape[1]

        cdef ITYPE_t idx_start, idx_end, n_points
        cdef DTYPE_t radius
        cdef ITYPE_t i, i_node, i_parent

        cdef DTYPE_t* centroid = node_centroid_arr
        cdef NodeInfo* node_info = node_info_arr
        cdef NodeInfo* parent_info
        cdef DTYPE_t* point

        #------------------------------------------------------------
        # take care of the root node
        node_info.idx_start = 0
        node_info.idx_end = n_samples
        n_points = n_samples

        # determine Node centroid
        compute_centroid(centroid, data, idx_array,
                         n_features, n_samples)

        # determine Node radius
        radius = 0
        for i from node_info.idx_start <= i < node_info.idx_end:
            radius = dmax(radius,
                          dist_p(centroid, data + n_features * idx_array[i],
                                 n_features, p))
        node_info.radius = dist_from_dist_p(radius, p)

        # check if this is a leaf
        if self.n_nodes == 1:
            node_info.is_leaf = 1

        else:
            node_info.is_leaf = 0

            # find dimension with largest spread
            i_max = find_split_dim(data, idx_array + node_info.idx_start,
                                   n_features, n_points)

            # sort idx_array along this dimension
            partition_indices(data,
                              idx_array + node_info.idx_start,
                              i_max,
                              n_points / 2,
                              n_features,
                              n_points)

        #------------------------------------------------------------
        # cycle through all child nodes
        for i_node from 1 <= i_node < self.n_nodes:
            i_parent = (i_node - 1) / 2
            parent_info = node_info_arr + i_parent

            node_info = node_info_arr + i_node

            if parent_info.is_leaf:
                raise ValueError("Fatal: parent is a leaf. Memory "
                                 "allocation is flawed")

            if i_node < self.n_nodes / 2:
                node_info.is_leaf = 0
            else:
                node_info.is_leaf = 1

            centroid = node_centroid_arr + i_node * n_features

            # find indices for this node
            idx_start = parent_info.idx_start
            idx_end = parent_info.idx_end

            if i_node % 2 == 1:
                idx_start = (idx_start + idx_end) / 2
            else:
                idx_end = (idx_start + idx_end) / 2

            node_info.idx_start = idx_start
            node_info.idx_end = idx_end

            n_points = idx_end - idx_start

            if n_points == 0:
                raise ValueError("zero-sized node")

            elif n_points == 1:
                #copy this point to centroid
                copy_array(centroid,
                           data + idx_array[idx_start] * n_features,
                           n_features)

                #store radius in array
                node_info.radius = 0

                #is a leaf
                node_info.is_leaf = 1

            else:
                # determine Node centroid
                compute_centroid(centroid, data, idx_array + idx_start,
                                 n_features, n_points)

                # determine Node radius
                radius = 0
                for i from idx_start <= i < idx_end:
                    radius = dmax(radius,
                                  dist_p(centroid,
                                         data + n_features * idx_array[i],
                                         n_features, p))
                node_info.radius = dist_from_dist_p(radius, p)

                if not node_info.is_leaf:
                    # find dimension with largest spread
                    i_max = find_split_dim(data, idx_array + idx_start,
                                           n_features, n_points)

                    # sort indices along this dimension
                    partition_indices(data,
                                      idx_array + idx_start,
                                      i_max,
                                      n_points / 2,
                                      n_features,
                                      n_points)

    cdef void query_one_(BallTree self,
                         DTYPE_t* pt,
                         ITYPE_t k,
                         DTYPE_t* near_set_dist,
                         ITYPE_t* near_set_indx,
                         stack* node_stack,
                         ITYPE_t use_max_heap):
        cdef DTYPE_t* data = <DTYPE_t*> self.data.data
        cdef ITYPE_t* idx_array = <ITYPE_t*> self.idx_array.data
        cdef DTYPE_t* node_centroid_arr = <DTYPE_t*>self.node_centroid_arr.data
        cdef NodeInfo* node_info_arr = <NodeInfo*> self.node_info_arr.data
        cdef NodeInfo* node_info = node_info_arr

        cdef DTYPE_t p = self.p
        cdef ITYPE_t n_features = self.data.shape[1]

        cdef DTYPE_t dmax, dist_pt, dist_p_LB, dist_p_LB_1, dist_p_LB_2
        cdef ITYPE_t i, i1, i2, i_node

        cdef stack_item item

        # This will keep track of any indices with distances values.  If at
        # the end of the tree traversal, this index is in the last position,
        # then the warning flag will be set.
        cdef ITYPE_t check_index = -1

        item.i_node = 0
        item.dist_p_LB = calc_dist_p_LB(pt, node_centroid_arr,
                                        node_info.radius,
                                        n_features, p)
        stack_push(node_stack, item)

        # create pointers to the priority-queue/max-heap functions.
        # they both can operate on near_set_dist and near_set_indx
        cdef DTYPE_t (*heapqueue_largest)(DTYPE_t*, ITYPE_t)
        cdef ITYPE_t (*heapqueue_idx_largest)(ITYPE_t*, ITYPE_t)
        cdef void (*heapqueue_insert)(DTYPE_t, ITYPE_t, DTYPE_t*,
                                      ITYPE_t*, ITYPE_t)

        if use_max_heap:
            heapqueue_largest = &max_heap_largest
            heapqueue_idx_largest = &max_heap_idx_largest
            heapqueue_insert = &max_heap_insert
        else:
            heapqueue_largest = &pqueue_largest
            heapqueue_idx_largest = &pqueue_idx_largest
            heapqueue_insert = &pqueue_insert

        while(node_stack.n > 0):
            item = stack_pop(node_stack)
            i_node = item.i_node
            dist_p_LB = item.dist_p_LB

            node_info = node_info_arr + i_node

            #------------------------------------------------------------
            # Case 0: query point is exactly on the boundary.  Set
            #         warning flag
            if dist_p_LB == heapqueue_largest(near_set_dist, k):
                # store index of point with same distance:
                # we'll check it later
                check_index = heapqueue_idx_largest(near_set_indx, k)
                continue

            #------------------------------------------------------------
            # Case 1: query point is outside node radius
            elif dist_p_LB > heapqueue_largest(near_set_dist, k):
                continue

            #------------------------------------------------------------
            # Case 2: this is a leaf node.  Update set of nearby points
            elif node_info.is_leaf:
                for i from node_info.idx_start <= i < node_info.idx_end:
                    dist_pt = dist_p(pt,
                                     data + n_features * idx_array[i],
                                     n_features, p)

                    dmax = heapqueue_largest(near_set_dist, k)

                    if dist_pt == dmax:
                        check_index = heapqueue_idx_largest(near_set_indx, k)

                    elif dist_pt < dmax:
                        heapqueue_insert(dist_pt, idx_array[i],
                                         near_set_dist, near_set_indx, k)
                        if dmax == heapqueue_largest(near_set_dist, k):
                            check_index = heapqueue_idx_largest(near_set_indx,
                                                                k)

            #------------------------------------------------------------
            # Case 3: Node is not a leaf.  Recursively query subnodes
            #         starting with the one whose centroid is closest
            else:
                i1 = 2 * i_node + 1
                i2 = i1 + 1
                dist_p_LB_1 = calc_dist_p_LB(pt, (node_centroid_arr
                                                  + i1 * n_features),
                                             node_info_arr[i1].radius,
                                             n_features, p)
                dist_p_LB_2 = calc_dist_p_LB(pt, (node_centroid_arr
                                                  + i2 * n_features),
                                             node_info_arr[i2].radius,
                                             n_features, p)

                # append children to stack: last-in-first-out
                if dist_p_LB_2 <= dist_p_LB_1:
                    item.i_node = i1
                    item.dist_p_LB = dist_p_LB_1
                    stack_push(node_stack, item)

                    item.i_node = i2
                    item.dist_p_LB = dist_p_LB_2
                    stack_push(node_stack, item)

                else:
                    item.i_node = i2
                    item.dist_p_LB = dist_p_LB_2
                    stack_push(node_stack, item)

                    item.i_node = i1
                    item.dist_p_LB = dist_p_LB_1
                    stack_push(node_stack, item)

        if check_index == heapqueue_idx_largest(near_set_indx, k):
            self.warning_flag = True

        for i from 0 <= i < k:
            near_set_dist[i] = dist_from_dist_p(near_set_dist[i], p)

    cdef ITYPE_t query_radius_count_(BallTree self,
                                     DTYPE_t* pt, DTYPE_t r,
                                     stack* node_stack):
        cdef DTYPE_t* data = <DTYPE_t*> self.data.data
        cdef ITYPE_t* idx_array = <ITYPE_t*> self.idx_array.data
        cdef DTYPE_t* node_centroid_arr = <DTYPE_t*>self.node_centroid_arr.data
        cdef NodeInfo* node_info_arr = <NodeInfo*> self.node_info_arr.data
        cdef NodeInfo* node_info = node_info_arr

        cdef DTYPE_t p = self.p
        cdef ITYPE_t n_features = self.data.shape[1]
        cdef ITYPE_t i, i_node
        cdef ITYPE_t count = 0
        cdef DTYPE_t r_p = dist_p_from_dist(r, p)
        cdef DTYPE_t dist_pt

        cdef stack_item item

        item.i_node = 0
        stack_push(node_stack, item)

        while(node_stack.n > 0):
            item = stack_pop(node_stack)
            i_node = item.i_node
            node_info = node_info_arr + i_node

            dist_pt = dist(pt, node_centroid_arr + n_features * i_node,
                           n_features, p)

            #------------------------------------------------------------
            # Case 1: all node points are outside distance r.
            #         prune this branch.
            if dist_pt - node_info.radius > r:
                continue

            #------------------------------------------------------------
            # Case 2: all node points are within distance r
            #         add all points
            elif dist_pt + node_info.radius < r:
                count += (node_info.idx_end - node_info.idx_start)

            #------------------------------------------------------------
            # Case 3: this is a leaf node.  Go through all points to
            #         determine if they fall within radius
            elif node_info.is_leaf:
                for i from node_info.idx_start <= i < node_info.idx_end:
                    dist_pt = dist_p(pt,
                                     data + idx_array[i] * n_features,
                                     n_features, p)
                    if dist_pt <= r_p:
                        count += 1

            #------------------------------------------------------------
            # Case 4: Node is not a leaf.  Recursively query subnodes
            else:
                item.i_node = 2 * i_node + 1
                stack_push(node_stack, item)

                item.i_node = i = 2 * i_node + 2
                stack_push(node_stack, item)

        return count

    cdef ITYPE_t query_radius_idx_only_(BallTree self,
                                        DTYPE_t* pt, DTYPE_t r,
                                        ITYPE_t* indices,
                                        stack* node_stack):
        cdef DTYPE_t* data = <DTYPE_t*> self.data.data
        cdef ITYPE_t* idx_array = <ITYPE_t*> self.idx_array.data
        cdef DTYPE_t* node_centroid_arr = <DTYPE_t*>self.node_centroid_arr.data
        cdef NodeInfo* node_info_arr = <NodeInfo*> self.node_info_arr.data
        cdef NodeInfo* node_info = node_info_arr

        cdef DTYPE_t p = self.p
        cdef ITYPE_t n_features = self.data.shape[1]
        cdef ITYPE_t i, i_node
        cdef ITYPE_t idx_i = 0
        cdef DTYPE_t r_p = dist_p_from_dist(r, p)
        cdef DTYPE_t dist_pt

        cdef stack_item item

        item.i_node = 0
        stack_push(node_stack, item)

        while(node_stack.n > 0):
            item = stack_pop(node_stack)
            i_node = item.i_node
            node_info = node_info_arr + i_node

            dist_pt = dist(pt, node_centroid_arr + n_features * i_node,
                           n_features, p)

            #------------------------------------------------------------
            # Case 1: all node points are outside distance r.
            #         prune this branch.
            if dist_pt - node_info.radius > r:
                continue

            #------------------------------------------------------------
            # Case 2: all node points are within distance r
            #         add all points
            elif dist_pt + node_info.radius < r:
                for i from node_info.idx_start <= i < node_info.idx_end:
                    indices[idx_i] = idx_array[i]
                    idx_i += 1

            #------------------------------------------------------------
            # Case 3: this is a leaf node.  Go through all points to
            #         determine if they fall within radius
            elif node_info.is_leaf:
                for i from node_info.idx_start <= i < node_info.idx_end:
                    dist_pt = dist_p(pt,
                                     data + idx_array[i] * n_features,
                                     n_features, p)
                    if dist_pt <= r_p:
                        indices[idx_i] = idx_array[i]
                        idx_i += 1

            #------------------------------------------------------------
            # Case 4: Node is not a leaf.  Recursively query subnodes
            else:
                item.i_node = 2 * i_node + 1
                stack_push(node_stack, item)

                item.i_node = i = 2 * i_node + 2
                stack_push(node_stack, item)

        return idx_i

    cdef ITYPE_t query_radius_distances_(BallTree self,
                                         DTYPE_t* pt, DTYPE_t r,
                                         ITYPE_t* indices,
                                         DTYPE_t* distances,
                                         stack* node_stack):
        cdef DTYPE_t* data = <DTYPE_t*> self.data.data
        cdef ITYPE_t* idx_array = <ITYPE_t*> self.idx_array.data
        cdef DTYPE_t* node_centroid_arr = <DTYPE_t*>self.node_centroid_arr.data
        cdef NodeInfo* node_info_arr = <NodeInfo*> self.node_info_arr.data
        cdef NodeInfo* node_info = node_info_arr

        cdef DTYPE_t p = self.p
        cdef ITYPE_t n_features = self.data.shape[1]
        cdef ITYPE_t i, i_node
        cdef ITYPE_t idx_i = 0
        cdef DTYPE_t r_p = dist_p_from_dist(r, p)
        cdef DTYPE_t dist_pt

        cdef stack_item item

        item.i_node = 0
        stack_push(node_stack, item)

        while(node_stack.n > 0):
            item = stack_pop(node_stack)
            i_node = item.i_node
            node_info = node_info_arr + i_node

            dist_pt = dist(pt, node_centroid_arr + n_features * i_node,
                           n_features, p)

            #------------------------------------------------------------
            # Case 1: all node points are outside distance r.
            #         prune this branch.
            if dist_pt - node_info.radius > r:
                continue

            #------------------------------------------------------------
            # Case 2: all node points are within distance r
            #         add all points
            elif dist_pt + node_info.radius < r:
                for i from node_info.idx_start <= i < node_info.idx_end:
                    dist_pt = dist(pt,
                                   data + idx_array[i] * n_features,
                                   n_features, p)
                    indices[idx_i] = idx_array[i]
                    distances[idx_i] = dist_pt
                    idx_i += 1

            #------------------------------------------------------------
            # Case 3: this is a leaf node.  Go through all points to
            #         determine if they fall within radius
            elif node_info.is_leaf:
                for i from node_info.idx_start <= i < node_info.idx_end:
                    dist_pt = dist_p(pt,
                                     data + idx_array[i] * n_features,
                                     n_features, p)
                    if dist_pt <= r_p:
                        indices[idx_i] = idx_array[i]
                        distances[idx_i] = dist_from_dist_p(dist_pt, p)
                        idx_i += 1

            #------------------------------------------------------------
            # Case 4: Node is not a leaf.  Recursively query subnodes
            else:
                item.i_node = 2 * i_node + 1
                stack_push(node_stack, item)

                item.i_node = i = 2 * i_node + 2
                stack_push(node_stack, item)

        return idx_i


######################################################################
# Helper functions for building and querying
#
@cython.profile(False)
cdef inline void copy_array(DTYPE_t* x, DTYPE_t* y, ITYPE_t n):
    # copy array y into array x
    cdef ITYPE_t i
    for i from 0 <= i < n:
        x[i] = y[i]


@cython.cdivision(True)
cdef void compute_centroid(DTYPE_t* centroid,
                           DTYPE_t* data,
                           ITYPE_t* node_indices,
                           ITYPE_t n_features,
                           ITYPE_t n_points):
    # `centroid` points to an array of length n_features
    # `data` points to an array of length n_samples * n_features
    # `node_indices` = idx_array + idx_start
    cdef DTYPE_t *this_pt
    cdef ITYPE_t i, j

    for j from 0 <= j < n_features:
        centroid[j] = 0

    for i from 0 <= i < n_points:
        this_pt = data + n_features * node_indices[i]
        for j from 0 <= j < n_features:
            centroid[j] += this_pt[j]

    for j from 0 <= j < n_features:
        centroid[j] /= n_points


cdef ITYPE_t find_split_dim(DTYPE_t* data,
                            ITYPE_t* node_indices,
                            ITYPE_t n_features,
                            ITYPE_t n_points):
    # this computes the following
    # j_max = np.argmax(np.max(data, 0) - np.min(data, 0))
    cdef DTYPE_t min_val, max_val, val, spread, max_spread
    cdef ITYPE_t i, j, j_max

    j_max = 0
    max_spread = 0

    for j from 0 <= j < n_features:
        max_val = data[node_indices[0] * n_features + j]
        min_val = max_val
        for i from 1 <= i < n_points:
            val = data[node_indices[i] * n_features + j]
            max_val = dmax(max_val, val)
            min_val = dmin(min_val, val)
        spread = max_val - min_val
        if spread > max_spread:
            max_spread = spread
            j_max = j
    return j_max


@cython.profile(False)
cdef inline void iswap(ITYPE_t* arr, ITYPE_t i1, ITYPE_t i2):
    cdef ITYPE_t tmp = arr[i1]
    arr[i1] = arr[i2]
    arr[i2] = tmp


@cython.profile(False)
cdef inline void dswap(DTYPE_t* arr, ITYPE_t i1, ITYPE_t i2):
    cdef DTYPE_t tmp = arr[i1]
    arr[i1] = arr[i2]
    arr[i2] = tmp


cdef void partition_indices(DTYPE_t* data,
                            ITYPE_t* node_indices,
                            ITYPE_t split_dim,
                            ITYPE_t split_index,
                            ITYPE_t n_features,
                            ITYPE_t n_points):
    # partition_indices will modify the array node_indices between
    # indices 0 and n_points.  Upon return (assuming numpy-style slicing)
    #   data[node_indices[0:split_index], split_dim]
    #     <= data[node_indices[split_index], split_dim]
    # and
    #   data[node_indices[split_index], split_dim]
    #     <= data[node_indices[split_index:n_points], split_dim]
    # will hold.  The algorithm amounts to a partial quicksort
    cdef ITYPE_t left, right, midindex, i
    cdef DTYPE_t d1, d2
    left = 0
    right = n_points - 1

    while True:
        midindex = left
        for i from left <= i < right:
            d1 = data[node_indices[i] * n_features + split_dim]
            d2 = data[node_indices[right] * n_features + split_dim]
            if d1 < d2:
                iswap(node_indices, i, midindex)
                midindex += 1
        iswap(node_indices, midindex, right)
        if midindex == split_index:
            break
        elif midindex < split_index:
            left = midindex + 1
        else:
            right = midindex - 1


######################################################################
# calc_dist_LB
# calc_dist_p_LB
#  This calculates the lower-bound distance between a point and a node
@cython.profile(False)
cdef inline DTYPE_t calc_dist_LB(DTYPE_t* pt,
                                 DTYPE_t* centroid,
                                 DTYPE_t radius,
                                 ITYPE_t n_features,
                                 DTYPE_t p):
    return dmax(0, (dist(pt, centroid, n_features, p)
                    - radius))


@cython.profile(False)
cdef inline DTYPE_t calc_dist_p_LB(DTYPE_t* pt,
                                   DTYPE_t* centroid,
                                   DTYPE_t radius,
                                   ITYPE_t n_features,
                                   DTYPE_t p):
    return dist_p_from_dist(dmax(0, (dist(pt, centroid, n_features, p)
                                     - radius)), p)


######################################################################
# priority queue
#  This is used to keep track of the neighbors as they are found.
#  It keeps the list of neighbors sorted, and inserts each new item
#  into the list.  In this fixed-size implementation, empty elements
#  are represented by infinities.
@cython.profile(False)
cdef inline DTYPE_t pqueue_largest(DTYPE_t* queue, ITYPE_t queue_size):
    return queue[queue_size - 1]


cdef inline ITYPE_t pqueue_idx_largest(ITYPE_t* idx_array, ITYPE_t queue_size):
    return idx_array[queue_size - 1]


cdef inline void pqueue_insert(DTYPE_t val, ITYPE_t i_val,
                               DTYPE_t* queue, ITYPE_t* idx_array,
                               ITYPE_t queue_size):
    cdef ITYPE_t i_lower = 0
    cdef ITYPE_t i_upper = queue_size - 1
    cdef ITYPE_t i_mid
    cdef ITYPE_t i

    if val >= queue[i_upper]:
        return
    elif val <= queue[i_lower]:
        i_mid = i_lower
    else:
        while True:
            if (i_upper - i_lower) < 2:
                i_mid = i_lower + 1
                break
            else:
                i_mid = (i_lower + i_upper) / 2

            if i_mid == i_lower:
                i_mid += 1
                break

            if val >= queue[i_mid]:
                i_lower = i_mid
            else:
                i_upper = i_mid

    for i from queue_size > i > i_mid:
        queue[i] = queue[i - 1]
        idx_array[i] = idx_array[i - 1]

    queue[i_mid] = val
    idx_array[i_mid] = i_val


######################################################################
# max_heap
#
#  This is a basic implementation of a fixed-size binary max-heap.
#  It can be used in place of priority_queue to keep track of the
#  k-nearest neighbors in a query.  The implementation is faster than
#  priority_queue for a very large number of neighbors (k > 50 or so).
#  The implementation is slower than priority_queue for fewer neighbors.
#  The other disadvantage is that for max_heap, the indices/distances must
#  be sorted upon completion of the query.  In priority_queue, the indices
#  and distances are sorted without an extra call.
#
#  The root node is at heap[0].  The two child nodes of node i are at
#  (2 * i + 1) and (2 * i + 2).
#  The parent node of node i is node floor((i-1)/2).  Node 0 has no parent.
#  A max heap has (heap[i] >= heap[2 * i + 1]) and (heap[i] >= heap[2 * i + 2])
#  for all valid indices.
#
#  In this implementation, an empty heap should be full of infinities
#
#  As part of this implementation, there is a quicksort provided with
#  `sort_dist_idx()`
@cython.profile(False)
cdef inline DTYPE_t max_heap_largest(DTYPE_t* heap, ITYPE_t k):
    return heap[0]


@cython.profile(False)
cdef inline ITYPE_t max_heap_idx_largest(ITYPE_t* idx_array, ITYPE_t k):
    return idx_array[0]


cdef void max_heap_insert(DTYPE_t val, ITYPE_t i_val,
                          DTYPE_t* heap,
                          ITYPE_t* idx_array,
                          ITYPE_t heap_size):
    cdef ITYPE_t i, ic1, ic2, i_tmp
    cdef DTYPE_t d_tmp

    # check if val should be in heap
    if val > heap[0]:
        return

    # insert val at position zero
    heap[0] = val
    idx_array[0] = i_val

    #descend the heap, swapping values until the max heap criterion is met
    i = 0
    while 1:
        ic1 = 2 * i + 1
        ic2 = ic1 + 1

        if ic1 >= heap_size:
            break
        elif ic2 >= heap_size:
            if heap[ic1] > val:
                i_swap = ic1
            else:
                break
        elif heap[ic1] >= heap[ic2]:
            if val < heap[ic1]:
                i_swap = ic1
            else:
                break
        else:
            if val < heap[ic2]:
                i_swap = ic2
            else:
                break

        heap[i] = heap[i_swap]
        idx_array[i] = idx_array[i_swap]

        i = i_swap

    heap[i] = val
    idx_array[i] = i_val


######################################################################
# sort_dist_idx :
#  this is a quicksort implementation which sorts `dist` and
#  simultaneously performs the same swaps on `idx`.
cdef void sort_dist_idx(DTYPE_t* dist, ITYPE_t* idx, ITYPE_t k):
    cdef ITYPE_t pivot_idx
    if k > 1:
        pivot_idx = partition_dist_idx(dist, idx, k)

        sort_dist_idx(dist, idx, pivot_idx)

        sort_dist_idx(dist + pivot_idx + 1,
                      idx + pivot_idx + 1,
                      k - pivot_idx - 1)


cdef ITYPE_t partition_dist_idx(DTYPE_t* dist, ITYPE_t* idx, ITYPE_t k):
    cdef ITYPE_t pivot_idx = k / 2
    cdef DTYPE_t pivot_val = dist[pivot_idx]
    cdef ITYPE_t store_idx = 0
    cdef ITYPE_t i

    dswap(dist, pivot_idx, k - 1)
    iswap(idx, pivot_idx, k - 1)

    for i from 0 <= i < k - 1:
        if dist[i] < pivot_val:
            dswap(dist, i, store_idx)
            iswap(idx, i, store_idx)
            store_idx += 1
    dswap(dist, store_idx, k - 1)
    iswap(idx, store_idx, k - 1)
    return store_idx
