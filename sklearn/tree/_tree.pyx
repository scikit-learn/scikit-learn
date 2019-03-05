# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
# cython: language_level=3

# Authors: Gilles Louppe <g.louppe@gmail.com>
#          Peter Prettenhofer <peter.prettenhofer@gmail.com>
#          Brian Holt <bdholt1@gmail.com>
#          Noel Dawe <noel@dawe.me>
#          Satrajit Gosh <satrajit.ghosh@gmail.com>
#          Lars Buitinck
#          Arnaud Joly <arnaud.v.joly@gmail.com>
#          Joel Nothman <joel.nothman@gmail.com>
#          Fares Hedayati <fares.hedayati@gmail.com>
#          Jacob Schreiber <jmschreiber91@gmail.com>
#          Nelson Liu <nelson@nelsonliu.me>
#
# License: BSD 3 clause

from cpython cimport Py_INCREF, PyObject

from libc.stdlib cimport free
from libc.math cimport fabs
from libc.string cimport memcpy
from libc.string cimport memset

import numpy as np
cimport numpy as np
np.import_array()

from scipy.sparse import issparse
from scipy.sparse import csc_matrix
from scipy.sparse import csr_matrix

from ._utils cimport Stack
from ._utils cimport StackRecord
from ._utils cimport PriorityHeap
from ._utils cimport PriorityHeapRecord
from ._utils cimport safe_realloc
from ._utils cimport sizet_ptr_to_ndarray

cdef extern from "numpy/arrayobject.h":
    object PyArray_NewFromDescr(object subtype, np.dtype descr,
                                int nd, np.npy_intp* dims,
                                np.npy_intp* strides,
                                void* data, int flags, object obj)

# =============================================================================
# Types and constants
# =============================================================================

from numpy import float32 as DTYPE
from numpy import float64 as DOUBLE

cdef double INFINITY = np.inf
cdef double EPSILON = np.finfo('double').eps

# Some handy constants (BestFirstTreeBuilder)
cdef int IS_FIRST = 1
cdef int IS_NOT_FIRST = 0
cdef int IS_LEFT = 1
cdef int IS_NOT_LEFT = 0

TREE_LEAF = -1
TREE_UNDEFINED = -2
cdef SIZE_t _TREE_LEAF = TREE_LEAF
cdef SIZE_t _TREE_UNDEFINED = TREE_UNDEFINED
cdef SIZE_t INITIAL_STACK_SIZE = 10

# Repeat struct definition for numpy
NODE_DTYPE = np.dtype({
    'names': ['left_child', 'right_child', 'feature', 'threshold', 'impurity',
              'n_node_samples', 'weighted_n_node_samples'],
    'formats': [np.intp, np.intp, np.intp, np.float64, np.float64, np.intp,
                np.float64],
    'offsets': [
        <Py_ssize_t> &(<Node*> NULL).left_child,
        <Py_ssize_t> &(<Node*> NULL).right_child,
        <Py_ssize_t> &(<Node*> NULL).feature,
        <Py_ssize_t> &(<Node*> NULL).threshold,
        <Py_ssize_t> &(<Node*> NULL).impurity,
        <Py_ssize_t> &(<Node*> NULL).n_node_samples,
        <Py_ssize_t> &(<Node*> NULL).weighted_n_node_samples
    ]
})

# =============================================================================
# TreeBuilder
# =============================================================================

cdef class TreeBuilder:
    """Interface for different tree building strategies."""

    cpdef build(self, Tree tree, object X, np.ndarray y,
                np.ndarray sample_weight=None,
                np.ndarray X_idx_sorted=None):
        """Build a decision tree from the training set (X, y)."""
        pass

    cdef inline _check_input(self, object X, np.ndarray y,
                             np.ndarray sample_weight):
        """Check input dtype, layout and format"""
        if issparse(X):
            X = X.tocsc()
            X.sort_indices()

            if X.data.dtype != DTYPE:
                X.data = np.ascontiguousarray(X.data, dtype=DTYPE)

            if X.indices.dtype != np.int32 or X.indptr.dtype != np.int32:
                raise ValueError("No support for np.int64 index based "
                                 "sparse matrices")

        elif X.dtype != DTYPE:
            # since we have to copy we will make it fortran for efficiency
            X = np.asfortranarray(X, dtype=DTYPE)

        if y.dtype != DOUBLE or not y.flags.contiguous:
            y = np.ascontiguousarray(y, dtype=DOUBLE)

        if (sample_weight is not None and
            (sample_weight.dtype != DOUBLE or
            not sample_weight.flags.contiguous)):
                sample_weight = np.asarray(sample_weight, dtype=DOUBLE,
                                           order="C")

        return X, y, sample_weight

# Depth first builder ---------------------------------------------------------

cdef class DepthFirstTreeBuilder(TreeBuilder):
    """Build a decision tree in depth-first fashion."""

    def __cinit__(self, Splitter splitter, SIZE_t min_samples_split,
                  SIZE_t min_samples_leaf, double min_weight_leaf,
                  SIZE_t max_depth, double min_impurity_decrease,
                  double min_impurity_split):
        self.splitter = splitter
        self.min_samples_split = min_samples_split
        self.min_samples_leaf = min_samples_leaf
        self.min_weight_leaf = min_weight_leaf
        self.max_depth = max_depth
        self.min_impurity_decrease = min_impurity_decrease
        self.min_impurity_split = min_impurity_split

    cpdef build(self, Tree tree, object X, np.ndarray y,
                np.ndarray sample_weight=None,
                np.ndarray X_idx_sorted=None):
        """Build a decision tree from the training set (X, y)."""

        # check input
        X, y, sample_weight = self._check_input(X, y, sample_weight)

        cdef DOUBLE_t* sample_weight_ptr = NULL
        if sample_weight is not None:
            sample_weight_ptr = <DOUBLE_t*> sample_weight.data

        # Initial capacity
        cdef int init_capacity

        if tree.max_depth <= 10:
            init_capacity = (2 ** (tree.max_depth + 1)) - 1
        else:
            init_capacity = 2047

        tree._resize(init_capacity)

        # Parameters
        cdef Splitter splitter = self.splitter
        cdef SIZE_t max_depth = self.max_depth
        cdef SIZE_t min_samples_leaf = self.min_samples_leaf
        cdef double min_weight_leaf = self.min_weight_leaf
        cdef SIZE_t min_samples_split = self.min_samples_split
        cdef double min_impurity_decrease = self.min_impurity_decrease
        cdef double min_impurity_split = self.min_impurity_split

        # Recursive partition (without actual recursion)
        splitter.init(X, y, sample_weight_ptr, X_idx_sorted)

        cdef SIZE_t start
        cdef SIZE_t end
        cdef SIZE_t depth
        cdef SIZE_t parent
        cdef bint is_left
        cdef SIZE_t n_node_samples = splitter.n_samples
        cdef double weighted_n_samples = splitter.weighted_n_samples
        cdef double weighted_n_node_samples
        cdef SplitRecord split
        cdef SIZE_t node_id

        cdef double impurity = INFINITY
        cdef SIZE_t n_constant_features
        cdef bint is_leaf
        cdef bint first = 1
        cdef SIZE_t max_depth_seen = -1
        cdef int rc = 0

        cdef Stack stack = Stack(INITIAL_STACK_SIZE)
        cdef StackRecord stack_record

        with nogil:
            # push root node onto stack
            rc = stack.push(0, n_node_samples, 0, _TREE_UNDEFINED, 0, INFINITY, 0)
            if rc == -1:
                # got return code -1 - out-of-memory
                with gil:
                    raise MemoryError()

            while not stack.is_empty():
                stack.pop(&stack_record)

                start = stack_record.start
                end = stack_record.end
                depth = stack_record.depth
                parent = stack_record.parent
                is_left = stack_record.is_left
                impurity = stack_record.impurity
                n_constant_features = stack_record.n_constant_features

                n_node_samples = end - start
                splitter.node_reset(start, end, &weighted_n_node_samples)

                is_leaf = (depth >= max_depth or
                           n_node_samples < min_samples_split or
                           n_node_samples < 2 * min_samples_leaf or
                           weighted_n_node_samples < 2 * min_weight_leaf)

                if first:
                    impurity = splitter.node_impurity()
                    first = 0

                is_leaf = (is_leaf or
                           (impurity <= min_impurity_split))

                if not is_leaf:
                    splitter.node_split(impurity, &split, &n_constant_features)
                    # If EPSILON=0 in the below comparison, float precision
                    # issues stop splitting, producing trees that are
                    # dissimilar to v0.18
                    is_leaf = (is_leaf or split.pos >= end or
                               (split.improvement + EPSILON <
                                min_impurity_decrease))

                node_id = tree._add_node(parent, is_left, is_leaf, split.feature,
                                         split.threshold, impurity, n_node_samples,
                                         weighted_n_node_samples)

                if node_id == <SIZE_t>(-1):
                    rc = -1
                    break

                # Store value for all nodes, to facilitate tree/model
                # inspection and interpretation
                splitter.node_value(tree.value + node_id * tree.value_stride)

                if not is_leaf:
                    # Push right child on stack
                    rc = stack.push(split.pos, end, depth + 1, node_id, 0,
                                    split.impurity_right, n_constant_features)
                    if rc == -1:
                        break

                    # Push left child on stack
                    rc = stack.push(start, split.pos, depth + 1, node_id, 1,
                                    split.impurity_left, n_constant_features)
                    if rc == -1:
                        break

                if depth > max_depth_seen:
                    max_depth_seen = depth

            if rc >= 0:
                rc = tree._resize_c(tree.node_count)

            if rc >= 0:
                tree.max_depth = max_depth_seen
        if rc == -1:
            raise MemoryError()


# Best first builder ----------------------------------------------------------

cdef inline int _add_to_frontier(PriorityHeapRecord* rec,
                                 PriorityHeap frontier) nogil except -1:
    """Adds record ``rec`` to the priority queue ``frontier``

    Returns -1 in case of failure to allocate memory (and raise MemoryError)
    or 0 otherwise.
    """
    return frontier.push(rec.node_id, rec.start, rec.end, rec.pos, rec.depth,
                         rec.is_leaf, rec.improvement, rec.impurity,
                         rec.impurity_left, rec.impurity_right)


cdef class BestFirstTreeBuilder(TreeBuilder):
    """Build a decision tree in best-first fashion.

    The best node to expand is given by the node at the frontier that has the
    highest impurity improvement.
    """
    cdef SIZE_t max_leaf_nodes

    def __cinit__(self, Splitter splitter, SIZE_t min_samples_split,
                  SIZE_t min_samples_leaf,  min_weight_leaf,
                  SIZE_t max_depth, SIZE_t max_leaf_nodes,
                  double min_impurity_decrease, double min_impurity_split):
        self.splitter = splitter
        self.min_samples_split = min_samples_split
        self.min_samples_leaf = min_samples_leaf
        self.min_weight_leaf = min_weight_leaf
        self.max_depth = max_depth
        self.max_leaf_nodes = max_leaf_nodes
        self.min_impurity_decrease = min_impurity_decrease
        self.min_impurity_split = min_impurity_split

    cpdef build(self, Tree tree, object X, np.ndarray y,
                np.ndarray sample_weight=None,
                np.ndarray X_idx_sorted=None):
        """Build a decision tree from the training set (X, y)."""

        # check input
        X, y, sample_weight = self._check_input(X, y, sample_weight)

        cdef DOUBLE_t* sample_weight_ptr = NULL
        if sample_weight is not None:
            sample_weight_ptr = <DOUBLE_t*> sample_weight.data

        # Parameters
        cdef Splitter splitter = self.splitter
        cdef SIZE_t max_leaf_nodes = self.max_leaf_nodes
        cdef SIZE_t min_samples_leaf = self.min_samples_leaf
        cdef double min_weight_leaf = self.min_weight_leaf
        cdef SIZE_t min_samples_split = self.min_samples_split

        # Recursive partition (without actual recursion)
        splitter.init(X, y, sample_weight_ptr, X_idx_sorted)

        cdef PriorityHeap frontier = PriorityHeap(INITIAL_STACK_SIZE)
        cdef PriorityHeapRecord record
        cdef PriorityHeapRecord split_node_left
        cdef PriorityHeapRecord split_node_right

        cdef SIZE_t n_node_samples = splitter.n_samples
        cdef SIZE_t max_split_nodes = max_leaf_nodes - 1
        cdef bint is_leaf
        cdef SIZE_t max_depth_seen = -1
        cdef int rc = 0
        cdef Node* node

        # Initial capacity
        cdef SIZE_t init_capacity = max_split_nodes + max_leaf_nodes
        tree._resize(init_capacity)

        with nogil:
            # add root to frontier
            rc = self._add_split_node(splitter, tree, 0, n_node_samples,
                                      INFINITY, IS_FIRST, IS_LEFT, NULL, 0,
                                      &split_node_left)
            if rc >= 0:
                rc = _add_to_frontier(&split_node_left, frontier)

            if rc == -1:
                with gil:
                    raise MemoryError()

            while not frontier.is_empty():
                frontier.pop(&record)

                node = &tree.nodes[record.node_id]
                is_leaf = (record.is_leaf or max_split_nodes <= 0)

                if is_leaf:
                    # Node is not expandable; set node as leaf
                    node.left_child = _TREE_LEAF
                    node.right_child = _TREE_LEAF
                    node.feature = _TREE_UNDEFINED
                    node.threshold = _TREE_UNDEFINED

                else:
                    # Node is expandable

                    # Decrement number of split nodes available
                    max_split_nodes -= 1

                    # Compute left split node
                    rc = self._add_split_node(splitter, tree,
                                              record.start, record.pos,
                                              record.impurity_left,
                                              IS_NOT_FIRST, IS_LEFT, node,
                                              record.depth + 1,
                                              &split_node_left)
                    if rc == -1:
                        break

                    # tree.nodes may have changed
                    node = &tree.nodes[record.node_id]

                    # Compute right split node
                    rc = self._add_split_node(splitter, tree, record.pos,
                                              record.end,
                                              record.impurity_right,
                                              IS_NOT_FIRST, IS_NOT_LEFT, node,
                                              record.depth + 1,
                                              &split_node_right)
                    if rc == -1:
                        break

                    # Add nodes to queue
                    rc = _add_to_frontier(&split_node_left, frontier)
                    if rc == -1:
                        break

                    rc = _add_to_frontier(&split_node_right, frontier)
                    if rc == -1:
                        break

                if record.depth > max_depth_seen:
                    max_depth_seen = record.depth

            if rc >= 0:
                rc = tree._resize_c(tree.node_count)

            if rc >= 0:
                tree.max_depth = max_depth_seen

        if rc == -1:
            raise MemoryError()

    cdef inline int _add_split_node(self, Splitter splitter, Tree tree,
                                    SIZE_t start, SIZE_t end, double impurity,
                                    bint is_first, bint is_left, Node* parent,
                                    SIZE_t depth,
                                    PriorityHeapRecord* res) nogil except -1:
        """Adds node w/ partition ``[start, end)`` to the frontier. """
        cdef SplitRecord split
        cdef SIZE_t node_id
        cdef SIZE_t n_node_samples
        cdef SIZE_t n_constant_features = 0
        cdef double weighted_n_samples = splitter.weighted_n_samples
        cdef double min_impurity_decrease = self.min_impurity_decrease
        cdef double min_impurity_split = self.min_impurity_split
        cdef double weighted_n_node_samples
        cdef bint is_leaf
        cdef SIZE_t n_left, n_right
        cdef double imp_diff

        splitter.node_reset(start, end, &weighted_n_node_samples)

        if is_first:
            impurity = splitter.node_impurity()

        n_node_samples = end - start
        is_leaf = (depth >= self.max_depth or
                   n_node_samples < self.min_samples_split or
                   n_node_samples < 2 * self.min_samples_leaf or
                   weighted_n_node_samples < 2 * self.min_weight_leaf or
                   impurity <= min_impurity_split)

        if not is_leaf:
            splitter.node_split(impurity, &split, &n_constant_features)
            # If EPSILON=0 in the below comparison, float precision issues stop
            # splitting early, producing trees that are dissimilar to v0.18
            is_leaf = (is_leaf or split.pos >= end or
                       split.improvement + EPSILON < min_impurity_decrease)

        node_id = tree._add_node(parent - tree.nodes
                                 if parent != NULL
                                 else _TREE_UNDEFINED,
                                 is_left, is_leaf,
                                 split.feature, split.threshold, impurity, n_node_samples,
                                 weighted_n_node_samples)
        if node_id == <SIZE_t>(-1):
            return -1

        # compute values also for split nodes (might become leafs later).
        splitter.node_value(tree.value + node_id * tree.value_stride)

        res.node_id = node_id
        res.start = start
        res.end = end
        res.depth = depth
        res.impurity = impurity

        if not is_leaf:
            # is split node
            res.pos = split.pos
            res.is_leaf = 0
            res.improvement = split.improvement
            res.impurity_left = split.impurity_left
            res.impurity_right = split.impurity_right

        else:
            # is leaf => 0 improvement
            res.pos = end
            res.is_leaf = 1
            res.improvement = 0.0
            res.impurity_left = impurity
            res.impurity_right = impurity

        return 0


# =============================================================================
# Tree
# =============================================================================

cdef class Tree:
    """Array-based representation of a binary decision tree.

    The binary tree is represented as a number of parallel arrays. The i-th
    element of each array holds information about the node `i`. Node 0 is the
    tree's root. You can find a detailed description of all arrays in
    `_tree.pxd`. NOTE: Some of the arrays only apply to either leaves or split
    nodes, resp. In this case the values of nodes of the other type are
    arbitrary!

    Attributes
    ----------
    node_count : int
        The number of nodes (internal nodes + leaves) in the tree.

    capacity : int
        The current capacity (i.e., size) of the arrays, which is at least as
        great as `node_count`.

    max_depth : int
        The depth of the tree, i.e. the maximum depth of its leaves.

    children_left : array of int, shape [node_count]
        children_left[i] holds the node id of the left child of node i.
        For leaves, children_left[i] == TREE_LEAF. Otherwise,
        children_left[i] > i. This child handles the case where
        X[:, feature[i]] <= threshold[i].

    children_right : array of int, shape [node_count]
        children_right[i] holds the node id of the right child of node i.
        For leaves, children_right[i] == TREE_LEAF. Otherwise,
        children_right[i] > i. This child handles the case where
        X[:, feature[i]] > threshold[i].

    feature : array of int, shape [node_count]
        feature[i] holds the feature to split on, for the internal node i.

    threshold : array of double, shape [node_count]
        threshold[i] holds the threshold for the internal node i.

    value : array of double, shape [node_count, n_outputs, max_n_classes]
        Contains the constant prediction value of each node.

    impurity : array of double, shape [node_count]
        impurity[i] holds the impurity (i.e., the value of the splitting
        criterion) at node i.

    n_node_samples : array of int, shape [node_count]
        n_node_samples[i] holds the number of training samples reaching node i.

    weighted_n_node_samples : array of int, shape [node_count]
        weighted_n_node_samples[i] holds the weighted number of training samples
        reaching node i.
    """
    # Wrap for outside world.
    # WARNING: these reference the current `nodes` and `value` buffers, which
    # must not be freed by a subsequent memory allocation.
    # (i.e. through `_resize` or `__setstate__`)
    property n_classes:
        def __get__(self):
            return sizet_ptr_to_ndarray(self.n_classes, self.n_outputs)

    property children_left:
        def __get__(self):
            return self._get_node_ndarray()['left_child'][:self.node_count]

    property children_right:
        def __get__(self):
            return self._get_node_ndarray()['right_child'][:self.node_count]

    property n_leaves:
        def __get__(self):
            return np.sum(np.logical_and(
                self.children_left == -1,
                self.children_right == -1))

    property feature:
        def __get__(self):
            return self._get_node_ndarray()['feature'][:self.node_count]

    property threshold:
        def __get__(self):
            return self._get_node_ndarray()['threshold'][:self.node_count]

    property impurity:
        def __get__(self):
            return self._get_node_ndarray()['impurity'][:self.node_count]

    property n_node_samples:
        def __get__(self):
            return self._get_node_ndarray()['n_node_samples'][:self.node_count]

    property weighted_n_node_samples:
        def __get__(self):
            return self._get_node_ndarray()['weighted_n_node_samples'][:self.node_count]

    property value:
        def __get__(self):
            return self._get_value_ndarray()[:self.node_count]

    def __cinit__(self, int n_features, np.ndarray[SIZE_t, ndim=1] n_classes,
                  int n_outputs):
        """Constructor."""
        # Input/Output layout
        self.n_features = n_features
        self.n_outputs = n_outputs
        self.n_classes = NULL
        safe_realloc(&self.n_classes, n_outputs)

        self.max_n_classes = np.max(n_classes)
        self.value_stride = n_outputs * self.max_n_classes

        cdef SIZE_t k
        for k in range(n_outputs):
            self.n_classes[k] = n_classes[k]

        # Inner structures
        self.max_depth = 0
        self.node_count = 0
        self.capacity = 0
        self.value = NULL
        self.nodes = NULL

    def __dealloc__(self):
        """Destructor."""
        # Free all inner structures
        free(self.n_classes)
        free(self.value)
        free(self.nodes)

    def __reduce__(self):
        """Reduce re-implementation, for pickling."""
        return (Tree, (self.n_features,
                       sizet_ptr_to_ndarray(self.n_classes, self.n_outputs),
                       self.n_outputs), self.__getstate__())

    def __getstate__(self):
        """Getstate re-implementation, for pickling."""
        d = {}
        # capacity is inferred during the __setstate__ using nodes
        d["max_depth"] = self.max_depth
        d["node_count"] = self.node_count
        d["nodes"] = self._get_node_ndarray()
        d["values"] = self._get_value_ndarray()
        return d

    def __setstate__(self, d):
        """Setstate re-implementation, for unpickling."""
        self.max_depth = d["max_depth"]
        self.node_count = d["node_count"]

        if 'nodes' not in d:
            raise ValueError('You have loaded Tree version which '
                             'cannot be imported')

        node_ndarray = d['nodes']
        value_ndarray = d['values']

        value_shape = (node_ndarray.shape[0], self.n_outputs,
                       self.max_n_classes)
        if (node_ndarray.ndim != 1 or
                node_ndarray.dtype != NODE_DTYPE or
                not node_ndarray.flags.c_contiguous or
                value_ndarray.shape != value_shape or
                not value_ndarray.flags.c_contiguous or
                value_ndarray.dtype != np.float64):
            raise ValueError('Did not recognise loaded array layout')

        self.capacity = node_ndarray.shape[0]
        if self._resize_c(self.capacity) != 0:
            raise MemoryError("resizing tree to %d" % self.capacity)
        nodes = memcpy(self.nodes, (<np.ndarray> node_ndarray).data,
                       self.capacity * sizeof(Node))
        value = memcpy(self.value, (<np.ndarray> value_ndarray).data,
                       self.capacity * self.value_stride * sizeof(double))

    cdef int _resize(self, SIZE_t capacity) nogil except -1:
        """Resize all inner arrays to `capacity`, if `capacity` == -1, then
           double the size of the inner arrays.

        Returns -1 in case of failure to allocate memory (and raise MemoryError)
        or 0 otherwise.
        """
        if self._resize_c(capacity) != 0:
            # Acquire gil only if we need to raise
            with gil:
                raise MemoryError()

    # XXX using (size_t)(-1) is ugly, but SIZE_MAX is not available in C89
    # (i.e., older MSVC).
    cdef int _resize_c(self, SIZE_t capacity=<SIZE_t>(-1)) nogil except -1:
        """Guts of _resize

        Returns -1 in case of failure to allocate memory (and raise MemoryError)
        or 0 otherwise.
        """
        if capacity == self.capacity and self.nodes != NULL:
            return 0

        if capacity == <SIZE_t>(-1):
            if self.capacity == 0:
                capacity = 3  # default initial value
            else:
                capacity = 2 * self.capacity

        safe_realloc(&self.nodes, capacity)
        safe_realloc(&self.value, capacity * self.value_stride)

        # value memory is initialised to 0 to enable classifier argmax
        if capacity > self.capacity:
            memset(<void*>(self.value + self.capacity * self.value_stride), 0,
                   (capacity - self.capacity) * self.value_stride *
                   sizeof(double))

        # if capacity smaller than node_count, adjust the counter
        if capacity < self.node_count:
            self.node_count = capacity

        self.capacity = capacity
        return 0

    cdef SIZE_t _add_node(self, SIZE_t parent, bint is_left, bint is_leaf,
                          SIZE_t feature, double threshold, double impurity,
                          SIZE_t n_node_samples,
                          double weighted_n_node_samples) nogil except -1:
        """Add a node to the tree.

        The new node registers itself as the child of its parent.

        Returns (size_t)(-1) on error.
        """
        cdef SIZE_t node_id = self.node_count

        if node_id >= self.capacity:
            if self._resize_c() != 0:
                return <SIZE_t>(-1)

        cdef Node* node = &self.nodes[node_id]
        node.impurity = impurity
        node.n_node_samples = n_node_samples
        node.weighted_n_node_samples = weighted_n_node_samples

        if parent != _TREE_UNDEFINED:
            if is_left:
                self.nodes[parent].left_child = node_id
            else:
                self.nodes[parent].right_child = node_id

        if is_leaf:
            node.left_child = _TREE_LEAF
            node.right_child = _TREE_LEAF
            node.feature = _TREE_UNDEFINED
            node.threshold = _TREE_UNDEFINED

        else:
            # left_child and right_child will be set later
            node.feature = feature
            node.threshold = threshold

        self.node_count += 1

        return node_id

    cpdef np.ndarray predict(self, object X):
        """Predict target for X."""
        out = self._get_value_ndarray().take(self.apply(X), axis=0,
                                             mode='clip')
        if self.n_outputs == 1:
            out = out.reshape(X.shape[0], self.max_n_classes)
        return out

    cpdef np.ndarray apply(self, object X):
        """Finds the terminal region (=leaf node) for each sample in X."""
        if issparse(X):
            return self._apply_sparse_csr(X)
        else:
            return self._apply_dense(X)

    cdef inline np.ndarray _apply_dense(self, object X):
        """Finds the terminal region (=leaf node) for each sample in X."""

        # Check input
        if not isinstance(X, np.ndarray):
            raise ValueError("X should be in np.ndarray format, got %s"
                             % type(X))

        if X.dtype != DTYPE:
            raise ValueError("X.dtype should be np.float32, got %s" % X.dtype)

        # Extract input
        cdef DTYPE_t[:, :] X_ndarray = X
        cdef SIZE_t n_samples = X.shape[0]

        # Initialize output
        cdef np.ndarray[SIZE_t] out = np.zeros((n_samples,), dtype=np.intp)
        cdef SIZE_t* out_ptr = <SIZE_t*> out.data

        # Initialize auxiliary data-structure
        cdef Node* node = NULL
        cdef SIZE_t i = 0

        with nogil:
            for i in range(n_samples):
                node = self.nodes
                # While node not a leaf
                while node.left_child != _TREE_LEAF:
                    # ... and node.right_child != _TREE_LEAF:
                    if X_ndarray[i, node.feature] <= node.threshold:
                        node = &self.nodes[node.left_child]
                    else:
                        node = &self.nodes[node.right_child]

                out_ptr[i] = <SIZE_t>(node - self.nodes)  # node offset

        return out

    cdef inline np.ndarray _apply_sparse_csr(self, object X):
        """Finds the terminal region (=leaf node) for each sample in sparse X.
        """
        # Check input
        if not isinstance(X, csr_matrix):
            raise ValueError("X should be in csr_matrix format, got %s"
                             % type(X))

        if X.dtype != DTYPE:
            raise ValueError("X.dtype should be np.float32, got %s" % X.dtype)

        # Extract input
        cdef np.ndarray[ndim=1, dtype=DTYPE_t] X_data_ndarray = X.data
        cdef np.ndarray[ndim=1, dtype=INT32_t] X_indices_ndarray  = X.indices
        cdef np.ndarray[ndim=1, dtype=INT32_t] X_indptr_ndarray  = X.indptr

        cdef DTYPE_t* X_data = <DTYPE_t*>X_data_ndarray.data
        cdef INT32_t* X_indices = <INT32_t*>X_indices_ndarray.data
        cdef INT32_t* X_indptr = <INT32_t*>X_indptr_ndarray.data

        cdef SIZE_t n_samples = X.shape[0]
        cdef SIZE_t n_features = X.shape[1]

        # Initialize output
        cdef np.ndarray[SIZE_t, ndim=1] out = np.zeros((n_samples,),
                                                       dtype=np.intp)
        cdef SIZE_t* out_ptr = <SIZE_t*> out.data

        # Initialize auxiliary data-structure
        cdef DTYPE_t feature_value = 0.
        cdef Node* node = NULL
        cdef DTYPE_t* X_sample = NULL
        cdef SIZE_t i = 0
        cdef INT32_t k = 0

        # feature_to_sample as a data structure records the last seen sample
        # for each feature; functionally, it is an efficient way to identify
        # which features are nonzero in the present sample.
        cdef SIZE_t* feature_to_sample = NULL

        safe_realloc(&X_sample, n_features)
        safe_realloc(&feature_to_sample, n_features)

        with nogil:
            memset(feature_to_sample, -1, n_features * sizeof(SIZE_t))

            for i in range(n_samples):
                node = self.nodes

                for k in range(X_indptr[i], X_indptr[i + 1]):
                    feature_to_sample[X_indices[k]] = i
                    X_sample[X_indices[k]] = X_data[k]

                # While node not a leaf
                while node.left_child != _TREE_LEAF:
                    # ... and node.right_child != _TREE_LEAF:
                    if feature_to_sample[node.feature] == i:
                        feature_value = X_sample[node.feature]

                    else:
                        feature_value = 0.

                    if feature_value <= node.threshold:
                        node = &self.nodes[node.left_child]
                    else:
                        node = &self.nodes[node.right_child]

                out_ptr[i] = <SIZE_t>(node - self.nodes)  # node offset

            # Free auxiliary arrays
            free(X_sample)
            free(feature_to_sample)

        return out

    cpdef object decision_path(self, object X):
        """Finds the decision path (=node) for each sample in X."""
        if issparse(X):
            return self._decision_path_sparse_csr(X)
        else:
            return self._decision_path_dense(X)

    cdef inline object _decision_path_dense(self, object X):
        """Finds the decision path (=node) for each sample in X."""

        # Check input
        if not isinstance(X, np.ndarray):
            raise ValueError("X should be in np.ndarray format, got %s"
                             % type(X))

        if X.dtype != DTYPE:
            raise ValueError("X.dtype should be np.float32, got %s" % X.dtype)

        # Extract input
        cdef DTYPE_t[:, :] X_ndarray = X
        cdef SIZE_t n_samples = X.shape[0]

        # Initialize output
        cdef np.ndarray[SIZE_t] indptr = np.zeros(n_samples + 1, dtype=np.intp)
        cdef SIZE_t* indptr_ptr = <SIZE_t*> indptr.data

        cdef np.ndarray[SIZE_t] indices = np.zeros(n_samples *
                                                   (1 + self.max_depth),
                                                   dtype=np.intp)
        cdef SIZE_t* indices_ptr = <SIZE_t*> indices.data

        # Initialize auxiliary data-structure
        cdef Node* node = NULL
        cdef SIZE_t i = 0

        with nogil:
            for i in range(n_samples):
                node = self.nodes
                indptr_ptr[i + 1] = indptr_ptr[i]

                # Add all external nodes
                while node.left_child != _TREE_LEAF:
                    # ... and node.right_child != _TREE_LEAF:
                    indices_ptr[indptr_ptr[i + 1]] = <SIZE_t>(node - self.nodes)
                    indptr_ptr[i + 1] += 1

                    if X_ndarray[i, node.feature] <= node.threshold:
                        node = &self.nodes[node.left_child]
                    else:
                        node = &self.nodes[node.right_child]

                # Add the leave node
                indices_ptr[indptr_ptr[i + 1]] = <SIZE_t>(node - self.nodes)
                indptr_ptr[i + 1] += 1

        indices = indices[:indptr[n_samples]]
        cdef np.ndarray[SIZE_t] data = np.ones(shape=len(indices),
                                               dtype=np.intp)
        out = csr_matrix((data, indices, indptr),
                         shape=(n_samples, self.node_count))

        return out

    cdef inline object _decision_path_sparse_csr(self, object X):
        """Finds the decision path (=node) for each sample in X."""

        # Check input
        if not isinstance(X, csr_matrix):
            raise ValueError("X should be in csr_matrix format, got %s"
                             % type(X))

        if X.dtype != DTYPE:
            raise ValueError("X.dtype should be np.float32, got %s" % X.dtype)

        # Extract input
        cdef np.ndarray[ndim=1, dtype=DTYPE_t] X_data_ndarray = X.data
        cdef np.ndarray[ndim=1, dtype=INT32_t] X_indices_ndarray  = X.indices
        cdef np.ndarray[ndim=1, dtype=INT32_t] X_indptr_ndarray  = X.indptr

        cdef DTYPE_t* X_data = <DTYPE_t*>X_data_ndarray.data
        cdef INT32_t* X_indices = <INT32_t*>X_indices_ndarray.data
        cdef INT32_t* X_indptr = <INT32_t*>X_indptr_ndarray.data

        cdef SIZE_t n_samples = X.shape[0]
        cdef SIZE_t n_features = X.shape[1]

        # Initialize output
        cdef np.ndarray[SIZE_t] indptr = np.zeros(n_samples + 1, dtype=np.intp)
        cdef SIZE_t* indptr_ptr = <SIZE_t*> indptr.data

        cdef np.ndarray[SIZE_t] indices = np.zeros(n_samples *
                                                   (1 + self.max_depth),
                                                   dtype=np.intp)
        cdef SIZE_t* indices_ptr = <SIZE_t*> indices.data

        # Initialize auxiliary data-structure
        cdef DTYPE_t feature_value = 0.
        cdef Node* node = NULL
        cdef DTYPE_t* X_sample = NULL
        cdef SIZE_t i = 0
        cdef INT32_t k = 0

        # feature_to_sample as a data structure records the last seen sample
        # for each feature; functionally, it is an efficient way to identify
        # which features are nonzero in the present sample.
        cdef SIZE_t* feature_to_sample = NULL

        safe_realloc(&X_sample, n_features)
        safe_realloc(&feature_to_sample, n_features)

        with nogil:
            memset(feature_to_sample, -1, n_features * sizeof(SIZE_t))

            for i in range(n_samples):
                node = self.nodes
                indptr_ptr[i + 1] = indptr_ptr[i]

                for k in range(X_indptr[i], X_indptr[i + 1]):
                    feature_to_sample[X_indices[k]] = i
                    X_sample[X_indices[k]] = X_data[k]

                # While node not a leaf
                while node.left_child != _TREE_LEAF:
                    # ... and node.right_child != _TREE_LEAF:

                    indices_ptr[indptr_ptr[i + 1]] = <SIZE_t>(node - self.nodes)
                    indptr_ptr[i + 1] += 1

                    if feature_to_sample[node.feature] == i:
                        feature_value = X_sample[node.feature]

                    else:
                        feature_value = 0.

                    if feature_value <= node.threshold:
                        node = &self.nodes[node.left_child]
                    else:
                        node = &self.nodes[node.right_child]

                # Add the leave node
                indices_ptr[indptr_ptr[i + 1]] = <SIZE_t>(node - self.nodes)
                indptr_ptr[i + 1] += 1

            # Free auxiliary arrays
            free(X_sample)
            free(feature_to_sample)

        indices = indices[:indptr[n_samples]]
        cdef np.ndarray[SIZE_t] data = np.ones(shape=len(indices),
                                               dtype=np.intp)
        out = csr_matrix((data, indices, indptr),
                         shape=(n_samples, self.node_count))

        return out


    cpdef compute_feature_importances(self, normalize=True):
        """Computes the importance of each feature (aka variable)."""
        cdef Node* left
        cdef Node* right
        cdef Node* nodes = self.nodes
        cdef Node* node = nodes
        cdef Node* end_node = node + self.node_count

        cdef double normalizer = 0.

        cdef np.ndarray[np.float64_t, ndim=1] importances
        importances = np.zeros((self.n_features,))
        cdef DOUBLE_t* importance_data = <DOUBLE_t*>importances.data

        with nogil:
            while node != end_node:
                if node.left_child != _TREE_LEAF:
                    # ... and node.right_child != _TREE_LEAF:
                    left = &nodes[node.left_child]
                    right = &nodes[node.right_child]

                    importance_data[node.feature] += (
                        node.weighted_n_node_samples * node.impurity -
                        left.weighted_n_node_samples * left.impurity -
                        right.weighted_n_node_samples * right.impurity)
                node += 1

        importances /= nodes[0].weighted_n_node_samples

        if normalize:
            normalizer = np.sum(importances)

            if normalizer > 0.0:
                # Avoid dividing by zero (e.g., when root is pure)
                importances /= normalizer

        return importances

    cdef np.ndarray _get_value_ndarray(self):
        """Wraps value as a 3-d NumPy array.

        The array keeps a reference to this Tree, which manages the underlying
        memory.
        """
        cdef np.npy_intp shape[3]
        shape[0] = <np.npy_intp> self.node_count
        shape[1] = <np.npy_intp> self.n_outputs
        shape[2] = <np.npy_intp> self.max_n_classes
        cdef np.ndarray arr
        arr = np.PyArray_SimpleNewFromData(3, shape, np.NPY_DOUBLE, self.value)
        Py_INCREF(self)
        arr.base = <PyObject*> self
        return arr

    cdef np.ndarray _get_node_ndarray(self):
        """Wraps nodes as a NumPy struct array.

        The array keeps a reference to this Tree, which manages the underlying
        memory. Individual fields are publicly accessible as properties of the
        Tree.
        """
        cdef np.npy_intp shape[1]
        shape[0] = <np.npy_intp> self.node_count
        cdef np.npy_intp strides[1]
        strides[0] = sizeof(Node)
        cdef np.ndarray arr
        Py_INCREF(NODE_DTYPE)
        arr = PyArray_NewFromDescr(np.ndarray, <np.dtype> NODE_DTYPE, 1, shape,
                                   strides, <void*> self.nodes,
                                   np.NPY_DEFAULT, None)
        Py_INCREF(self)
        arr.base = <PyObject*> self
        return arr
