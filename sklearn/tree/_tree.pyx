# cython: language_level=3
# cython: boundscheck=False, wraparound=False, initializedcheck=False, cdivision=True

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
#          Haoyin Xu <haoyinxu@gmail.com>
#
# License: BSD 3 clause

from cpython cimport Py_INCREF, PyObject, PyTypeObject
from cython.operator cimport dereference as deref
from libc.math cimport isnan
from libc.stdint cimport INTPTR_MAX
from libc.stdlib cimport free, malloc
from libc.string cimport memcpy, memset
from libcpp cimport bool
from libcpp.algorithm cimport pop_heap, push_heap
from libcpp.vector cimport vector

import struct

import numpy as np

cimport numpy as cnp

cnp.import_array()

from scipy.sparse import issparse
from scipy.sparse import csr_matrix

from ._utils cimport safe_realloc
from ._utils cimport sizet_ptr_to_ndarray


cdef extern from "numpy/arrayobject.h":
    object PyArray_NewFromDescr(PyTypeObject* subtype, cnp.dtype descr,
                                int nd, cnp.npy_intp* dims,
                                cnp.npy_intp* strides,
                                void* data, int flags, object obj)
    int PyArray_SetBaseObject(cnp.ndarray arr, PyObject* obj)

cdef extern from "<stack>" namespace "std" nogil:
    cdef cppclass stack[T]:
        ctypedef T value_type
        stack() except +
        bint empty()
        void pop()
        void push(T&) except +  # Raise c++ exception for bad_alloc -> MemoryError
        T& top()

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

# Build the corresponding numpy dtype for Node.
# This works by casting `dummy` to an array of Node of length 1, which numpy
# can construct a `dtype`-object for. See https://stackoverflow.com/q/62448946
# for a more detailed explanation.
cdef Node dummy
NODE_DTYPE = np.asarray(<Node[:1]>(&dummy)).dtype

# =============================================================================
# TreeBuilder
# =============================================================================

cdef class TreeBuilder:
    """Interface for different tree building strategies."""

    cpdef initialize_node_queue(
        self,
        Tree tree,
        object X,
        const DOUBLE_t[:, ::1] y,
        const DOUBLE_t[:] sample_weight=None,
        const unsigned char[::1] missing_values_in_feature_mask=None,
    ):
        """Build a decision tree from the training set (X, y)."""
        pass

    cpdef build(
        self,
        Tree tree,
        object X,
        const DOUBLE_t[:, ::1] y,
        const DOUBLE_t[:] sample_weight=None,
        const unsigned char[::1] missing_values_in_feature_mask=None,
    ):
        """Build a decision tree from the training set (X, y)."""
        pass

    cdef inline _check_input(
        self,
        object X,
        const DOUBLE_t[:, ::1] y,
        const DOUBLE_t[:] sample_weight,
    ):
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

        # TODO: This check for y seems to be redundant, as it is also
        #  present in the BaseDecisionTree's fit method, and therefore
        #  can be removed.
        if y.base.dtype != DOUBLE or not y.base.flags.contiguous:
            y = np.ascontiguousarray(y, dtype=DOUBLE)

        if (
            sample_weight is not None and
            (
                sample_weight.base.dtype != DOUBLE or
                not sample_weight.base.flags.contiguous
            )
        ):
            sample_weight = np.asarray(sample_weight, dtype=DOUBLE, order="C")

        return X, y, sample_weight

# Depth first builder ---------------------------------------------------------
# A record on the stack for depth-first tree growing
cdef struct StackRecord:
    SIZE_t start
    SIZE_t end
    SIZE_t depth
    SIZE_t parent
    bint is_left
    double impurity
    SIZE_t n_constant_features
    double lower_bound
    double upper_bound

cdef class DepthFirstTreeBuilder(TreeBuilder):
    """Build a decision tree in depth-first fashion."""

    def __cinit__(
        self,
        Splitter splitter,
        SIZE_t min_samples_split,
        SIZE_t min_samples_leaf,
        double min_weight_leaf,
        SIZE_t max_depth,
        double min_impurity_decrease,
        unsigned char store_leaf_values=False,
        cnp.ndarray initial_roots=None,
    ):
        self.splitter = splitter
        self.min_samples_split = min_samples_split
        self.min_samples_leaf = min_samples_leaf
        self.min_weight_leaf = min_weight_leaf
        self.max_depth = max_depth
        self.min_impurity_decrease = min_impurity_decrease
        self.store_leaf_values = store_leaf_values
        self.initial_roots = initial_roots

    def __reduce__(self):
        """Reduce re-implementation, for pickling."""
        return(DepthFirstTreeBuilder, (self.splitter,
                                       self.min_samples_split,
                                       self.min_samples_leaf,
                                       self.min_weight_leaf,
                                       self.max_depth,
                                       self.min_impurity_decrease,
                                       self.store_leaf_values,
                                       self.initial_roots))

    cpdef initialize_node_queue(
        self,
        Tree tree,
        object X,
        const DOUBLE_t[:, ::1] y,
        const DOUBLE_t[:] sample_weight=None,
        const unsigned char[::1] missing_values_in_feature_mask=None,
    ):
        """Initialize a list of roots"""
        X, y, sample_weight = self._check_input(X, y, sample_weight)

        # organize samples by decision paths
        paths = tree.decision_path(X)
        cdef int PARENT
        cdef int CHILD
        cdef int i
        false_roots = {}
        X_copy = {}
        y_copy = {}
        for i in range(X.shape[0]):
            # collect depths from the node paths
            depth_i = paths[i].indices.shape[0] - 1
            PARENT = depth_i - 1
            CHILD = depth_i

            # find leaf node's & their parent node's IDs
            if PARENT < 0:
                parent_i = 0
            else:
                parent_i = paths[i].indices[PARENT]
            child_i = paths[i].indices[CHILD]
            left = 0
            if tree.children_left[parent_i] == child_i:
                left = 1  # leaf node is left child

            # organize samples by the leaf they fall into (false root)
            # leaf nodes are marked by parent node and
            # their relative position (left or right child)
            if (parent_i, left) in false_roots:
                false_roots[(parent_i, left)][0] += 1
                X_copy[(parent_i, left)].append(X[i])
                y_copy[(parent_i, left)].append(y[i])
            else:
                false_roots[(parent_i, left)] = [1, depth_i]
                X_copy[(parent_i, left)] = [X[i]]
                y_copy[(parent_i, left)] = [y[i]]

        X_list = []
        y_list = []

        # reorder the samples according to parent node IDs
        for key, value in reversed(sorted(X_copy.items())):
            X_list = X_list + value
            y_list = y_list + y_copy[key]
        cdef object X_new = np.array(X_list)
        cdef cnp.ndarray y_new = np.array(y_list)

        # initialize the splitter using sorted samples
        cdef Splitter splitter = self.splitter
        splitter.init(X_new, y_new, sample_weight, missing_values_in_feature_mask)

        # convert dict to numpy array and store value
        self.initial_roots = np.array(list(false_roots.items()))

    cpdef build(
        self,
        Tree tree,
        object X,
        const DOUBLE_t[:, ::1] y,
        const DOUBLE_t[:] sample_weight=None,
        const unsigned char[::1] missing_values_in_feature_mask=None,
    ):
        """Build a decision tree from the training set (X, y)."""

        # check input
        X, y, sample_weight = self._check_input(X, y, sample_weight)

        # Parameters
        cdef Splitter splitter = self.splitter
        cdef SIZE_t max_depth = self.max_depth
        cdef SIZE_t min_samples_leaf = self.min_samples_leaf
        cdef double min_weight_leaf = self.min_weight_leaf
        cdef SIZE_t min_samples_split = self.min_samples_split
        cdef double min_impurity_decrease = self.min_impurity_decrease
        cdef unsigned char store_leaf_values = self.store_leaf_values
        cdef cnp.ndarray initial_roots = self.initial_roots

        # Initial capacity
        cdef int init_capacity
        cdef bint first = 0
        if initial_roots is None:
            # Recursive partition (without actual recursion)
            splitter.init(X, y, sample_weight, missing_values_in_feature_mask)

            if tree.max_depth <= 10:
                init_capacity = <int> (2 ** (tree.max_depth + 1)) - 1
            else:
                init_capacity = 2047

            tree._resize(init_capacity)
            first = 1
        else:
            # convert numpy array back to dict
            false_roots = {}
            for key_value_pair in initial_roots:
                false_roots[tuple(key_value_pair[0])] = key_value_pair[1]

            # reset the root array
            self.initial_roots = None

        cdef SIZE_t start = 0
        cdef SIZE_t end = 0
        cdef SIZE_t depth
        cdef SIZE_t parent
        cdef bint is_left
        cdef SIZE_t n_node_samples = splitter.n_samples
        cdef double weighted_n_node_samples
        cdef SIZE_t node_id

        cdef SplitRecord split
        cdef SplitRecord* split_ptr = <SplitRecord *>malloc(splitter.pointer_size())

        cdef double impurity = INFINITY
        cdef double lower_bound
        cdef double upper_bound
        cdef double middle_value
        cdef SIZE_t n_constant_features
        cdef bint is_leaf
        cdef SIZE_t max_depth_seen = -1 if first else tree.max_depth
        cdef int rc = 0

        cdef stack[StackRecord] builder_stack
        cdef stack[StackRecord] update_stack
        cdef StackRecord stack_record

        if not first:
            # push reached leaf nodes onto stack
            for key, value in reversed(sorted(false_roots.items())):
                end += value[0]
                update_stack.push({
                    "start": start,
                    "end": end,
                    "depth": value[1],
                    "parent": key[0],
                    "is_left": key[1],
                    "impurity": tree.impurity[key[0]],
                    "n_constant_features": 0,
                    "lower_bound": -INFINITY,
                    "upper_bound": INFINITY,
                })
                start += value[0]
        else:
            # push root node onto stack
            builder_stack.push({
                "start": 0,
                "end": n_node_samples,
                "depth": 0,
                "parent": _TREE_UNDEFINED,
                "is_left": 0,
                "impurity": INFINITY,
                "n_constant_features": 0,
                "lower_bound": -INFINITY,
                "upper_bound": INFINITY,
            })

        with nogil:
            while not update_stack.empty():
                stack_record = update_stack.top()
                update_stack.pop()

                start = stack_record.start
                end = stack_record.end
                depth = stack_record.depth
                parent = stack_record.parent
                is_left = stack_record.is_left
                impurity = stack_record.impurity
                n_constant_features = stack_record.n_constant_features
                lower_bound = stack_record.lower_bound
                upper_bound = stack_record.upper_bound

                n_node_samples = end - start
                splitter.node_reset(start, end, &weighted_n_node_samples)

                is_leaf = (depth >= max_depth or
                           n_node_samples < min_samples_split or
                           n_node_samples < 2 * min_samples_leaf or
                           weighted_n_node_samples < 2 * min_weight_leaf)

                # impurity == 0 with tolerance due to rounding errors
                is_leaf = is_leaf or impurity <= EPSILON

                if not is_leaf:
                    splitter.node_split(
                        impurity,
                        split_ptr,
                        &n_constant_features,
                        lower_bound,
                        upper_bound
                    )

                    # assign local copy of SplitRecord to assign
                    # pos, improvement, and impurity scores
                    split = deref(split_ptr)

                    # If EPSILON=0 in the below comparison, float precision
                    # issues stop splitting, producing trees that are
                    # dissimilar to v0.18
                    is_leaf = (is_leaf or split.pos >= end or
                               (split.improvement + EPSILON <
                                min_impurity_decrease))

                node_id = tree._update_node(parent, is_left, is_leaf,
                                            split_ptr, impurity, n_node_samples,
                                            weighted_n_node_samples,
                                            split.missing_go_to_left)

                if node_id == INTPTR_MAX:
                    rc = -1
                    break

                # Store value for all nodes, to facilitate tree/model
                # inspection and interpretation
                splitter.node_value(tree.value + node_id * tree.value_stride)
                if splitter.with_monotonic_cst:
                    splitter.clip_node_value(tree.value + node_id * tree.value_stride, lower_bound, upper_bound)

                if not is_leaf:
                    if (
                        not splitter.with_monotonic_cst or
                        splitter.monotonic_cst[split.feature] == 0
                    ):
                        # Split on a feature with no monotonicity constraint

                        # Current bounds must always be propagated to both children.
                        # If a monotonic constraint is active, bounds are used in
                        # node value clipping.
                        left_child_min = right_child_min = lower_bound
                        left_child_max = right_child_max = upper_bound
                    elif splitter.monotonic_cst[split.feature] == 1:
                        # Split on a feature with monotonic increase constraint
                        left_child_min = lower_bound
                        right_child_max = upper_bound

                        # Lower bound for right child and upper bound for left child
                        # are set to the same value.
                        middle_value = splitter.criterion.middle_value()
                        right_child_min = middle_value
                        left_child_max = middle_value
                    else:  # i.e. splitter.monotonic_cst[split.feature] == -1
                        # Split on a feature with monotonic decrease constraint
                        right_child_min = lower_bound
                        left_child_max = upper_bound

                        # Lower bound for left child and upper bound for right child
                        # are set to the same value.
                        middle_value = splitter.criterion.middle_value()
                        left_child_min = middle_value
                        right_child_max = middle_value

                    # Push right child on stack
                    builder_stack.push({
                        "start": split.pos,
                        "end": end,
                        "depth": depth + 1,
                        "parent": node_id,
                        "is_left": 0,
                        "impurity": split.impurity_right,
                        "n_constant_features": n_constant_features,
                        "lower_bound": right_child_min,
                        "upper_bound": right_child_max,
                    })

                    # Push left child on stack
                    builder_stack.push({
                        "start": start,
                        "end": split.pos,
                        "depth": depth + 1,
                        "parent": node_id,
                        "is_left": 1,
                        "impurity": split.impurity_left,
                        "n_constant_features": n_constant_features,
                        "lower_bound": left_child_min,
                        "upper_bound": left_child_max,
                    })
                elif store_leaf_values and is_leaf:
                    # copy leaf values to leaf_values array
                    splitter.node_samples(tree.value_samples[node_id])

                if depth > max_depth_seen:
                    max_depth_seen = depth

            while not builder_stack.empty():
                stack_record = builder_stack.top()
                builder_stack.pop()

                start = stack_record.start
                end = stack_record.end
                depth = stack_record.depth
                parent = stack_record.parent
                is_left = stack_record.is_left
                impurity = stack_record.impurity
                n_constant_features = stack_record.n_constant_features
                lower_bound = stack_record.lower_bound
                upper_bound = stack_record.upper_bound

                n_node_samples = end - start
                splitter.node_reset(start, end, &weighted_n_node_samples)

                is_leaf = (depth >= max_depth or
                           n_node_samples < min_samples_split or
                           n_node_samples < 2 * min_samples_leaf or
                           weighted_n_node_samples < 2 * min_weight_leaf)

                if first:
                    impurity = splitter.node_impurity()
                    first=0

                # impurity == 0 with tolerance due to rounding errors
                is_leaf = is_leaf or impurity <= EPSILON

                if not is_leaf:
                    splitter.node_split(
                        impurity,
                        split_ptr,
                        &n_constant_features,
                        lower_bound,
                        upper_bound
                    )

                    # assign local copy of SplitRecord to assign
                    # pos, improvement, and impurity scores
                    split = deref(split_ptr)

                    # If EPSILON=0 in the below comparison, float precision
                    # issues stop splitting, producing trees that are
                    # dissimilar to v0.18
                    is_leaf = (is_leaf or split.pos >= end or
                               (split.improvement + EPSILON <
                                min_impurity_decrease))

                node_id = tree._add_node(parent, is_left, is_leaf, split_ptr,
                                         impurity, n_node_samples,
                                         weighted_n_node_samples, split.missing_go_to_left)

                if node_id == INTPTR_MAX:
                    rc = -1
                    break

                # Store value for all nodes, to facilitate tree/model
                # inspection and interpretation
                splitter.node_value(tree.value + node_id * tree.value_stride)
                if splitter.with_monotonic_cst:
                    splitter.clip_node_value(tree.value + node_id * tree.value_stride, lower_bound, upper_bound)

                if not is_leaf:
                    if (
                        not splitter.with_monotonic_cst or
                        splitter.monotonic_cst[split.feature] == 0
                    ):
                        # Split on a feature with no monotonicity constraint

                        # Current bounds must always be propagated to both children.
                        # If a monotonic constraint is active, bounds are used in
                        # node value clipping.
                        left_child_min = right_child_min = lower_bound
                        left_child_max = right_child_max = upper_bound
                    elif splitter.monotonic_cst[split.feature] == 1:
                        # Split on a feature with monotonic increase constraint
                        left_child_min = lower_bound
                        right_child_max = upper_bound

                        # Lower bound for right child and upper bound for left child
                        # are set to the same value.
                        middle_value = splitter.criterion.middle_value()
                        right_child_min = middle_value
                        left_child_max = middle_value
                    else:  # i.e. splitter.monotonic_cst[split.feature] == -1
                        # Split on a feature with monotonic decrease constraint
                        right_child_min = lower_bound
                        left_child_max = upper_bound

                        # Lower bound for left child and upper bound for right child
                        # are set to the same value.
                        middle_value = splitter.criterion.middle_value()
                        left_child_min = middle_value
                        right_child_max = middle_value

                    # Push right child on stack
                    builder_stack.push({
                        "start": split.pos,
                        "end": end,
                        "depth": depth + 1,
                        "parent": node_id,
                        "is_left": 0,
                        "impurity": split.impurity_right,
                        "n_constant_features": n_constant_features,
                        "lower_bound": right_child_min,
                        "upper_bound": right_child_max,
                    })

                    # Push left child on stack
                    builder_stack.push({
                        "start": start,
                        "end": split.pos,
                        "depth": depth + 1,
                        "parent": node_id,
                        "is_left": 1,
                        "impurity": split.impurity_left,
                        "n_constant_features": n_constant_features,
                        "lower_bound": left_child_min,
                        "upper_bound": left_child_max,
                    })
                elif store_leaf_values and is_leaf:
                    # copy leaf values to leaf_values array
                    splitter.node_samples(tree.value_samples[node_id])

                if depth > max_depth_seen:
                    max_depth_seen = depth

            if rc >= 0:
                rc = tree._resize_c(tree.node_count)

            if rc >= 0:
                tree.max_depth = max_depth_seen

        # free the memory created for the SplitRecord pointer
        free(split_ptr)

        if rc == -1:
            raise MemoryError()

# Best first builder ----------------------------------------------------------
cdef struct FrontierRecord:
    # Record of information of a Node, the frontier for a split. Those records are
    # maintained in a heap to access the Node with the best improvement in impurity,
    # allowing growing trees greedily on this improvement.
    SIZE_t node_id
    SIZE_t start
    SIZE_t end
    SIZE_t pos
    SIZE_t depth
    bint is_leaf
    double impurity
    double impurity_left
    double impurity_right
    double improvement
    double lower_bound
    double upper_bound
    double middle_value

cdef inline bool _compare_records(
    const FrontierRecord& left,
    const FrontierRecord& right,
):
    return left.improvement < right.improvement

cdef inline void _add_to_frontier(
    FrontierRecord rec,
    vector[FrontierRecord]& frontier,
) noexcept nogil:
    """Adds record `rec` to the priority queue `frontier`."""
    frontier.push_back(rec)
    push_heap(frontier.begin(), frontier.end(), &_compare_records)


cdef class BestFirstTreeBuilder(TreeBuilder):
    """Build a decision tree in best-first fashion.

    The best node to expand is given by the node at the frontier that has the
    highest impurity improvement.
    """
    cdef SIZE_t max_leaf_nodes

    def __cinit__(
        self,
        Splitter splitter,
        SIZE_t min_samples_split,
        SIZE_t min_samples_leaf,
        double min_weight_leaf,
        SIZE_t max_depth,
        SIZE_t max_leaf_nodes,
        double min_impurity_decrease,
        unsigned char store_leaf_values=False,
        cnp.ndarray initial_roots=None,
    ):
        self.splitter = splitter
        self.min_samples_split = min_samples_split
        self.min_samples_leaf = min_samples_leaf
        self.min_weight_leaf = min_weight_leaf
        self.max_depth = max_depth
        self.max_leaf_nodes = max_leaf_nodes
        self.min_impurity_decrease = min_impurity_decrease
        self.store_leaf_values = store_leaf_values
        self.initial_roots = initial_roots

    def __reduce__(self):
        """Reduce re-implementation, for pickling."""
        return(BestFirstTreeBuilder, (self.splitter,
                                      self.min_samples_split,
                                      self.min_samples_leaf,
                                      self.min_weight_leaf,
                                      self.max_depth,
                                      self.max_leaf_nodes,
                                      self.min_impurity_decrease,
                                      self.store_leaf_values,
                                      self.initial_roots))

    cpdef build(
        self,
        Tree tree,
        object X,
        const DOUBLE_t[:, ::1] y,
        const DOUBLE_t[:] sample_weight=None,
        const unsigned char[::1] missing_values_in_feature_mask=None,
    ):
        """Build a decision tree from the training set (X, y)."""

        # check input
        X, y, sample_weight = self._check_input(X, y, sample_weight)

        # Parameters
        cdef Splitter splitter = self.splitter
        cdef SIZE_t max_leaf_nodes = self.max_leaf_nodes
        cdef unsigned char store_leaf_values = self.store_leaf_values

        # Recursive partition (without actual recursion)
        splitter.init(X, y, sample_weight, missing_values_in_feature_mask)

        cdef vector[FrontierRecord] frontier
        cdef FrontierRecord record
        cdef FrontierRecord split_node_left
        cdef FrontierRecord split_node_right
        cdef double left_child_min
        cdef double left_child_max
        cdef double right_child_min
        cdef double right_child_max

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
            rc = self._add_split_node(
                splitter=splitter,
                tree=tree,
                start=0,
                end=n_node_samples,
                impurity=INFINITY,
                is_first=IS_FIRST,
                is_left=IS_LEFT,
                parent=NULL,
                depth=0,
                lower_bound=-INFINITY,
                upper_bound=INFINITY,
                res=&split_node_left,
            )
            if rc >= 0:
                _add_to_frontier(split_node_left, frontier)

            while not frontier.empty():
                pop_heap(frontier.begin(), frontier.end(), &_compare_records)
                record = frontier.back()
                frontier.pop_back()

                node = &tree.nodes[record.node_id]
                is_leaf = (record.is_leaf or max_split_nodes <= 0)

                if is_leaf:
                    # Node is not expandable; set node as leaf
                    node.left_child = _TREE_LEAF
                    node.right_child = _TREE_LEAF
                    node.feature = _TREE_UNDEFINED
                    node.threshold = _TREE_UNDEFINED

                    if store_leaf_values:
                        # copy leaf values to leaf_values array
                        splitter.node_samples(tree.value_samples[record.node_id])
                else:
                    # Node is expandable

                    if (
                        not splitter.with_monotonic_cst or
                        splitter.monotonic_cst[node.feature] == 0
                    ):
                        # Split on a feature with no monotonicity constraint

                        # Current bounds must always be propagated to both children.
                        # If a monotonic constraint is active, bounds are used in
                        # node value clipping.
                        left_child_min = right_child_min = record.lower_bound
                        left_child_max = right_child_max = record.upper_bound
                    elif splitter.monotonic_cst[node.feature] == 1:
                        # Split on a feature with monotonic increase constraint
                        left_child_min = record.lower_bound
                        right_child_max = record.upper_bound

                        # Lower bound for right child and upper bound for left child
                        # are set to the same value.
                        right_child_min = record.middle_value
                        left_child_max = record.middle_value
                    else:  # i.e. splitter.monotonic_cst[split.feature] == -1
                        # Split on a feature with monotonic decrease constraint
                        right_child_min = record.lower_bound
                        left_child_max = record.upper_bound

                        # Lower bound for left child and upper bound for right child
                        # are set to the same value.
                        left_child_min = record.middle_value
                        right_child_max = record.middle_value

                    # Decrement number of split nodes available
                    max_split_nodes -= 1

                    # Compute left split node
                    rc = self._add_split_node(
                        splitter=splitter,
                        tree=tree,
                        start=record.start,
                        end=record.pos,
                        impurity=record.impurity_left,
                        is_first=IS_NOT_FIRST,
                        is_left=IS_LEFT,
                        parent=node,
                        depth=record.depth + 1,
                        lower_bound=left_child_min,
                        upper_bound=left_child_max,
                        res=&split_node_left,
                    )
                    if rc == -1:
                        break

                    # tree.nodes may have changed
                    node = &tree.nodes[record.node_id]

                    # Compute right split node
                    rc = self._add_split_node(
                        splitter=splitter,
                        tree=tree,
                        start=record.pos,
                        end=record.end,
                        impurity=record.impurity_right,
                        is_first=IS_NOT_FIRST,
                        is_left=IS_NOT_LEFT,
                        parent=node,
                        depth=record.depth + 1,
                        lower_bound=right_child_min,
                        upper_bound=right_child_max,
                        res=&split_node_right,
                    )
                    if rc == -1:
                        break

                    # Add nodes to queue
                    _add_to_frontier(split_node_left, frontier)
                    _add_to_frontier(split_node_right, frontier)

                if record.depth > max_depth_seen:
                    max_depth_seen = record.depth

            if rc >= 0:
                rc = tree._resize_c(tree.node_count)

            if rc >= 0:
                tree.max_depth = max_depth_seen

        if rc == -1:
            raise MemoryError()

    cdef inline int _add_split_node(
        self,
        Splitter splitter,
        Tree tree,
        SIZE_t start,
        SIZE_t end,
        double impurity,
        bint is_first,
        bint is_left,
        Node* parent,
        SIZE_t depth,
        double lower_bound,
        double upper_bound,
        FrontierRecord* res
    ) except -1 nogil:
        """Adds node w/ partition ``[start, end)`` to the frontier. """
        cdef SplitRecord split
        cdef SplitRecord* split_ptr = <SplitRecord *>malloc(splitter.pointer_size())

        cdef SIZE_t node_id
        cdef SIZE_t n_node_samples
        cdef SIZE_t n_constant_features = 0
        cdef double min_impurity_decrease = self.min_impurity_decrease
        cdef double weighted_n_node_samples
        cdef bint is_leaf

        splitter.node_reset(start, end, &weighted_n_node_samples)

        if is_first:
            impurity = splitter.node_impurity()

        n_node_samples = end - start
        is_leaf = (depth >= self.max_depth or
                   n_node_samples < self.min_samples_split or
                   n_node_samples < 2 * self.min_samples_leaf or
                   weighted_n_node_samples < 2 * self.min_weight_leaf or
                   impurity <= EPSILON  # impurity == 0 with tolerance
                   )

        if not is_leaf:
            splitter.node_split(
                impurity,
                split_ptr,
                &n_constant_features,
                lower_bound,
                upper_bound
            )
            # assign local copy of SplitRecord to assign
            # pos, improvement, and impurity scores
            split = deref(split_ptr)

            # If EPSILON=0 in the below comparison, float precision issues stop
            # splitting early, producing trees that are dissimilar to v0.18
            is_leaf = (is_leaf or split.pos >= end or
                       split.improvement + EPSILON < min_impurity_decrease)

        node_id = tree._add_node(parent - tree.nodes
                                 if parent != NULL
                                 else _TREE_UNDEFINED,
                                 is_left, is_leaf,
                                 split_ptr, impurity, n_node_samples,
                                 weighted_n_node_samples, split.missing_go_to_left)
        if node_id == INTPTR_MAX:
            return -1

        # compute values also for split nodes (might become leafs later).
        splitter.node_value(tree.value + node_id * tree.value_stride)
        if splitter.with_monotonic_cst:
            splitter.clip_node_value(tree.value + node_id * tree.value_stride, lower_bound, upper_bound)

        res.node_id = node_id
        res.start = start
        res.end = end
        res.depth = depth
        res.impurity = impurity
        res.lower_bound = lower_bound
        res.upper_bound = upper_bound
        res.middle_value = splitter.criterion.middle_value()

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

        free(split_ptr)
        return 0


# =============================================================================
# Tree
# =============================================================================

cdef class BaseTree:
    """Base class for Cython tree models.

    Downstream classes must implement methods to actually traverse the tree.
    """
    cdef int _resize(
        self,
        SIZE_t capacity
    ) except -1 nogil:
        """Resize all inner arrays to `capacity`, if `capacity` == -1, then
           double the size of the inner arrays.

        Returns -1 in case of failure to allocate memory (and raise MemoryError)
        or 0 otherwise.
        """
        if self._resize_c(capacity) != 0:
            # Acquire gil only if we need to raise
            with gil:
                raise MemoryError()

    cdef int _resize_c(
        self,
        SIZE_t capacity=INTPTR_MAX
    ) except -1 nogil:
        """Guts of _resize

        Returns -1 in case of failure to allocate memory (and raise MemoryError)
        or 0 otherwise.
        """
        if capacity == self.capacity and self.nodes != NULL:
            return 0

        if capacity == INTPTR_MAX:
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

    cdef int _set_split_node(
        self,
        SplitRecord* split_node,
        Node* node,
        SIZE_t node_id,
    ) except -1 nogil:
        """Set split node data.

        Parameters
        ----------
        split_node : SplitRecord*
            The pointer to the record of the split node data.
        node : Node*
            The pointer to the node that will hold the split node.
        node_id : SIZE_t
            The index of the node.
        """
        # left_child and right_child will be set later for a split node
        node.feature = split_node.feature
        node.threshold = split_node.threshold
        return 1

    cdef int _set_leaf_node(
        self,
        SplitRecord* split_node,
        Node* node,
        SIZE_t node_id,
    ) except -1 nogil:
        """Set leaf node data.

        Parameters
        ----------
        split_node : SplitRecord*
            The pointer to the record of the leaf node data.
        node : Node*
            The pointer to the node that will hold the leaf node.
        node_id : SIZE_t
            The index of the node.
        """
        node.left_child = _TREE_LEAF
        node.right_child = _TREE_LEAF
        node.feature = _TREE_UNDEFINED
        node.threshold = _TREE_UNDEFINED
        return 1

    cdef DTYPE_t _compute_feature(
        self,
        const DTYPE_t[:, :] X_ndarray,
        SIZE_t sample_index,
        Node *node
    ) noexcept nogil:
        """Compute feature from a given data matrix, X.

        In axis-aligned trees, this is simply the value in the column of X
        for this specific feature.
        """
        # the feature index
        cdef DTYPE_t feature = X_ndarray[sample_index, node.feature]
        return feature

    cdef SIZE_t _add_node(
        self,
        SIZE_t parent,
        bint is_left,
        bint is_leaf,
        SplitRecord* split_node,
        double impurity,
        SIZE_t n_node_samples,
        double weighted_n_node_samples,
        unsigned char missing_go_to_left
    ) except -1 nogil:
        """Add a node to the tree.

        The new node registers itself as the child of its parent.

        Parameters
        ----------
        parent : SIZE_t
            The index of the parent. If '_TREE_UNDEFINED', then the current
            node is a root node.
        is_left : bint
            Whether or not the current node is to the left of the parent node.
        is_leaf : bint
            Whether or not the current node is a leaf node.
        split_node : SplitRecord*
            A pointer to a SplitRecord pointer address.
        impurity : double
            The impurity of the node to be added.
        n_node_samples : SIZE_t
            The number of samples in the node.
        weighted_n_node_samples : double
            The weight of the samples in the node.

        Returns (size_t)(-1) on error.
        """
        cdef SIZE_t node_id = self.node_count

        if node_id >= self.capacity:
            if self._resize_c() != 0:
                return INTPTR_MAX

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
            if self._set_leaf_node(split_node, node, node_id) != 1:
                with gil:
                    raise RuntimeError
        else:
            if self._set_split_node(split_node, node, node_id) != 1:
                with gil:
                    raise RuntimeError
            node.missing_go_to_left = missing_go_to_left

        self.node_count += 1

        return node_id

    cdef inline SIZE_t _update_node(
        self,
        SIZE_t parent,
        bint is_left,
        bint is_leaf,
        SplitRecord* split_node,
        double impurity,
        SIZE_t n_node_samples,
        double weighted_n_node_samples,
        unsigned char missing_go_to_left
    ) except -1 nogil:
        """Update a node on the tree.

        The updated node remains on the same position.
        Returns (size_t)(-1) on error.
        """
        cdef SIZE_t node_id
        if is_left:
            node_id = self.nodes[parent].left_child
        else:
            node_id = self.nodes[parent].right_child

        if node_id >= self.capacity:
            if self._resize_c() != 0:
                return INTPTR_MAX

        cdef Node* node = &self.nodes[node_id]
        node.impurity = impurity
        node.n_node_samples = n_node_samples
        node.weighted_n_node_samples = weighted_n_node_samples

        if is_leaf:
            if self._set_leaf_node(split_node, node, node_id) != 1:
                with gil:
                    raise RuntimeError
        else:
            if self._set_split_node(split_node, node, node_id) != 1:
                with gil:
                    raise RuntimeError
            node.missing_go_to_left = missing_go_to_left

        return node_id

    cpdef cnp.ndarray apply(self, object X):
        """Finds the terminal region (=leaf node) for each sample in X."""
        if issparse(X):
            return self._apply_sparse_csr(X)
        else:
            return self._apply_dense(X)

    cdef inline cnp.ndarray _apply_dense(self, object X):
        """Finds the terminal region (=leaf node) for each sample in X."""

        # Check input
        if not isinstance(X, np.ndarray):
            raise ValueError("X should be in np.ndarray format, got %s"
                             % type(X))

        if X.dtype != DTYPE:
            raise ValueError("X.dtype should be np.float32, got %s" % X.dtype)

        # Extract input
        cdef const DTYPE_t[:, :] X_ndarray = X
        cdef SIZE_t n_samples = X.shape[0]
        cdef DTYPE_t X_i_node_feature

        # Initialize output
        cdef SIZE_t[:] out = np.zeros(n_samples, dtype=np.intp)

        # Initialize auxiliary data-structure
        cdef Node* node = NULL
        cdef SIZE_t i = 0

        with nogil:
            for i in range(n_samples):
                node = self.nodes

                # While node not a leaf
                while node.left_child != _TREE_LEAF:
                    X_i_node_feature = self._compute_feature(X_ndarray, i, node)
                    # ... and node.right_child != _TREE_LEAF:
                    if isnan(X_i_node_feature):
                        if node.missing_go_to_left:
                            node = &self.nodes[node.left_child]
                        else:
                            node = &self.nodes[node.right_child]
                    elif X_i_node_feature <= node.threshold:
                        node = &self.nodes[node.left_child]
                    else:
                        node = &self.nodes[node.right_child]

                out[i] = <SIZE_t>(node - self.nodes)  # node offset

        return np.asarray(out)

    cdef inline cnp.ndarray _apply_sparse_csr(self, object X):
        """Finds the terminal region (=leaf node) for each sample in sparse X.
        """
        # Check input
        if not (issparse(X) and X.format == 'csr'):
            raise ValueError("X should be in csr_matrix format, got %s"
                             % type(X))

        if X.dtype != DTYPE:
            raise ValueError("X.dtype should be np.float32, got %s" % X.dtype)

        # Extract input
        cdef const DTYPE_t[:] X_data = X.data
        cdef const INT32_t[:] X_indices = X.indices
        cdef const INT32_t[:] X_indptr = X.indptr

        cdef SIZE_t n_samples = X.shape[0]
        cdef SIZE_t n_features = X.shape[1]

        # Initialize output
        cdef SIZE_t[:] out = np.zeros(n_samples, dtype=np.intp)

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

                out[i] = <SIZE_t>(node - self.nodes)  # node offset

            # Free auxiliary arrays
            free(X_sample)
            free(feature_to_sample)

        return np.asarray(out)

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
        cdef const DTYPE_t[:, :] X_ndarray = X
        cdef SIZE_t n_samples = X.shape[0]

        # Initialize output
        cdef SIZE_t[:] indptr = np.zeros(n_samples + 1, dtype=np.intp)
        cdef SIZE_t[:] indices = np.zeros(
            n_samples * (1 + self.max_depth), dtype=np.intp
        )

        # Initialize auxiliary data-structure
        cdef Node* node = NULL
        cdef SIZE_t i = 0

        # the feature index
        cdef DOUBLE_t feature

        with nogil:
            for i in range(n_samples):
                node = self.nodes
                indptr[i + 1] = indptr[i]

                # Add all external nodes
                while node.left_child != _TREE_LEAF:
                    # ... and node.right_child != _TREE_LEAF:
                    indices[indptr[i + 1]] = <SIZE_t>(node - self.nodes)
                    indptr[i + 1] += 1

                    # compute the feature value to compare against threshold
                    feature = self._compute_feature(X_ndarray, i, node)
                    if feature <= node.threshold:
                        node = &self.nodes[node.left_child]
                    else:
                        node = &self.nodes[node.right_child]

                # Add the leave node
                indices[indptr[i + 1]] = <SIZE_t>(node - self.nodes)
                indptr[i + 1] += 1

        indices = indices[:indptr[n_samples]]
        cdef SIZE_t[:] data = np.ones(shape=len(indices), dtype=np.intp)
        out = csr_matrix((data, indices, indptr),
                         shape=(n_samples, self.node_count))

        return out

    cdef inline object _decision_path_sparse_csr(self, object X):
        """Finds the decision path (=node) for each sample in X."""

        # Check input
        if not (issparse(X) and X.format == "csr"):
            raise ValueError("X should be in csr_matrix format, got %s"
                             % type(X))

        if X.dtype != DTYPE:
            raise ValueError("X.dtype should be np.float32, got %s" % X.dtype)

        # Extract input
        cdef const DTYPE_t[:] X_data = X.data
        cdef const INT32_t[:] X_indices = X.indices
        cdef const INT32_t[:] X_indptr = X.indptr

        cdef SIZE_t n_samples = X.shape[0]
        cdef SIZE_t n_features = X.shape[1]

        # Initialize output
        cdef SIZE_t[:] indptr = np.zeros(n_samples + 1, dtype=np.intp)
        cdef SIZE_t[:] indices = np.zeros(
            n_samples * (1 + self.max_depth), dtype=np.intp
        )

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
                indptr[i + 1] = indptr[i]

                for k in range(X_indptr[i], X_indptr[i + 1]):
                    feature_to_sample[X_indices[k]] = i
                    X_sample[X_indices[k]] = X_data[k]

                # While node not a leaf
                while node.left_child != _TREE_LEAF:
                    # ... and node.right_child != _TREE_LEAF:

                    indices[indptr[i + 1]] = <SIZE_t>(node - self.nodes)
                    indptr[i + 1] += 1

                    if feature_to_sample[node.feature] == i:
                        feature_value = X_sample[node.feature]

                    else:
                        feature_value = 0.

                    if feature_value <= node.threshold:
                        node = &self.nodes[node.left_child]
                    else:
                        node = &self.nodes[node.right_child]

                # Add the leave node
                indices[indptr[i + 1]] = <SIZE_t>(node - self.nodes)
                indptr[i + 1] += 1

            # Free auxiliary arrays
            free(X_sample)
            free(feature_to_sample)

        indices = indices[:indptr[n_samples]]
        cdef SIZE_t[:] data = np.ones(shape=len(indices), dtype=np.intp)
        out = csr_matrix((data, indices, indptr),
                         shape=(n_samples, self.node_count))

        return out

    cpdef compute_node_depths(self):
        """Compute the depth of each node in a tree.

        .. versionadded:: 1.3

        Returns
        -------
        depths : ndarray of shape (self.node_count,), dtype=np.int64
            The depth of each node in the tree.
        """
        cdef:
            cnp.int64_t[::1] depths = np.empty(self.node_count, dtype=np.int64)
            cnp.npy_intp[:] children_left = self.children_left
            cnp.npy_intp[:] children_right = self.children_right
            cnp.npy_intp node_id
            cnp.npy_intp node_count = self.node_count
            cnp.int64_t depth

        depths[0] = 1  # init root node
        for node_id in range(node_count):
            if children_left[node_id] != _TREE_LEAF:
                depth = depths[node_id] + 1
                depths[children_left[node_id]] = depth
                depths[children_right[node_id]] = depth

        return depths.base

    cpdef compute_feature_importances(self, normalize=True):
        """Computes the importance of each feature (aka variable)."""
        cdef Node* nodes = self.nodes
        cdef Node* node = nodes
        cdef Node* end_node = node + self.node_count

        cdef double normalizer = 0.
        cdef int i = 0

        cdef cnp.float64_t[:] importances = np.zeros(self.n_features)

        with nogil:
            while node != end_node:
                if node.left_child != _TREE_LEAF:
                    # ... and node.right_child != _TREE_LEAF:
                    self._compute_feature_importances(
                        importances, node)

                node += 1

        for i in range(self.n_features):
            importances[i] /= nodes[0].weighted_n_node_samples

        if normalize:
            normalizer = np.sum(importances)

            if normalizer > 0.0:
                # Avoid dividing by zero (e.g., when root is pure)
                for i in range(self.n_features):
                    importances[i] /= normalizer

        return np.asarray(importances)

    cdef void _compute_feature_importances(
        self,
        cnp.float64_t[:] importances,
        Node* node
    ) noexcept nogil:
        """Compute feature importances from a Node in the Tree.

        Wrapped in a private function to allow subclassing that
        computes feature importances.
        """
        cdef Node* nodes = self.nodes
        cdef Node* left
        cdef Node* right

        left = &nodes[node.left_child]
        right = &nodes[node.right_child]

        importances[node.feature] += (
                        node.weighted_n_node_samples * node.impurity -
                        left.weighted_n_node_samples * left.impurity -
                        right.weighted_n_node_samples * right.impurity)

    def compute_partial_dependence(self, DTYPE_t[:, ::1] X,
                                   int[::1] target_features,
                                   double[::1] out):
        """Partial dependence of the response on the ``target_feature`` set.

        For each sample in ``X`` a tree traversal is performed.
        Each traversal starts from the root with weight 1.0.

        At each non-leaf node that splits on a target feature, either
        the left child or the right child is visited based on the feature
        value of the current sample, and the weight is not modified.
        At each non-leaf node that splits on a complementary feature,
        both children are visited and the weight is multiplied by the fraction
        of training samples which went to each child.

        At each leaf, the value of the node is multiplied by the current
        weight (weights sum to 1 for all visited terminal nodes).

        Parameters
        ----------
        X : view on 2d ndarray, shape (n_samples, n_target_features)
            The grid points on which the partial dependence should be
            evaluated.
        target_features : view on 1d ndarray, shape (n_target_features)
            The set of target features for which the partial dependence
            should be evaluated.
        out : view on 1d ndarray, shape (n_samples)
            The value of the partial dependence function on each grid
            point.
        """
        cdef:
            double[::1] weight_stack = np.zeros(self.node_count,
                                                dtype=np.float64)
            SIZE_t[::1] node_idx_stack = np.zeros(self.node_count,
                                                  dtype=np.intp)
            SIZE_t sample_idx
            SIZE_t feature_idx
            int stack_size
            double left_sample_frac
            double current_weight
            double total_weight  # used for sanity check only
            Node *current_node  # use a pointer to avoid copying attributes
            SIZE_t current_node_idx
            bint is_target_feature
            SIZE_t _TREE_LEAF = TREE_LEAF  # to avoid python interactions

        for sample_idx in range(X.shape[0]):
            # init stacks for current sample
            stack_size = 1
            node_idx_stack[0] = 0  # root node
            weight_stack[0] = 1  # all the samples are in the root node
            total_weight = 0

            while stack_size > 0:
                # pop the stack
                stack_size -= 1
                current_node_idx = node_idx_stack[stack_size]
                current_node = &self.nodes[current_node_idx]

                if current_node.left_child == _TREE_LEAF:
                    # leaf node
                    out[sample_idx] += (weight_stack[stack_size] *
                                        self.value[current_node_idx])
                    total_weight += weight_stack[stack_size]
                else:
                    # non-leaf node

                    # determine if the split feature is a target feature
                    is_target_feature = False
                    for feature_idx in range(target_features.shape[0]):
                        if target_features[feature_idx] == current_node.feature:
                            is_target_feature = True
                            break

                    if is_target_feature:
                        # In this case, we push left or right child on stack
                        if X[sample_idx, feature_idx] <= current_node.threshold:
                            node_idx_stack[stack_size] = current_node.left_child
                        else:
                            node_idx_stack[stack_size] = current_node.right_child
                        stack_size += 1
                    else:
                        # In this case, we push both children onto the stack,
                        # and give a weight proportional to the number of
                        # samples going through each branch.

                        # push left child
                        node_idx_stack[stack_size] = current_node.left_child
                        left_sample_frac = (
                            self.nodes[current_node.left_child].weighted_n_node_samples /
                            current_node.weighted_n_node_samples)
                        current_weight = weight_stack[stack_size]
                        weight_stack[stack_size] = current_weight * left_sample_frac
                        stack_size += 1

                        # push right child
                        node_idx_stack[stack_size] = current_node.right_child
                        weight_stack[stack_size] = (
                            current_weight * (1 - left_sample_frac))
                        stack_size += 1

            # Sanity check. Should never happen.
            if not (0.999 < total_weight < 1.001):
                raise ValueError("Total weight should be 1.0 but was %.9f" %
                                 total_weight)


cdef class Tree(BaseTree):
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

    weighted_n_node_samples : array of double, shape [node_count]
        weighted_n_node_samples[i] holds the weighted number of training samples
        reaching node i.

    leaf_node_samples : dict of node id to numpy array of shapes (n_samples_node, n_features)
        A dictionary mapping leaf nodes to the samples of data that are used
        to fit the prediction at each leaf.
    """
    # Wrap for outside world.
    # WARNING: these reference the current `nodes` and `value` buffers, which
    # must not be freed by a subsequent memory allocation.
    # (i.e. through `_resize` or `__setstate__`)
    @property
    def n_classes(self):
        return sizet_ptr_to_ndarray(self.n_classes, self.n_outputs)

    @property
    def children_left(self):
        return self._get_node_ndarray()['left_child'][:self.node_count]

    @property
    def children_right(self):
        return self._get_node_ndarray()['right_child'][:self.node_count]

    @property
    def n_leaves(self):
        return np.sum(np.logical_and(
            self.children_left == -1,
            self.children_right == -1))

    @property
    def feature(self):
        return self._get_node_ndarray()['feature'][:self.node_count]

    @property
    def threshold(self):
        return self._get_node_ndarray()['threshold'][:self.node_count]

    @property
    def impurity(self):
        return self._get_node_ndarray()['impurity'][:self.node_count]

    @property
    def n_node_samples(self):
        return self._get_node_ndarray()['n_node_samples'][:self.node_count]

    @property
    def weighted_n_node_samples(self):
        return self._get_node_ndarray()['weighted_n_node_samples'][:self.node_count]

    @property
    def missing_go_to_left(self):
        return self._get_node_ndarray()['missing_go_to_left'][:self.node_count]

    @property
    def value(self):
        return self._get_value_ndarray()[:self.node_count]

    @property
    def leaf_nodes_samples(self):
        leaf_node_samples = dict()
        keys = self._get_value_samples_keys()
        for node_id in keys:
            leaf_node_samples[node_id] = self._get_value_samples_ndarray(node_id)
        return leaf_node_samples

    # TODO: Convert n_classes to cython.integral memory view once
    #  https://github.com/cython/cython/issues/5243 is fixed
    def __cinit__(self, int n_features, cnp.ndarray n_classes, int n_outputs):
        """Constructor."""
        cdef SIZE_t dummy = 0
        size_t_dtype = np.array(dummy).dtype

        n_classes = _check_n_classes(n_classes, size_t_dtype)

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

        # initialize the hash map for the value samples
        self.value_samples = unordered_map[SIZE_t, vector[vector[DOUBLE_t]]]()

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
        d['value_samples'] = self.leaf_nodes_samples
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

        node_ndarray = _check_node_ndarray(node_ndarray, expected_dtype=NODE_DTYPE)
        value_ndarray = _check_value_ndarray(
            value_ndarray,
            expected_dtype=np.dtype(np.float64),
            expected_shape=value_shape
        )

        self.capacity = node_ndarray.shape[0]
        if self._resize_c(self.capacity) != 0:
            raise MemoryError("resizing tree to %d" % self.capacity)

        memcpy(self.nodes, cnp.PyArray_DATA(node_ndarray),
               self.capacity * sizeof(Node))
        memcpy(self.value, cnp.PyArray_DATA(value_ndarray),
               self.capacity * self.value_stride * sizeof(double))

        # store the leaf node samples if they exist
        value_samples_dict = d['value_samples']
        for node_id, leaf_samples in value_samples_dict.items():
            self.value_samples[node_id].resize(leaf_samples.shape[0])
            for idx in range(leaf_samples.shape[0]):
                for jdx in range(leaf_samples.shape[1]):
                    self.value_samples[node_id][idx].push_back(leaf_samples[idx, jdx])

    cdef cnp.ndarray _get_value_samples_ndarray(self, SIZE_t node_id):
        """Wraps value_samples as a 2-d NumPy array per node_id."""
        cdef int i, j
        cdef int n_samples = self.value_samples[node_id].size()
        cdef cnp.ndarray[DOUBLE_t, ndim=2, mode='c'] leaf_node_samples = np.empty(shape=(n_samples, self.n_outputs), dtype=np.float64)

        for i in range(n_samples):
            for j in range(self.n_outputs):
                leaf_node_samples[i, j] = self.value_samples[node_id][i][j]
        return leaf_node_samples

    cdef cnp.ndarray _get_value_samples_keys(self):
        """Wraps value_samples keys as a 1-d NumPy array of keys."""
        cdef cnp.ndarray[SIZE_t, ndim=1, mode='c'] keys = np.empty(len(self.value_samples), dtype=np.intp)
        cdef unsigned int i = 0

        for key in self.value_samples:
            keys[i] = key.first
            i += 1
        return keys

    cdef cnp.ndarray _get_value_ndarray(self):
        """Wraps value as a 3-d NumPy array.

        The array keeps a reference to this Tree, which manages the underlying
        memory.
        """
        cdef cnp.npy_intp shape[3]
        shape[0] = <cnp.npy_intp> self.node_count
        shape[1] = <cnp.npy_intp> self.n_outputs
        shape[2] = <cnp.npy_intp> self.max_n_classes
        cdef cnp.ndarray arr
        arr = cnp.PyArray_SimpleNewFromData(3, shape, cnp.NPY_DOUBLE, self.value)
        Py_INCREF(self)
        if PyArray_SetBaseObject(arr, <PyObject*> self) < 0:
            raise ValueError("Can't initialize array.")
        return arr

    cdef cnp.ndarray _get_node_ndarray(self):
        """Wraps nodes as a NumPy struct array.

        The array keeps a reference to this Tree, which manages the underlying
        memory. Individual fields are publicly accessible as properties of the
        Tree.
        """
        cdef cnp.npy_intp shape[1]
        shape[0] = <cnp.npy_intp> self.node_count
        cdef cnp.npy_intp strides[1]
        strides[0] = sizeof(Node)
        cdef cnp.ndarray arr
        Py_INCREF(NODE_DTYPE)
        arr = PyArray_NewFromDescr(<PyTypeObject *> cnp.ndarray,
                                   <cnp.dtype> NODE_DTYPE, 1, shape,
                                   strides, <void*> self.nodes,
                                   cnp.NPY_ARRAY_DEFAULT, None)
        Py_INCREF(self)
        if PyArray_SetBaseObject(arr, <PyObject*> self) < 0:
            raise ValueError("Can't initialize array.")
        return arr

    cpdef cnp.ndarray predict(self, object X):
        """Predict target for X."""
        out = self._get_value_ndarray().take(self.apply(X), axis=0,
                                             mode='clip')
        if self.n_outputs == 1:
            out = out.reshape(X.shape[0], self.max_n_classes)
        return out


def _check_n_classes(n_classes, expected_dtype):
    if n_classes.ndim != 1:
        raise ValueError(
            f"Wrong dimensions for n_classes from the pickle: "
            f"expected 1, got {n_classes.ndim}"
        )

    if n_classes.dtype == expected_dtype:
        return n_classes

    # Handles both different endianness and different bitness
    if n_classes.dtype.kind == "i" and n_classes.dtype.itemsize in [4, 8]:
        return n_classes.astype(expected_dtype, casting="same_kind")

    raise ValueError(
        "n_classes from the pickle has an incompatible dtype:\n"
        f"- expected: {expected_dtype}\n"
        f"- got:      {n_classes.dtype}"
    )


def _check_value_ndarray(value_ndarray, expected_dtype, expected_shape):
    if value_ndarray.shape != expected_shape:
        raise ValueError(
            "Wrong shape for value array from the pickle: "
            f"expected {expected_shape}, got {value_ndarray.shape}"
        )

    if not value_ndarray.flags.c_contiguous:
        raise ValueError(
            "value array from the pickle should be a C-contiguous array"
        )

    if value_ndarray.dtype == expected_dtype:
        return value_ndarray

    # Handles different endianness
    if value_ndarray.dtype.str.endswith('f8'):
        return value_ndarray.astype(expected_dtype, casting='equiv')

    raise ValueError(
        "value array from the pickle has an incompatible dtype:\n"
        f"- expected: {expected_dtype}\n"
        f"- got:      {value_ndarray.dtype}"
    )


def _dtype_to_dict(dtype):
    return {name: dt.str for name, (dt, *rest) in dtype.fields.items()}


def _dtype_dict_with_modified_bitness(dtype_dict):
    # field names in Node struct with SIZE_t types (see sklearn/tree/_tree.pxd)
    indexing_field_names = ["left_child", "right_child", "feature", "n_node_samples"]

    expected_dtype_size = str(struct.calcsize("P"))
    allowed_dtype_size = "8" if expected_dtype_size == "4" else "4"

    allowed_dtype_dict = dtype_dict.copy()
    for name in indexing_field_names:
        allowed_dtype_dict[name] = allowed_dtype_dict[name].replace(
            expected_dtype_size, allowed_dtype_size
        )

    return allowed_dtype_dict


def _all_compatible_dtype_dicts(dtype):
    # The Cython code for decision trees uses platform-specific SIZE_t
    # typed indexing fields that correspond to either i4 or i8 dtypes for
    # the matching fields in the numpy array depending on the bitness of
    # the platform (32 bit or 64 bit respectively).
    #
    # We need to cast the indexing fields of the NODE_DTYPE-dtyped array at
    # pickle load time to enable cross-bitness deployment scenarios. We
    # typically want to make it possible to run the expensive fit method of
    # a tree estimator on a 64 bit server platform, pickle the estimator
    # for deployment and run the predict method of a low power 32 bit edge
    # platform.
    #
    # A similar thing happens for endianness, the machine where the pickle was
    # saved can have a different endianness than the machine where the pickle
    # is loaded

    dtype_dict = _dtype_to_dict(dtype)
    dtype_dict_with_modified_bitness = _dtype_dict_with_modified_bitness(dtype_dict)
    dtype_dict_with_modified_endianness = _dtype_to_dict(dtype.newbyteorder())
    dtype_dict_with_modified_bitness_and_endianness = _dtype_dict_with_modified_bitness(
        dtype_dict_with_modified_endianness
    )

    return [
        dtype_dict,
        dtype_dict_with_modified_bitness,
        dtype_dict_with_modified_endianness,
        dtype_dict_with_modified_bitness_and_endianness,
    ]


def _check_node_ndarray(node_ndarray, expected_dtype):
    if node_ndarray.ndim != 1:
        raise ValueError(
            "Wrong dimensions for node array from the pickle: "
            f"expected 1, got {node_ndarray.ndim}"
        )

    if not node_ndarray.flags.c_contiguous:
        raise ValueError(
            "node array from the pickle should be a C-contiguous array"
        )

    node_ndarray_dtype = node_ndarray.dtype
    if node_ndarray_dtype == expected_dtype:
        return node_ndarray

    node_ndarray_dtype_dict = _dtype_to_dict(node_ndarray_dtype)
    all_compatible_dtype_dicts = _all_compatible_dtype_dicts(expected_dtype)

    if node_ndarray_dtype_dict not in all_compatible_dtype_dicts:
        raise ValueError(
            "node array from the pickle has an incompatible dtype:\n"
            f"- expected: {expected_dtype}\n"
            f"- got     : {node_ndarray_dtype}"
        )

    return node_ndarray.astype(expected_dtype, casting="same_kind")


# =============================================================================
# Build Pruned Tree
# =============================================================================


cdef class _CCPPruneController:
    """Base class used by build_pruned_tree_ccp and ccp_pruning_path
    to control pruning.
    """
    cdef bint stop_pruning(self, DOUBLE_t effective_alpha) noexcept nogil:
        """Return 1 to stop pruning and 0 to continue pruning"""
        return 0

    cdef void save_metrics(self, DOUBLE_t effective_alpha,
                           DOUBLE_t subtree_impurities) noexcept nogil:
        """Save metrics when pruning"""
        pass

    cdef void after_pruning(self, unsigned char[:] in_subtree) noexcept nogil:
        """Called after pruning"""
        pass


cdef class _AlphaPruner(_CCPPruneController):
    """Use alpha to control when to stop pruning."""
    cdef DOUBLE_t ccp_alpha
    cdef SIZE_t capacity

    def __cinit__(self, DOUBLE_t ccp_alpha):
        self.ccp_alpha = ccp_alpha
        self.capacity = 0

    cdef bint stop_pruning(self, DOUBLE_t effective_alpha) noexcept nogil:
        # The subtree on the previous iteration has the greatest ccp_alpha
        # less than or equal to self.ccp_alpha
        return self.ccp_alpha < effective_alpha

    cdef void after_pruning(self, unsigned char[:] in_subtree) noexcept nogil:
        """Updates the number of leaves in subtree"""
        for i in range(in_subtree.shape[0]):
            if in_subtree[i]:
                self.capacity += 1


cdef class _PathFinder(_CCPPruneController):
    """Record metrics used to return the cost complexity path."""
    cdef DOUBLE_t[:] ccp_alphas
    cdef DOUBLE_t[:] impurities
    cdef UINT32_t count

    def __cinit__(self,  int node_count):
        self.ccp_alphas = np.zeros(shape=(node_count), dtype=np.float64)
        self.impurities = np.zeros(shape=(node_count), dtype=np.float64)
        self.count = 0

    cdef void save_metrics(self,
                           DOUBLE_t effective_alpha,
                           DOUBLE_t subtree_impurities) noexcept nogil:
        self.ccp_alphas[self.count] = effective_alpha
        self.impurities[self.count] = subtree_impurities
        self.count += 1


cdef struct CostComplexityPruningRecord:
    SIZE_t node_idx
    SIZE_t parent

cdef _cost_complexity_prune(unsigned char[:] leaves_in_subtree,  # OUT
                            Tree orig_tree,
                            _CCPPruneController controller):
    """Perform cost complexity pruning.

    This function takes an already grown tree, `orig_tree` and outputs a
    boolean mask `leaves_in_subtree` which are the leaves in the pruned tree.
    During the pruning process, the controller is passed the effective alpha and
    the subtree impurities. Furthermore, the controller signals when to stop
    pruning.

    Parameters
    ----------
    leaves_in_subtree : unsigned char[:]
        Output for leaves of subtree
    orig_tree : Tree
        Original tree
    ccp_controller : _CCPPruneController
        Cost complexity controller
    """

    cdef:
        SIZE_t i
        SIZE_t n_nodes = orig_tree.node_count
        # prior probability using weighted samples
        DOUBLE_t[:] weighted_n_node_samples = orig_tree.weighted_n_node_samples
        DOUBLE_t total_sum_weights = weighted_n_node_samples[0]
        DOUBLE_t[:] impurity = orig_tree.impurity
        # weighted impurity of each node
        DOUBLE_t[:] r_node = np.empty(shape=n_nodes, dtype=np.float64)

        SIZE_t[:] child_l = orig_tree.children_left
        SIZE_t[:] child_r = orig_tree.children_right
        SIZE_t[:] parent = np.zeros(shape=n_nodes, dtype=np.intp)

        stack[CostComplexityPruningRecord] ccp_stack
        CostComplexityPruningRecord stack_record
        SIZE_t node_idx
        stack[SIZE_t] node_indices_stack

        SIZE_t[:] n_leaves = np.zeros(shape=n_nodes, dtype=np.intp)
        DOUBLE_t[:] r_branch = np.zeros(shape=n_nodes, dtype=np.float64)
        DOUBLE_t current_r
        SIZE_t leaf_idx
        SIZE_t parent_idx

        # candidate nodes that can be pruned
        unsigned char[:] candidate_nodes = np.zeros(shape=n_nodes,
                                                    dtype=np.uint8)
        # nodes in subtree
        unsigned char[:] in_subtree = np.ones(shape=n_nodes, dtype=np.uint8)
        SIZE_t pruned_branch_node_idx
        DOUBLE_t subtree_alpha
        DOUBLE_t effective_alpha
        SIZE_t n_pruned_leaves
        DOUBLE_t r_diff
        DOUBLE_t max_float64 = np.finfo(np.float64).max

    # find parent node ids and leaves
    with nogil:

        for i in range(r_node.shape[0]):
            r_node[i] = (
                weighted_n_node_samples[i] * impurity[i] / total_sum_weights)

        # Push the root node
        ccp_stack.push({"node_idx": 0, "parent": _TREE_UNDEFINED})

        while not ccp_stack.empty():
            stack_record = ccp_stack.top()
            ccp_stack.pop()

            node_idx = stack_record.node_idx
            parent[node_idx] = stack_record.parent

            if child_l[node_idx] == _TREE_LEAF:
                # ... and child_r[node_idx] == _TREE_LEAF:
                leaves_in_subtree[node_idx] = 1
            else:
                ccp_stack.push({"node_idx": child_l[node_idx], "parent": node_idx})
                ccp_stack.push({"node_idx": child_r[node_idx], "parent": node_idx})

        # computes number of leaves in all branches and the overall impurity of
        # the branch. The overall impurity is the sum of r_node in its leaves.
        for leaf_idx in range(leaves_in_subtree.shape[0]):
            if not leaves_in_subtree[leaf_idx]:
                continue
            r_branch[leaf_idx] = r_node[leaf_idx]

            # bubble up values to ancestor nodes
            current_r = r_node[leaf_idx]
            while leaf_idx != 0:
                parent_idx = parent[leaf_idx]
                r_branch[parent_idx] += current_r
                n_leaves[parent_idx] += 1
                leaf_idx = parent_idx

        for i in range(leaves_in_subtree.shape[0]):
            candidate_nodes[i] = not leaves_in_subtree[i]

        # save metrics before pruning
        controller.save_metrics(0.0, r_branch[0])

        # while root node is not a leaf
        while candidate_nodes[0]:

            # computes ccp_alpha for subtrees and finds the minimal alpha
            effective_alpha = max_float64
            for i in range(n_nodes):
                if not candidate_nodes[i]:
                    continue
                subtree_alpha = (r_node[i] - r_branch[i]) / (n_leaves[i] - 1)
                if subtree_alpha < effective_alpha:
                    effective_alpha = subtree_alpha
                    pruned_branch_node_idx = i

            if controller.stop_pruning(effective_alpha):
                break

            node_indices_stack.push(pruned_branch_node_idx)

            # descendants of branch are not in subtree
            while not node_indices_stack.empty():
                node_idx = node_indices_stack.top()
                node_indices_stack.pop()

                if not in_subtree[node_idx]:
                    continue  # branch has already been marked for pruning
                candidate_nodes[node_idx] = 0
                leaves_in_subtree[node_idx] = 0
                in_subtree[node_idx] = 0

                if child_l[node_idx] != _TREE_LEAF:
                    # ... and child_r[node_idx] != _TREE_LEAF:
                    node_indices_stack.push(child_l[node_idx])
                    node_indices_stack.push(child_r[node_idx])
            leaves_in_subtree[pruned_branch_node_idx] = 1
            in_subtree[pruned_branch_node_idx] = 1

            # updates number of leaves
            n_pruned_leaves = n_leaves[pruned_branch_node_idx] - 1
            n_leaves[pruned_branch_node_idx] = 0

            # computes the increase in r_branch to bubble up
            r_diff = r_node[pruned_branch_node_idx] - r_branch[pruned_branch_node_idx]
            r_branch[pruned_branch_node_idx] = r_node[pruned_branch_node_idx]

            # bubble up values to ancestors
            node_idx = parent[pruned_branch_node_idx]
            while node_idx != _TREE_UNDEFINED:
                n_leaves[node_idx] -= n_pruned_leaves
                r_branch[node_idx] += r_diff
                node_idx = parent[node_idx]

            controller.save_metrics(effective_alpha, r_branch[0])

        controller.after_pruning(in_subtree)


def _build_pruned_tree_ccp(
    Tree tree,  # OUT
    Tree orig_tree,
    DOUBLE_t ccp_alpha
):
    """Build a pruned tree from the original tree using cost complexity
    pruning.

    The values and nodes from the original tree are copied into the pruned
    tree.

    Parameters
    ----------
    tree : Tree
        Location to place the pruned tree
    orig_tree : Tree
        Original tree
    ccp_alpha : positive double
        Complexity parameter. The subtree with the largest cost complexity
        that is smaller than ``ccp_alpha`` will be chosen. By default,
        no pruning is performed.
    """

    cdef:
        SIZE_t n_nodes = orig_tree.node_count
        unsigned char[:] leaves_in_subtree = np.zeros(
            shape=n_nodes, dtype=np.uint8)

    pruning_controller = _AlphaPruner(ccp_alpha=ccp_alpha)

    _cost_complexity_prune(leaves_in_subtree, orig_tree, pruning_controller)

    _build_pruned_tree(tree, orig_tree, leaves_in_subtree,
                       pruning_controller.capacity)


def ccp_pruning_path(Tree orig_tree):
    """Computes the cost complexity pruning path.

    Parameters
    ----------
    tree : Tree
        Original tree.

    Returns
    -------
    path_info : dict
        Information about pruning path with attributes:

        ccp_alphas : ndarray
            Effective alphas of subtree during pruning.

        impurities : ndarray
            Sum of the impurities of the subtree leaves for the
            corresponding alpha value in ``ccp_alphas``.
    """
    cdef:
        unsigned char[:] leaves_in_subtree = np.zeros(
            shape=orig_tree.node_count, dtype=np.uint8)

    path_finder = _PathFinder(orig_tree.node_count)

    _cost_complexity_prune(leaves_in_subtree, orig_tree, path_finder)

    cdef:
        UINT32_t total_items = path_finder.count
        DOUBLE_t[:] ccp_alphas = np.empty(shape=total_items, dtype=np.float64)
        DOUBLE_t[:] impurities = np.empty(shape=total_items, dtype=np.float64)
        UINT32_t count = 0

    while count < total_items:
        ccp_alphas[count] = path_finder.ccp_alphas[count]
        impurities[count] = path_finder.impurities[count]
        count += 1

    return {
        'ccp_alphas': np.asarray(ccp_alphas),
        'impurities': np.asarray(impurities),
    }


cdef struct BuildPrunedRecord:
    SIZE_t start
    SIZE_t depth
    SIZE_t parent
    bint is_left

cdef _build_pruned_tree(
    Tree tree,  # OUT
    Tree orig_tree,
    const unsigned char[:] leaves_in_subtree,
    SIZE_t capacity
):
    """Build a pruned tree.

    Build a pruned tree from the original tree by transforming the nodes in
    ``leaves_in_subtree`` into leaves.

    Parameters
    ----------
    tree : Tree
        Location to place the pruned tree
    orig_tree : Tree
        Original tree
    leaves_in_subtree : unsigned char memoryview, shape=(node_count, )
        Boolean mask for leaves to include in subtree
    capacity : SIZE_t
        Number of nodes to initially allocate in pruned tree
    """
    tree._resize(capacity)

    cdef:
        SIZE_t orig_node_id
        SIZE_t new_node_id
        SIZE_t depth
        SIZE_t parent
        bint is_left
        bint is_leaf

        # value_stride for original tree and new tree are the same
        SIZE_t value_stride = orig_tree.value_stride
        SIZE_t max_depth_seen = -1
        int rc = 0
        Node* node
        double* orig_value_ptr
        double* new_value_ptr

        stack[BuildPrunedRecord] prune_stack
        BuildPrunedRecord stack_record

        SplitRecord split

    with nogil:
        # push root node onto stack
        prune_stack.push({"start": 0, "depth": 0, "parent": _TREE_UNDEFINED, "is_left": 0})

        while not prune_stack.empty():
            stack_record = prune_stack.top()
            prune_stack.pop()

            orig_node_id = stack_record.start
            depth = stack_record.depth
            parent = stack_record.parent
            is_left = stack_record.is_left

            is_leaf = leaves_in_subtree[orig_node_id]
            node = &orig_tree.nodes[orig_node_id]

            # redefine to a SplitRecord to pass into _add_node
            split.feature = node.feature
            split.threshold = node.threshold

            new_node_id = tree._add_node(
                parent, is_left, is_leaf, &split,
                node.impurity, node.n_node_samples,
                node.weighted_n_node_samples, node.missing_go_to_left)

            if new_node_id == INTPTR_MAX:
                rc = -1
                break

            # copy value from original tree to new tree
            orig_value_ptr = orig_tree.value + value_stride * orig_node_id
            new_value_ptr = tree.value + value_stride * new_node_id
            memcpy(new_value_ptr, orig_value_ptr, sizeof(double) * value_stride)

            if not is_leaf:
                # Push right child on stack
                prune_stack.push({"start": node.right_child, "depth": depth + 1,
                                  "parent": new_node_id, "is_left": 0})
                # push left child on stack
                prune_stack.push({"start": node.left_child, "depth": depth + 1,
                                  "parent": new_node_id, "is_left": 1})

            if depth > max_depth_seen:
                max_depth_seen = depth

        if rc >= 0:
            tree.max_depth = max_depth_seen
    if rc == -1:
        raise MemoryError("pruning tree")
