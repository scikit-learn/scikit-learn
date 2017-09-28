# Authors: Gilles Louppe <g.louppe@gmail.com>
#          Peter Prettenhofer <peter.prettenhofer@gmail.com>
#          Brian Holt <bdholt1@gmail.com>
#          Joel Nothman <joel.nothman@gmail.com>
#          Arnaud Joly <arnaud.v.joly@gmail.com>
#          Jacob Schreiber <jmschreiber91@gmail.com>
#          Nelson Liu <nelson@nelsonliu.me>
#
# License: BSD 3 clause

# See _tree.pyx for details.

import numpy as np
cimport numpy as np

from ._splitter cimport Splitter

ctypedef np.npy_float32 DTYPE_t          # Type of X
ctypedef np.npy_float64 DOUBLE_t         # Type of y, sample_weight
ctypedef np.npy_intp SIZE_t              # Type for indices and counters
ctypedef np.npy_int32 INT32_t            # Signed 32 bit integer
ctypedef np.npy_uint32 UINT32_t          # Unsigned 32 bit integer

cdef:
    int TREE_UNDEFINED = -2
    int FEAT_UNKNOWN = -3
    int TREE_LEAF = -1
    bint TREE_NOT_LEAF = 0
    SIZE_t _TREE_LEAF = TREE_LEAF
    SIZE_t _TREE_UNDEFINED = TREE_UNDEFINED
    SIZE_t INITIAL_STACK_SIZE = 10


cdef struct Node:
    # Base storage structure for the nodes in a Tree object

    SIZE_t left_child                    # id of the left child of the node
    SIZE_t right_child                   # id of the right child of the node
    SIZE_t feature                       # Feature used for splitting the node
    DOUBLE_t threshold                   # Threshold value at the node
    DOUBLE_t impurity                    # Impurity of the node (i.e., the value of the criterion)
    SIZE_t n_node_samples                # Number of samples at the node
    DOUBLE_t weighted_n_node_samples     # Weighted number of samples at the node


cdef class Tree:
    # The Tree object is a binary tree structure constructed by the
    # TreeBuilder. The tree structure is used for predictions and
    # feature importances.

    # Input/Output layout
    cdef public SIZE_t n_features        # Number of features in X
    cdef SIZE_t* n_classes               # Number of classes in y[:, k]
    cdef public SIZE_t n_outputs         # Number of outputs in y
    cdef public SIZE_t max_n_classes     # max(n_classes)

    # Inner structures: values are stored separately from node structure,
    # since size is determined at runtime.
    cdef public SIZE_t max_depth         # Max depth of the tree
    cdef public SIZE_t node_count        # Counter for node IDs
    cdef public SIZE_t capacity          # Capacity of tree, in terms of nodes
    cdef Node* nodes                     # Array of nodes
    cdef double* value                   # (capacity, n_outputs, max_n_classes) array of values
    cdef SIZE_t value_stride             # = n_outputs * max_n_classes

    # Methods
    cdef inline SIZE_t _update_node(self, SIZE_t node_id, SIZE_t left_child,
                                    SIZE_t right_child, double threshold,
                                    double impurity, SIZE_t feature,
                                    SIZE_t n_node_samples,
                                    double weighted_n_node_samples):
        cdef Node* node = &self.nodes[node_id]
        node.left_child = left_child
        node.right_child = right_child
        node.threshold = threshold
        node.feature = feature
        node.impurity = impurity
        node.n_node_samples = n_node_samples
        node.weighted_n_node_samples = weighted_n_node_samples

    cdef inline void _set_node_as_leaf(self, Splitter* splitter):
        self._update_node(
            splitter[0].split_record.nid,
            TREE_LEAF, TREE_LEAF,
            TREE_UNDEFINED,
            splitter[0].original_split_record.impurity,
            TREE_UNDEFINED,
            splitter[0].original_split_record.c_stats.n_samples,
            splitter[0].original_split_record.c_stats.sum_weighted_samples)

    cdef inline void _update_parent_splitter(self, Splitter* splitter,
                                             int left_nid, int right_nid):
        self._update_node(
            splitter[0].split_record.nid,
            left_nid, right_nid,
            splitter[0].best_split_record.threshold,
            splitter[0].best_split_record.impurity,
            splitter[0].best_split_record.feature,
            splitter[0].best_split_record.c_stats.n_samples,
            splitter[0].best_split_record.c_stats.sum_weighted_samples)

    cdef inline SIZE_t _add_node(self, SIZE_t parent, bint is_left,
                                 bint is_leaf,
                                 SIZE_t feature, double threshold,
                                 double impurity,
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
            node.left_child = _TREE_UNDEFINED
            node.right_child = _TREE_UNDEFINED

        self.node_count += 1

        return node_id

    cdef inline SIZE_t _add_node_with_value(self, SIZE_t parent, bint is_left,
                                            bint is_leaf,
                                            SIZE_t feature, double threshold,
                                            double impurity,
                                            SIZE_t n_node_samples,
                                            double weighted_n_node_samples,
                                            # FIXME must be an array for multioutput
                                            double node_value):
        node_id = self._add_node(parent, is_left, is_leaf, feature, threshold,
                                 impurity, n_node_samples,
                                 weighted_n_node_samples)
        self.value[node_id * self.value_stride] = node_value
        return node_id

    cdef inline SIZE_t _add_node_set_id(self, int parent_nid,
                                        Splitter* splitter, bint is_left):
        nid = self._add_node_with_value(
            parent_nid,
            is_left, TREE_NOT_LEAF, FEAT_UNKNOWN, TREE_UNDEFINED,
            splitter[0].split_record.impurity,
            splitter[0].split_record.c_stats.n_samples,
            splitter[0].split_record.c_stats.sum_weighted_samples,
            splitter[0].split_record.c_stats.sum_y /
            splitter[0].split_record.c_stats.sum_weighted_samples)

        splitter[0].split_record.nid = nid
        splitter[0].original_split_record.nid = nid
        splitter[0].best_split_record.nid = nid

        return nid

    cdef int _resize(self, SIZE_t capacity) nogil except -1
    cdef int _resize_c(self, SIZE_t capacity=*) nogil except -1

    cdef np.ndarray _get_value_ndarray(self)
    cdef np.ndarray _get_node_ndarray(self)

    cpdef np.ndarray predict(self, object X)

    cpdef np.ndarray apply(self, object X)
    cdef np.ndarray _apply_dense(self, object X)
    cdef np.ndarray _apply_sparse_csr(self, object X)

    cpdef object decision_path(self, object X)
    cdef object _decision_path_dense(self, object X)
    cdef object _decision_path_sparse_csr(self, object X)

    cpdef compute_feature_importances(self, normalize=*)


# =============================================================================
# Tree builder
# =============================================================================

cdef class TreeBuilder:
    # The TreeBuilder recursively builds a Tree object from training samples,
    # using a Splitter object for splitting internal nodes and assigning
    # values to leaves.
    #
    # This class controls the various stopping criteria and the node splitting
    # evaluation order, e.g. depth-first or best-first.

    # cdef Splitter splitter          # Splitting algorithm

    # cdef SIZE_t min_samples_split   # Minimum number of samples in an internal node
    # cdef SIZE_t min_samples_leaf    # Minimum number of samples in a leaf
    # cdef double min_weight_leaf     # Minimum weight in a leaf
    # cdef SIZE_t max_depth           # Maximal tree depth
    # cdef double min_impurity_split  # Impurity threshold for early stopping

    cpdef build(self, Tree tree, float[:, ::1] X, int[::1, :] X_idx_sorted,
                float[::1] y, float[::1] sample_weight,
                float sum_total_weighted_samples)
    # cpdef build(self, Tree tree, object X, np.ndarray y,
    #             np.ndarray sample_weight=*,
    #             np.ndarray X_idx_sorted=*)
    # cdef _check_input(self, object X, np.ndarray y, np.ndarray sample_weight)
