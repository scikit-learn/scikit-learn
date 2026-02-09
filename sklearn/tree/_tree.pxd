# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# See _tree.pyx for details.

import numpy as np
cimport numpy as cnp

from sklearn.utils._typedefs cimport float32_t, float64_t, intp_t, int32_t, uint8_t, uint32_t

from sklearn.tree._splitter cimport Splitter
from sklearn.tree._splitter cimport SplitRecord

cdef struct Node:
    # Base storage structure for the nodes in a Tree object

    intp_t left_child                    # id of the left child of the node
    intp_t right_child                   # id of the right child of the node
    intp_t feature                       # Feature used for splitting the node
    float64_t threshold                  # Threshold value at the node
    float64_t impurity                   # Impurity of the node (i.e., the value of the criterion)
    intp_t n_node_samples                # Number of samples at the node
    float64_t weighted_n_node_samples    # Weighted number of samples at the node
    uint8_t missing_go_to_left     # Whether features have missing values


cdef struct ParentInfo:
    # Structure to store information about the parent of a node
    # This is passed to the splitter, to provide information about the previous split

    float64_t lower_bound           # the lower bound of the parent's impurity
    float64_t upper_bound           # the upper bound of the parent's impurity
    float64_t impurity              # the impurity of the parent
    intp_t n_constant_features      # the number of constant features found in parent
    intp_t n_active_interaction_groups  # active interaction groups inherited by child
    intp_t n_allowed_features       # number of allowed features inherited by child

cdef class Tree:
    # The Tree object is a binary tree structure constructed by the
    # TreeBuilder. The tree structure is used for predictions and
    # feature importances.

    # Input/Output layout
    cdef public intp_t n_features        # Number of features in X
    cdef intp_t* n_classes               # Number of classes in y[:, k]
    cdef public intp_t n_outputs         # Number of outputs in y
    cdef public intp_t max_n_classes     # max(n_classes)

    # Inner structures: values are stored separately from node structure,
    # since size is determined at runtime.
    cdef public intp_t max_depth         # Max depth of the tree
    cdef public intp_t node_count        # Counter for node IDs
    cdef public intp_t capacity          # Capacity of tree, in terms of nodes
    cdef Node* nodes                     # Array of nodes
    cdef float64_t* value                # (capacity, n_outputs, max_n_classes) array of values
    cdef intp_t value_stride             # = n_outputs * max_n_classes

    # Methods
    cdef intp_t _add_node(self, intp_t parent, bint is_left, bint is_leaf,
                          intp_t feature, float64_t threshold, float64_t impurity,
                          intp_t n_node_samples,
                          float64_t weighted_n_node_samples,
                          uint8_t missing_go_to_left) except -1 nogil
    cdef int _resize(self, intp_t capacity) except -1 nogil
    cdef int _resize_c(self, intp_t capacity=*) except -1 nogil

    cdef cnp.ndarray _get_value_ndarray(self)
    cdef cnp.ndarray _get_node_ndarray(self)

    cpdef cnp.ndarray predict(self, object X)

    cpdef cnp.ndarray apply(self, object X)
    cdef cnp.ndarray _apply_dense(self, object X)
    cdef cnp.ndarray _apply_sparse_csr(self, object X)

    cpdef object decision_path(self, object X)
    cdef object _decision_path_dense(self, object X)
    cdef object _decision_path_sparse_csr(self, object X)

    cpdef compute_node_depths(self)
    cpdef compute_feature_importances(self, normalize=*)


# =============================================================================
# Tree builder
# =============================================================================

cdef class InteractionConstraints:
    cdef const intp_t[:] feature_to_groups_indptr
    cdef const intp_t[:] feature_to_groups_indices
    cdef const intp_t[:] group_to_features_indptr
    cdef const intp_t[:] group_to_features_indices
    cdef intp_t[::1] interaction_groups
    cdef int32_t[::1] group_marks
    cdef int32_t[::1] feature_marks
    cdef int32_t group_mark_token
    cdef int32_t feature_mark_token
    cdef intp_t n_interaction_groups

    cdef void init_fit_state(self, Splitter splitter)
    cdef void reset_node_state(
        self,
        ParentInfo* parent_record,
        intp_t n_features,
    ) noexcept nogil
    cdef void update_after_split(
        self,
        intp_t[::1] features,
        ParentInfo* parent_record,
        intp_t split_feature,
    ) noexcept nogil
    cdef void restore_state_from_path(
        self,
        intp_t[::1] features,
        ParentInfo* parent_record,
        Node* nodes,
        intp_t[::1] parent_node_ids,
        intp_t parent_node_id,
        intp_t n_features,
    ) noexcept nogil


cdef class TreeBuilder:
    # The TreeBuilder recursively builds a Tree object from training samples,
    # using a Splitter object for splitting internal nodes and assigning
    # values to leaves.
    #
    # This class controls the various stopping criteria and the node splitting
    # evaluation order, e.g. depth-first or best-first.

    cdef Splitter splitter              # Splitting algorithm

    cdef intp_t min_samples_split       # Minimum number of samples in an internal node
    cdef intp_t min_samples_leaf        # Minimum number of samples in a leaf
    cdef float64_t min_weight_leaf         # Minimum weight in a leaf
    cdef intp_t max_depth               # Maximal tree depth
    cdef float64_t min_impurity_decrease   # Impurity threshold for early stopping
    cdef InteractionConstraints interaction_constraints

    cpdef build(
        self,
        Tree tree,
        object X,
        const float64_t[:, ::1] y,
        const float64_t[:] sample_weight=*,
        const uint8_t[::1] missing_values_in_feature_mask=*,
    )

    cdef _check_input(
        self,
        object X,
        const float64_t[:, ::1] y,
        const float64_t[:] sample_weight,
    )


# =============================================================================
# Tree pruning
# =============================================================================

# The private function allows any external caller to prune the tree and return
# a new tree with the pruned nodes. The pruned tree is a new tree object.
#
# .. warning:: this function is not backwards compatible and may change without
#              notice.
cdef void _build_pruned_tree(
    Tree tree,  # OUT
    Tree orig_tree,
    const uint8_t[:] leaves_in_subtree,
    intp_t capacity
)
