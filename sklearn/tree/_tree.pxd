# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# See _tree.pyx for details.

import numpy as np

cimport numpy as cnp

from ..ensemble._hist_gradient_boosting.common cimport BITSET_INNER_DTYPE_C
from ..utils._typedefs cimport float32_t, float64_t, int32_t, intp_t, uint32_t, BITSET_t
from ._utils cimport ParentInfo
from ._partitioner cimport SplitValue
from ._splitter cimport SplitRecord, Splitter


cdef struct Node:
    # Base storage structure for the nodes in a Tree object

    intp_t left_child                    # id of the left child of the node
    intp_t right_child                   # id of the right child of the node
    intp_t feature                       # Feature used for splitting the node
    # SplitValue split_value             # Generalized threshold for categorical and
    #                                    # non-categorical features
    float64_t threshold
    BITSET_t cat_split
    float64_t impurity                   # Impurity of the node (i.e., the value of the criterion)
    intp_t n_node_samples                # Number of samples at the node
    float64_t weighted_n_node_samples    # Weighted number of samples at the node
    unsigned char missing_go_to_left     # Whether features have missing values


cdef class Tree:
    # The Tree object is a binary tree structure constructed by the
    # TreeBuilder. The tree structure is used for predictions and
    # feature importances.

    # Input/Output layout
    cdef public intp_t n_features        # Number of features in X
    cdef intp_t* n_classes               # Number of classes in y[:, k]
    cdef public intp_t n_outputs         # Number of outputs in y
    cdef public intp_t max_n_classes     # max(n_classes)
    cdef int32_t *n_categories           # (n_features,) array giving number of
    #                                    # categories (<0 for non-categorical)

    # Inner structures: values are stored separately from node structure,
    # since size is determined at runtime.
    cdef public intp_t max_depth         # Max depth of the tree
    cdef public intp_t node_count        # Counter for node IDs
    cdef public intp_t capacity          # Capacity of tree, in terms of nodes
    cdef Node* nodes                     # Array of nodes
    cdef float64_t* value                # (capacity, n_outputs, max_n_classes) array of values
    cdef intp_t value_stride             # = n_outputs * max_n_classes

    # Methods
    cdef intp_t _add_node(
        self,
        intp_t parent,
        bint is_left,
        bint is_leaf,
        intp_t feature,
        SplitValue split_value,
        # float64_t threshold,
        float64_t impurity,
        intp_t n_node_samples,
        float64_t weighted_n_node_samples,
        unsigned char missing_go_to_left
    ) except -1 nogil
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

    cpdef build(
        self,
        Tree tree,
        object X,
        const float64_t[:, ::1] y,
        const float64_t[:] sample_weight=*,
        const unsigned char[::1] missing_values_in_feature_mask=*,
        const int32_t[::1] n_categories=*,
    )

    cdef _check_input(
        self,
        object X,
        const float64_t[:, ::1] y,
        const float64_t[:] sample_weight,
    )


# cdef class CategoryCacheMgr:
#     # Class to manage the category cache memory during Tree.apply()

#     cdef intp_t n_nodes
#     cdef BITSET_INNER_DTYPE_C **bits

#     cdef void populate(
#         self,
#         Node *nodes,
#         intp_t n_nodes,
#         int32_t *n_categories
#     )
