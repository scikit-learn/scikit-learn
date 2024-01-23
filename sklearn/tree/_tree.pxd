# Authors: Gilles Louppe <g.louppe@gmail.com>
#          Peter Prettenhofer <peter.prettenhofer@gmail.com>
#          Brian Holt <bdholt1@gmail.com>
#          Joel Nothman <joel.nothman@gmail.com>
#          Arnaud Joly <arnaud.v.joly@gmail.com>
#          Jacob Schreiber <jmschreiber91@gmail.com>
#          Nelson Liu <nelson@nelsonliu.me>
#          Haoyin Xu <haoyinxu@gmail.com>
#
# License: BSD 3 clause

# See _tree.pyx for details.

import numpy as np

cimport numpy as cnp
from libcpp.unordered_map cimport unordered_map
from libcpp.vector cimport vector

from ..utils._typedefs cimport float32_t, float64_t, intp_t, int32_t, uint32_t

from ._utils cimport UINT32_t
from ._splitter cimport SplitRecord, Splitter


cdef struct Node:
    # Base storage structure for the nodes in a Tree object

    intp_t left_child                    # id of the left child of the node
    intp_t right_child                   # id of the right child of the node
    intp_t feature                       # Feature used for splitting the node
    float64_t threshold                  # Threshold value at the node
    float64_t impurity                   # Impurity of the node (i.e., the value of the criterion)
    intp_t n_node_samples                # Number of samples at the node
    float64_t weighted_n_node_samples    # Weighted number of samples at the node
    unsigned char missing_go_to_left     # Whether features have missing values


cdef class BaseTree:

    # Inner structures: values are stored separately from node structure,
    # since size is determined at runtime.
    cdef public intp_t max_depth         # Max depth of the tree
    cdef public intp_t node_count        # Counter for node IDs
    cdef public intp_t capacity          # Capacity of tree, in terms of nodes
    cdef Node* nodes                     # Array of nodes
    cdef float64_t* value                   # (capacity, n_outputs, max_n_classes) array of values
    cdef intp_t value_stride             # = n_outputs * max_n_classes

    # Methods
    cdef int _add_node(
        self,
        intp_t parent,
        bint is_left,
        bint is_leaf,
        SplitRecord* split_node,
        float64_t impurity,
        intp_t n_node_samples,
        float64_t weighted_n_node_samples,
        unsigned char missing_go_to_left
    ) except -1 nogil
    cdef int _resize(self, intp_t capacity) except -1 nogil
    cdef int _resize_c(self, intp_t capacity=*) except -1 nogil

    cdef int _update_node(
        self,
        intp_t parent,
        bint is_left,
        bint is_leaf,
        SplitRecord* split_node,
        float64_t impurity,
        intp_t n_node_samples,
        float64_t weighted_n_node_samples,
        unsigned char missing_go_to_left
    ) except -1 nogil

    # Python API methods: These are methods exposed to Python
    cpdef cnp.ndarray apply(self, object X)
    cdef cnp.ndarray _apply_dense(self, object X)
    cdef cnp.ndarray _apply_sparse_csr(self, object X)

    cpdef object decision_path(self, object X)
    cdef object _decision_path_dense(self, object X)
    cdef object _decision_path_sparse_csr(self, object X)

    cpdef compute_node_depths(self)
    cpdef compute_feature_importances(self, normalize=*)

    # Abstract methods: these functions must be implemented by any decision tree
    cdef int _set_split_node(
        self,
        SplitRecord* split_node,
        Node* node,
        intp_t node_id,
    ) except -1 nogil
    cdef int _set_leaf_node(
        self,
        SplitRecord* split_node,
        Node* node,
        intp_t node_id,
    ) except -1 nogil
    cdef float32_t _compute_feature(
        self,
        const float32_t[:, :] X_ndarray,
        intp_t sample_index,
        Node *node
    ) noexcept nogil
    cdef void _compute_feature_importances(
        self,
        float64_t[:] importances,
        Node* node,
    ) noexcept nogil

cdef class Tree(BaseTree):
    # The Tree object is a binary tree structure constructed by the
    # TreeBuilder. The tree structure is used for predictions and
    # feature importances.

    # The Supervised Tree object is a binary tree structure constructed by the
    # TreeBuilder. The tree structure is used for predictions and
    # feature importances.
    #
    # Value of upstream properties:
    # - value_stride = n_outputs * max_n_classes
    # - value = (capacity, n_outputs, max_n_classes) array of values

    # Input/Output layout
    cdef public intp_t n_features        # Number of features in X
    cdef intp_t* n_classes               # Number of classes in y[:, k]
    cdef public intp_t n_outputs         # Number of outputs in y
    cdef public intp_t max_n_classes     # max(n_classes)

    # Enables the use of tree to store distributions of the output to allow
    # arbitrary usage of the the leaves. This is used in the quantile
    # estimators for example.
    # for storing samples at each leaf node with leaf's node ID as the key and
    # the sample values as the value
    cdef unordered_map[intp_t, vector[vector[float64_t]]] value_samples

    # Methods
    cdef cnp.ndarray _get_value_ndarray(self)
    cdef cnp.ndarray _get_node_ndarray(self)
    cdef cnp.ndarray _get_value_samples_ndarray(self, intp_t node_id)
    cdef cnp.ndarray _get_value_samples_keys(self)

    cpdef cnp.ndarray predict(self, object X)

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

    cdef intp_t min_samples_split           # Minimum number of samples in an internal node
    cdef intp_t min_samples_leaf            # Minimum number of samples in a leaf
    cdef float64_t min_weight_leaf          # Minimum weight in a leaf
    cdef intp_t max_depth                   # Maximal tree depth
    cdef float64_t min_impurity_decrease    # Impurity threshold for early stopping
    cdef cnp.ndarray initial_roots          # Leaf nodes for streaming updates

    cdef unsigned char store_leaf_values    # Whether to store leaf values

    cpdef initialize_node_queue(
      self,
      Tree tree,
      object X,
      const float64_t[:, ::1] y,
      const float64_t[:] sample_weight=*,
      const unsigned char[::1] missing_values_in_feature_mask=*,
    )

    cpdef build(
        self,
        Tree tree,
        object X,
        const float64_t[:, ::1] y,
        const float64_t[:] sample_weight=*,
        const unsigned char[::1] missing_values_in_feature_mask=*,
    )

    cdef _check_input(
        self,
        object X,
        const float64_t[:, ::1] y,
        const float64_t[:] sample_weight,
    )
