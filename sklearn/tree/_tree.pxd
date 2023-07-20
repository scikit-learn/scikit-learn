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

cimport numpy as cnp
from libcpp.unordered_map cimport unordered_map
from libcpp.vector cimport vector

from ._splitter cimport Splitter
from ._utils cimport DOUBLE_t  # Type of y, sample_weight
from ._utils cimport DTYPE_t  # Type of X
from ._utils cimport INT32_t  # Signed 32 bit integer
from ._utils cimport SIZE_t  # Type for indices and counters
from ._utils cimport UINT32_t  # Unsigned 32 bit integer
from ._utils cimport UINT64_t  # Unsigned 64 bit integer
from ._utils cimport SplitValue, SplitRecord, Node


cdef class CategoryCacheMgr:
    # Class to manage the category cache memory during Tree.apply()

    cdef SIZE_t n_nodes
    cdef vector[vector[UINT64_t]] bits

    cdef void populate(
        self,
        Node* nodes,
        SIZE_t n_nodes,
        INT32_t[:] n_categories
    ) noexcept


cdef class BaseTree:
    # Inner structures: values are stored separately from node structure,
    # since size is determined at runtime.
    cdef public SIZE_t max_depth         # Max depth of the tree
    cdef public SIZE_t node_count        # Counter for node IDs
    cdef public SIZE_t capacity          # Capacity of tree, in terms of nodes
    cdef Node* nodes                     # Array of nodes

    cdef SIZE_t value_stride             # The dimensionality of a vectorized output per sample
    cdef double* value                   # Array of values prediction values for each node

    # Generic Methods: These are generic methods used by any tree.
    cdef int _resize(self, SIZE_t capacity) except -1 nogil
    cdef int _resize_c(self, SIZE_t capacity=*) except -1 nogil
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
        Node* node
    ) except -1 nogil
    cdef int _set_leaf_node(
        self,
        SplitRecord* split_node,
        Node* node
    ) except -1 nogil
    cdef DTYPE_t _compute_feature(
        self,
        const DTYPE_t[:, :] X_ndarray,
        SIZE_t sample_index,
        Node *node
    ) noexcept nogil
    cdef void _compute_feature_importances(
        self,
        cnp.float64_t[:] importances,
        Node* node,
    ) noexcept nogil

cdef class Tree(BaseTree):
    # The Supervised Tree object is a binary tree structure constructed by the
    # TreeBuilder. The tree structure is used for predictions and
    # feature importances.
    #
    # Value of upstream properties:
    # - value_stride = n_outputs * max_n_classes
    # - value = (capacity, n_outputs, max_n_classes) array of values

    # Input/Output layout for supervised tree
    cdef public SIZE_t n_features       # Number of features in X
    cdef SIZE_t* n_classes              # Number of classes in y[:, k]
    cdef public SIZE_t n_outputs        # Number of outputs in y
    cdef public SIZE_t max_n_classes    # max(n_classes)

    cdef INT32_t* n_categories          # (n_features,) array of number of categories per feature
    #                                   # is <0 for non-categorial (i.e. -1)

    # Enables the use of tree to store distributions of the output to allow
    # arbitrary usage of the the leaves. This is used in the quantile
    # estimators for example.
    # for storing samples at each leaf node with leaf's node ID as the key and
    # the sample values as the value
    cdef unordered_map[SIZE_t, vector[vector[DOUBLE_t]]] value_samples

    # Methods
    cdef cnp.ndarray _get_value_ndarray(self)
    cdef cnp.ndarray _get_node_ndarray(self)
    cdef cnp.ndarray _get_value_samples_ndarray(self, SIZE_t node_id)
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

    cdef SIZE_t min_samples_split       # Minimum number of samples in an internal node
    cdef SIZE_t min_samples_leaf        # Minimum number of samples in a leaf
    cdef double min_weight_leaf         # Minimum weight in a leaf
    cdef SIZE_t max_depth               # Maximal tree depth
    cdef double min_impurity_decrease   # Impurity threshold for early stopping

    cdef unsigned char store_leaf_values    # Whether to store leaf values

    cpdef build(
        self,
        Tree tree,
        object X,
        const DOUBLE_t[:, ::1] y,
        const DOUBLE_t[:] sample_weight=*,
        const unsigned char[::1] missing_values_in_feature_mask=*,
        const INT32_t[:] n_categories=*,
    )

    cdef _check_input(
        self,
        object X,
        const DOUBLE_t[:, ::1] y,
        const DOUBLE_t[:] sample_weight,
    )
