# Authors: Gilles Louppe <g.louppe@gmail.com>
#          Peter Prettenhofer <peter.prettenhofer@gmail.com>
#          Brian Holt <bdholt1@gmail.com>
#          Joel Nothman <joel.nothman@gmail.com>
#          Arnaud Joly <arnaud.v.joly@gmail.com>
#
# Licence: BSD 3 clause

# See _tree.pyx for details.

import numpy as np
cimport numpy as np

ctypedef np.npy_float32 DTYPE_t          # Type of X
ctypedef np.npy_float64 DOUBLE_t         # Type of y, sample_weight
ctypedef np.npy_intp SIZE_t              # Type for indices and counters
ctypedef np.npy_int32 INT32_t            # Signed 32 bit integer
ctypedef np.npy_uint32 UINT32_t          # Unsigned 32 bit integer

from sklearn.tree._utils cimport SplitRecord

# =============================================================================
# Criterion
# =============================================================================

cdef class Criterion:
    # The criterion computes the impurity of a node and the reduction of
    # impurity of a split on that node. It also computes the output statistics
    # such as the mean in regression and class probabilities in classification.

    cdef SIZE_t n_outputs

    cdef DOUBLE_t* y             
    cdef DOUBLE_t* w        

    cdef DOUBLE_t* node_value_left
    cdef DOUBLE_t* node_value_right

    cdef DOUBLE_t* yw_cl
    cdef DOUBLE_t* yw_cr

    cdef public SIZE_t max_classes
    cdef public SIZE_t y_stride

    cdef SIZE_t min_leaf_samples
    cdef DOUBLE_t min_leaf_weight
    cdef SIZE_t n_samples
    cdef DOUBLE_t weighted_n_samples

    # The criterion object is maintained such that left and right collected
    # statistics correspond to samples[start:pos] and samples[pos:end].

    # Methods
    cdef void init(self, DOUBLE_t* y, DOUBLE_t* w, SIZE_t n_samples, 
        SIZE_t min_leaf_samples, DOUBLE_t min_leaf_weight, DOUBLE_t* w_sum,
        DOUBLE_t* yw_sq_sum, DOUBLE_t** node_value)

    cdef SplitRecord best_split(self, DTYPE_t* X_i, SIZE_t* samples, 
        SIZE_t start, SIZE_t end, SIZE_t feature, DOUBLE_t w_sum, 
        DOUBLE_t yw_sq_sum, DOUBLE_t* node_value) nogil

    cdef SplitRecord random_split(self, DTYPE_t* X_i, SIZE_t* samples, 
        SIZE_t start, SIZE_t end, SIZE_t feature, DOUBLE_t w_sum, 
        DOUBLE_t yw_sq_sum, DOUBLE_t* node_value, UINT32_t* rand_r) nogil

# =============================================================================
# Splitter
# =============================================================================

cdef class Splitter:
    # The splitter searches in the input space for a feature and a threshold
    # to split the samples samples[start:end].
    #
    # The impurity computations are delegated to a criterion object.

    # Internal structures
    cdef public Criterion criterion      # Impurity criterion
    cdef public SIZE_t max_features      # Number of features to test
    cdef public SIZE_t min_samples_leaf  # Min samples in a leaf
    cdef public double min_weight_leaf   # Minimum weight in a leaf

    cdef object random_state             # Random state
    cdef UINT32_t rand_r_state           # sklearn_rand_r random number state

    cdef SIZE_t* samples                 # Sample indices in X, y
    cdef SIZE_t n_samples
    cdef SIZE_t* sample_mask
    cdef SIZE_t* features                # Feature indices in X
    cdef SIZE_t n_features               # X.shape[1]

    cdef DOUBLE_t* y
    cdef SIZE_t y_stride
    cdef DOUBLE_t* sample_weight

    cdef np.ndarray X_idx_sorted
    cdef INT32_t* X_idx_sorted_ptr
    cdef SIZE_t X_idx_sorted_stride

    cdef DTYPE_t* X
    cdef DTYPE_t* X_i
    cdef SIZE_t X_sample_stride
    cdef SIZE_t X_feature_stride

    cdef SIZE_t best
    cdef bint presort

    # The samples vector `samples` is maintained by the Splitter object such
    # that the samples contained in a node are contiguous. With this setting,
    # `node_split` reorganizes the node samples `samples[start:end]` in two
    # subsets `samples[start:pos]` and `samples[pos:end]`.

    # The 1-d  `features` array of size n_features contains the features
    # indices and allows fast sampling without replacement of features.

    # The 1-d `constant_features` array of size n_features holds in
    # `constant_features[:n_constant_features]` the feature ids with
    # constant values for all the samples that reached a specific node.
    # The value `n_constant_features` is given by the the parent node to its
    # child nodes.  The content of the range `[n_constant_features:]` is left
    # undefined, but preallocated for performance reasons
    # This allows optimization with depth-based tree building.

    # Methods
    cdef void init(self, object X, np.ndarray X_idx_sorted, bint presort,
                   np.ndarray y, DOUBLE_t* sample_weight, DOUBLE_t* w_sum, 
                   DOUBLE_t* yw_sq_sum, DOUBLE_t** node_value, 
                   SIZE_t* n_node_samples)

    cdef SplitRecord _split(self, SIZE_t start, SIZE_t end,
        SIZE_t feature, DOUBLE_t w_sum, DOUBLE_t yw_sq, DOUBLE_t* node_value, 
        SIZE_t best) nogil

    cdef SplitRecord split(self, SIZE_t start, SIZE_t end, DOUBLE_t w_sum, 
        DOUBLE_t yw_sq, DOUBLE_t* node_value) nogil


# =============================================================================
# Tree
# =============================================================================

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
    cdef SIZE_t _add_node(self, SIZE_t parent, bint is_left, bint is_leaf,
                          SIZE_t feature, double threshold, double impurity,
                          SIZE_t n_node_samples,
                          double weighted_n_samples, DOUBLE_t* value) nogil
    cdef void _resize(self, SIZE_t capacity) except *
    cdef int _resize_c(self, SIZE_t capacity=*) nogil

    cdef np.ndarray _get_value_ndarray(self)
    cdef np.ndarray _get_node_ndarray(self)

    cpdef np.ndarray predict(self, object X)
    cpdef np.ndarray apply(self, object X)
    cdef np.ndarray _apply_dense(self, object X)
    cdef np.ndarray _apply_sparse_csr(self, object X)

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

    cdef Splitter splitter          # Splitting algorithm

    cdef SIZE_t min_samples_split   # Minimum number of samples in an internal node
    cdef SIZE_t min_samples_leaf    # Minimum number of samples in a leaf
    cdef double min_weight_leaf     # Minimum weight in a leaf
    cdef SIZE_t max_depth           # Maximal tree depth
    cdef SIZE_t max_leaf_nodes      # Maximal number of leaf nodes the tree
                                    # can have, used in best first splitting

    cdef _check_input(self, object X, np.ndarray y, np.ndarray sample_weight)

    cpdef depth_first(self, Tree tree, object X, np.ndarray y,
                np.ndarray sample_weight, bint presort,
                np.ndarray X_idx_sorted)

    cpdef best_first(self, Tree tree, object X, np.ndarray y,
                np.ndarray sample_weight, bint presort,
                np.ndarray X_idx_sorted)
