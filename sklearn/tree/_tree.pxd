# Authors: Gilles Louppe <g.louppe@gmail.com>
#          Peter Prettenhofer <peter.prettenhofer@gmail.com>
#          Brian Holt <bdholt1@gmail.com>
# Licence: BSD 3 clause

# See _tree.pyx for details.

import numpy as np
cimport numpy as np

ctypedef np.npy_float32 DTYPE_t          # Type of X
ctypedef np.npy_float64 DOUBLE_t         # Type of y, sample_weight
ctypedef np.npy_intp SIZE_t              # Type for indices and counters
ctypedef np.npy_int32 INT32_t            # Signed 32 bit integer
ctypedef np.npy_uint32 UINT32_t          # Unsigned 32 bit integer


# =============================================================================
# Criterion
# =============================================================================

cdef class Criterion:
    # Internal structures
    cdef DOUBLE_t* y                     # Values of y
    cdef SIZE_t y_stride                 # Stride in y (since n_outputs >= 1)
    cdef DOUBLE_t* sample_weight         # Sample weights

    cdef SIZE_t* samples                 # Sample indices in X, y
    cdef SIZE_t start                    # samples[start:pos] are the samples in the left node
    cdef SIZE_t pos                      # samples[pos:end] are the samples in the right node
    cdef SIZE_t end

    cdef SIZE_t n_outputs                # Number of outputs
    cdef SIZE_t n_node_samples           # Number of samples in the node (end-start)
    cdef double weighted_n_node_samples  # Weighted number of samples
    cdef double weighted_n_left          # Weighted number of samples in the left node
    cdef double weighted_n_right         # Weighted number of samples in the right node

    # The criterion object is maintained such that left and right collected
    # statistics correspond to samples[start:pos] and samples[pos:end].

    # Methods
    cdef void init(self, DOUBLE_t* y,
                         SIZE_t y_stride,
                         DOUBLE_t* sample_weight,
                         SIZE_t* samples,
                         SIZE_t start,
                         SIZE_t end) nogil
    cdef void reset(self) nogil
    cdef void update(self, SIZE_t new_pos) nogil
    cdef double node_impurity(self) nogil
    cdef double children_impurity(self) nogil
    cdef void node_value(self, double* dest) nogil


# =============================================================================
# Splitter
# =============================================================================
cdef class Splitter:
    # Internal structures
    cdef public Criterion criterion      # Impurity criterion
    cdef public SIZE_t max_features      # Number of features to test
    cdef public SIZE_t min_samples_leaf  # Min samples in a leaf

    cdef object random_state             # Random state
    cdef UINT32_t rand_r_state            # sklearn_rand_r random number state

    cdef SIZE_t* samples                 # Sample indices in X, y
    cdef SIZE_t n_samples                # X.shape[0]
    cdef SIZE_t* features                # Feature indices in X
    cdef SIZE_t n_features               # X.shape[1]
    cdef SIZE_t start                    # Start position for the current node
    cdef SIZE_t end                      # End position for the current ndoe

    cdef DTYPE_t* X
    cdef SIZE_t X_stride
    cdef DOUBLE_t* y
    cdef SIZE_t y_stride
    cdef DOUBLE_t* sample_weight

    # The samples vector `samples` is maintained by the Splitter object such
    # that the samples contained in a node are contiguous. With this setting,
    # split reorganizes the node samples `samples[start:end]` in two
    # subsets `samples[start:pos]` and `start[pos:end]`.

    # Methods
    cdef void init(self, np.ndarray X,
                         np.ndarray y,
                         DOUBLE_t* sample_weight)

    cdef void node_reset(self, SIZE_t start, SIZE_t end, double* impurity) nogil

    cdef void node_split(self, SIZE_t* pos, # Set to >= end if the node is a leaf
                               SIZE_t* feature,
                               double* threshold) nogil

    cdef void node_value(self, double* dest) nogil


# =============================================================================
# Tree
# =============================================================================

cdef class Tree:
    # Input/Output layout
    cdef public SIZE_t n_features        # Number of features in X
    cdef SIZE_t* n_classes               # Number of classes in y[:, k]
    cdef public SIZE_t n_outputs         # Number of outputs in y
    cdef public SIZE_t max_n_classes     # max(n_classes)
    cdef public SIZE_t value_stride      # n_outputs * max_n_classes

    # Parameters
    cdef public Splitter splitter        # Splitting algorithm
    cdef public SIZE_t max_depth         # Max depth of the tree
    cdef public SIZE_t min_samples_split # Minimum number of samples in an internal node
    cdef public SIZE_t min_samples_leaf  # Minimum number of samples in a leaf
    cdef public object random_state      # Random state

    # Inner structures
    cdef public SIZE_t node_count        # Counter for node IDs
    cdef public SIZE_t capacity          # Capacity of tree, in terms of nodes
    cdef SIZE_t* children_left           # children_left[i] is the left child of node i
    cdef SIZE_t* children_right          # children_right[i] is the right child of node i
    cdef SIZE_t* feature                 # features[i] is the feature used for splitting node i
    cdef double* threshold               # threshold[i] is the threshold value at node i
    cdef double* value                   # value[i * value_stride:(i+1) * value_stride] are the values contained at node i
    cdef double* impurity                # impurity[i] is the impurity of node i (i.e., the value of the criterion)
    cdef SIZE_t* n_node_samples          # n_node_samples[i] is the number of samples at node i

    # Methods
    cdef SIZE_t _add_node(self, SIZE_t parent,
                                bint is_left,
                                bint is_leaf,
                                SIZE_t feature,
                                double threshold,
                                double impurity,
                                SIZE_t n_node_samples) nogil
    cdef void _resize(self, SIZE_t capacity=*) nogil

    cpdef build(self, np.ndarray X,
                      np.ndarray y,
                      np.ndarray sample_weight=*)

    cpdef predict(self, np.ndarray[DTYPE_t, ndim=2] X)
    cpdef apply(self, np.ndarray[DTYPE_t, ndim=2] X)
    cpdef compute_feature_importances(self, normalize=*)
