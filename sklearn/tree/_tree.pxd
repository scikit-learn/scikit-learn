# Author: Peter Prettenhofer, Brian Holt, Gilles Louppe
# Licence: BSD 3 clause

# See _tree.pyx for details.

import numpy as np
cimport numpy as np

from cpython cimport bool

ctypedef np.npy_float32 DTYPE_t          # Type of X
ctypedef np.npy_float64 DOUBLE_t         # Type of y, sample_weight
ctypedef np.npy_intp SIZE_t              # Type for indices and counters


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
    cdef SIZE_t n_left                   # Number of samples in the left node (pos-start)
    cdef SIZE_t n_right                  # Number of samples in the right node (end-pos)
    cdef double weighted_n_left          # Weighted number of samples in the left node
    cdef double weighted_n_right         # Weighted number of samples in the right node

    # The criterion object is maintained such left and right collected
    # statistics correspond to samples[start:pos] and samples[pos:end].

    # Methods
    cdef void init(self, DOUBLE_t* y,
                         SIZE_t y_stride,
                         DOUBLE_t* sample_weight,
                         SIZE_t* samples,
                         SIZE_t start,
                         SIZE_t pos,
                         SIZE_t end)
    cdef void reset(self)
    cdef void update(self, SIZE_t new_pos)
    cdef double node_impurity(self)
    cdef double children_impurity(self)


"""
# =============================================================================
# Splitter
# =============================================================================

cdef class Splitter:
    # Internal structures
    cdef Criterion criterion             # Impurity criterion
    cdef SIZE_t* samples                 # Sample indices in X, y
    cdef SIZE_t n_samples                # Length of samples

    # The samples vector `samples` is maintained by the Splitter object such
    # that the samples contained in a node are contiguous. With this setting,
    # find_split reorganizes the node samples `samples[start:end]` in two
    # subsets `samples[start:pos]` and `start[pos:end]`.

    # Methods
    cdef void init(self, np.ndarray[DTYPE_t, ndim=2] X,
                         np.ndarray[DOUBLE_t, ndim=2, mode="c"] y,
                         np.ndarray[DOUBLE_t, ndim=1, mode="c"] sample_weight,
                         Criterion criterion)

    cdef SIZE_t find_split(self, SIZE_t start, SIZE_t end)


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
    cdef public SIZE_t capacity          # Capacity
    cdef int* children_left              # children_left[i] is the left child of node i
    cdef int* children_right             # children_right[i] is the right child of node i
    cdef SIZE_t* feature                 # features[i] is the feature used for splitting node i
    cdef double* threshold               # threshold[i] is the threshold value at node i
    cdef double* value                   # value[i] is the values contained at node i
    cdef double* impurity                # impurity[i] is the impurity of node i
    cdef SIZE_t* n_node_samples          # n_node_samples[i] is the number of samples at node i

    # Methods
    cpdef build(self, np.ndarray X,
                      np.ndarray y,
                      np.ndarray sample_weight=*)

    cdef int add_split_node(self, int parent, bool is_left_child, SIZE_t feature,
                                  double threshold, double* value,
                                  double best_error, double init_error,
                                  SIZE_t n_node_samples)

    cdef int add_leaf(self, int parent, int is_left_child, double* value,
                            double error, SIZE_t n_node_samples)

    cdef void resize(self, SIZE_t capacity=*)
    cpdef predict(self, np.ndarray[DTYPE_t, ndim=2] X)
    cpdef apply(self, np.ndarray[DTYPE_t, ndim=2] X)
    cpdef compute_feature_importances(self)
"""
