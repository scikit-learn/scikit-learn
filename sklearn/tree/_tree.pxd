# Author: Peter Prettenhofer, Brian Holt, Gilles Louppe
# Licence: BSD 3 clause

# See _tree.pyx for details.

import numpy as np
cimport numpy as np

from cpython cimport bool

ctypedef np.float32_t DTYPE_t               # Type of X
ctypedef np.float64_t DOUBLE_t              # Type of y, sample_weight


# =============================================================================
# Criterion
# =============================================================================

cdef class Criterion:
    # Internal structures
    cdef DOUBLE_t* y                        # Values of y
    cdef Py_ssize_t y_stride                # Stride in y (since n_outputs >= 1)
    cdef DOUBLE_t* sample_weight            # Sample weights

    cdef Py_ssize_t* samples                # Sample indices in X, y
    cdef Py_ssize_t start                   # samples[start:i] are the samples in the left node
    cdef Py_ssize_t i                       # samples[i+1:end] are the samples in the right node
    cdef Py_ssize_t end

    cdef Py_ssize_t n_outputs               # Number of outputs
    cdef Py_ssize_t n_node_samples          # Number of samples in the node (end-start+1)
    cdef double weighted_n_node_samples     # Weighted number of samples
    cdef Py_ssize_t n_left                  # Number of samples in the left node (i-start+1)
    cdef Py_ssize_t n_right                 # Number of samples in the right node (end-(i+1)+1)
    cdef double weighted_n_left             # Weighted number of samples in the left node
    cdef double weighted_n_right            # Weighted number of samples in the right node

    # Methods
    cdef void init(self, DOUBLE_t* y,
                         Py_ssize_t y_stride,
                         DOUBLE_t* sample_weight,
                         Py_ssize_t* samples,
                         Py_ssize_t start,
                         Py_ssize_t i,
                         Py_ssize_t end)
    cdef doubled init_value(self)
    cdef double eval(self)
    cdef bool update(self, Py_ssize_t j)


# =============================================================================
# Splitter
# =============================================================================

cdef class Splitter:
    # Internal structures
    cdef Criterion criterion                # Impurity criterion
    cdef Py_ssize_t* samples                # Sample indices in X, y
    cdef Py_ssize_t n_samples               # Length of samples

    # The samples vector `samples` is maintained by the Splitter object such
    # that the samples contained in a node are contiguous. With this setting,
    # find_split reorganizes the node samples `samples[start:end]` in two
    # subsets `samples[start:i]` and `start[i+1:end]`.

    # Methods
    cdef void init(self, np.ndarray[DTYPE_t, ndim=2] X,
                         np.ndarray[DOUBLE_t, ndim=2, mode="c"] y,
                         np.ndarray[DOUBLE_t, ndim=1, mode="c"] sample_weight
                         Criterion criterion)

    cdef Py_ssize_t find_split(self, Py_ssize_t start, Py_ssize_t end)


# =============================================================================
# Tree
# =============================================================================

cdef class Tree:
    # Input/Output layout
    cdef public Py_ssize_t n_features       # Number of features in X
    cdef Py_ssize_t* n_classes              # Number of classes in y[:, k]
    cdef public Py_ssize_t n_outputs        # Number of outputs in y

    cdef public Py_ssize_t max_n_classes    # max(n_classes)
    cdef public Py_ssize_t value_stride     # n_outputs * max_n_classes

    # Parameters
    cdef public Splitter splitter             # Splitting algorithm
    cdef public Py_ssize_t max_depth          # Max depth of the tree
    cdef public Py_ssize_t min_samples_split  # Minimum number of samples in an internal node
    cdef public Py_ssize_t min_samples_leaf   # Minimum number of samples in a leaf
    cdef public object random_state           # Random state

    # Inner structures
    cdef public Py_ssize_t node_count       # Counter for node IDs
    cdef public Py_ssize_t capacity         # Capacity
    cdef int* children_left                 # children_left[i] is the left child of node i
    cdef int* children_right                # children_right[i] is the right child of node i
    cdef Py_ssize_t* feature                # features[i] is the feature used for splitting node i
    cdef double* threshold                  # threshold[i] is the threshold value at node i
    cdef double* value                      # value[i] is the values contained at node i
    cdef double* impurity                   # impurity[i] is the impurity of node i
    cdef Py_ssize_t* n_node_samples         # n_node_samples[i] is the number of samples at node i

    # Methods
    cpdef build(self, np.ndarray X,
                      np.ndarray y,
                      np.ndarray sample_weight=*)

    cdef int add_split_node(self, int parent, bool is_left_child, Py_ssize_t feature,
                                  double threshold, double* value,
                                  double best_error, double init_error,
                                  Py_ssize_t n_node_samples)

    cdef int add_leaf(self, int parent, int is_left_child, double* value,
                            double error, Py_ssize_t n_node_samples)

    cdef void resize(self, Py_ssize_t capacity=*)
    cpdef predict(self, np.ndarray[DTYPE_t, ndim=2] X)
    cpdef apply(self, np.ndarray[DTYPE_t, ndim=2] X)
    cpdef compute_feature_importances(self)
