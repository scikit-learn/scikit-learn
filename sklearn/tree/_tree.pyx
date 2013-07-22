# encoding: utf-8
# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
#
# Author: Peter Prettenhofer, Brian Holt, Gilles Louppe, Noel Dawe
#
# Licence: BSD 3 clause


# =============================================================================
# Imports
# =============================================================================

cimport cython
from cpython cimport bool
from libc.float cimport DBL_MAX
from libc.math cimport log, pow
from libc.stdlib cimport calloc, free, malloc, realloc
from libc.string cimport memcpy

import numpy as np
cimport numpy as np
np.import_array()

from numpy import zeros as np_zeros
from numpy import ones as np_ones
from numpy import bool as np_bool
from numpy import float32 as np_float32
from numpy import float64 as np_float64


# =============================================================================
# Types and constants
# =============================================================================

# Dtype
DTYPE = np_float32
DOUBLE = np_float64
# ctypedef np.float32_t DTYPE_t
# ctypedef np.float64_t DOUBLE_t
# ctypedef np.int8_t BOOL_t

# Constants
cdef double INFINITY = np.inf

TREE_LEAF = -1
TREE_UNDEFINED = -2
cdef int _TREE_LEAF = TREE_LEAF
cdef int _TREE_UNDEFINED = TREE_UNDEFINED

TREE_SPLIT_BEST = 1
TREE_SPLIT_RANDOM = 2
cdef int _TREE_SPLIT_BEST = TREE_SPLIT_BEST
cdef int _TREE_SPLIT_RANDOM = TREE_SPLIT_RANDOM


# =============================================================================
# Tree
# =============================================================================

cdef class Tree:
    """Struct-of-arrays representation of a binary decision tree.

    The binary tree is represented as a number of parallel arrays.
    The i-th element of each array holds information about the
    node `i`. Node 0 is the tree's root.
    You can find a detailed description of all arrays
    below. NOTE: Some of the arrays only apply to either leaves or
    split nodes, resp. In this case the values of nodes of the other
    type are arbitrary!

    Parameters
    ----------
    n_features : int
        The number of features

    n_classes : array-like
        n_classes[k] is the number of classes for output k.

    n_outputs : int
        The number of outputs.

    criterion : Criterion
    max_depth : double
    min_samples_split : int
    min_samples_leaf : int
    min_density : double
    max_features : int
    find_split_algorithm : int

    Attributes
    ----------
    node_count : int
        The number of nodes (internal nodes + leaves) in the tree.

    capacity : int
        The current capacity (i.e., size) of the arrays.

    children_left : int*
        children_left[i] holds the node id of the left child of node i.
        For leaves, children_left[i] == TREE_LEAF. Otherwise,
        children_left[i] > i. This child handles the case where
        X[:, feature[i]] <= threshold[i].

    children_right : int*
        children_right[i] holds the node id of the right child of node i.
        For leaves, children_right[i] == TREE_LEAF. Otherwise,
        children_right[i] > i. This child handles the case where
        X[:, feature[i]] > threshold[i].

    feature : int*
        feature[i] holds the feature to split on, for the internal node i.

    threshold : double*
        threshold[i] holds the threshold for the internal node i.

    value : double*
        Contains the constant prediction value of each node.

    best_error : double*
        best_error[i] holds the error of the (best) split at node i.
        For leaves init_error[i] == best_error[i].

    init_error : double*
        init_error[i] holds the initial error at node i (before splitting).
        For leaves init_error[i] == best_error[i].

    n_samples : int*
        n_samples[i] holds the number of training samples reaching node i.
    """

    # # Input/Output layout
    # cdef public int n_features
    # cdef int* n_classes
    # cdef public int n_outputs

    # cdef public int max_n_classes
    # cdef public Py_ssize_t value_stride

    # # Parameters
    # cdef public Criterion criterion
    # cdef public double max_depth
    # cdef public int min_samples_split
    # cdef public int min_samples_leaf
    # cdef public double min_density
    # cdef public int max_features
    # cdef public int find_split_algorithm
    # cdef public object random_state

    # # Inner structures
    # cdef public int node_count
    # cdef public int capacity
    # cdef int* children_left
    # cdef int* children_right
    # cdef int* feature
    # cdef double* threshold
    # cdef double* value
    # cdef double* best_error
    # cdef double* init_error
    # cdef int* n_samples

    # Wrap for outside world
    property n_classes:
        def __get__(self):
            return intp_to_ndarray(self.n_classes, self.n_outputs)

    property children_left:
        def __get__(self):
            return intp_to_ndarray(self.children_left, self.node_count)

    property children_right:
        def __get__(self):
            return intp_to_ndarray(self.children_right, self.node_count)

    property feature:
        def __get__(self):
            return intp_to_ndarray(self.feature, self.node_count)

    property threshold:
        def __get__(self):
            return doublep_to_ndarray(self.threshold, self.node_count)

    property value:
        def __get__(self):
            cdef np.npy_intp shape[3]

            shape[0] = <np.npy_intp> self.node_count
            shape[1] = <np.npy_intp> self.n_outputs
            shape[2] = <np.npy_intp> self.max_n_classes

            return np.PyArray_SimpleNewFromData(
                3, shape, np.NPY_DOUBLE, self.value)

    property best_error:
        def __get__(self):
            return doublep_to_ndarray(self.best_error, self.node_count)

    property init_error:
        def __get__(self):
            return doublep_to_ndarray(self.init_error, self.node_count)

    property n_samples:
        def __get__(self):
            return intp_to_ndarray(self.n_samples, self.node_count)

    def __cinit__(self, int n_features, object n_classes, int n_outputs,
                 Criterion criterion, double max_depth, int min_samples_split,
                 int min_samples_leaf, double min_density, int max_features,
                 int find_split_algorithm, object random_state):
        """Constructor."""
        # Input/Output layout
        cdef int k

        self.n_features = n_features
        self.n_outputs = n_outputs
        self.n_classes = <int*> malloc(n_outputs * sizeof(int))

        if self.n_classes == NULL:
            raise MemoryError()

        self.max_n_classes = np.max(n_classes)
        self.value_stride = self.n_outputs * self.max_n_classes

        for k from 0 <= k < n_outputs:
            self.n_classes[k] = n_classes[k]

        # Parameters
        self.criterion = criterion
        self.max_depth = max_depth
        self.min_samples_split = min_samples_split
        self.min_samples_leaf = min_samples_leaf
        self.min_density = min_density
        self.max_features = max_features
        self.find_split_algorithm = find_split_algorithm
        self.random_state = random_state

        # Inner structures
        self.node_count = 0
        self.capacity = 0

        self.children_left = NULL
        self.children_right = NULL
        self.feature = NULL
        self.threshold = NULL
        self.value = NULL
        self.best_error = NULL
        self.init_error = NULL
        self.n_samples = NULL

        self.features = np.arange(n_features, dtype=np.int32)

    def __dealloc__(self):
        """Destructor."""
        # Free all inner structures
        free(self.n_classes)
        free(self.children_left)
        free(self.children_right)
        free(self.feature)
        free(self.threshold)
        free(self.value)
        free(self.best_error)
        free(self.init_error)
        free(self.n_samples)

    def __reduce__(self):
        """Reduce re-implementation, for pickling."""
        return (Tree, (self.n_features,
                       intp_to_ndarray(self.n_classes, self.n_outputs),
                       self.n_outputs,
                       self.criterion,
                       self.max_depth,
                       self.min_samples_split,
                       self.min_samples_leaf,
                       self.min_density,
                       self.max_features,
                       self.find_split_algorithm,
                       self.random_state), self.__getstate__())

    def __getstate__(self):
        """Getstate re-implementation, for pickling."""
        d = {}

        d["node_count"] = self.node_count
        d["capacity"] = self.capacity
        d["children_left"] = intp_to_ndarray(self.children_left, self.capacity)
        d["children_right"] = intp_to_ndarray(self.children_right, self.capacity)
        d["feature"] = intp_to_ndarray(self.feature, self.capacity)
        d["threshold"] = doublep_to_ndarray(self.threshold, self.capacity)
        d["value"] = doublep_to_ndarray(self.value, self.capacity * self.value_stride)
        d["best_error"] = doublep_to_ndarray(self.best_error, self.capacity)
        d["init_error"] = doublep_to_ndarray(self.init_error, self.capacity)
        d["n_samples"] = intp_to_ndarray(self.n_samples, self.capacity)

        return d

    def __setstate__(self, d):
        """Setstate re-implementation, for unpickling."""
        self.resize(d["capacity"])
        self.node_count = d["node_count"]

        cdef int* children_left = <int*> (<np.ndarray> d["children_left"]).data
        cdef int* children_right =  <int*> (<np.ndarray> d["children_right"]).data
        cdef int* feature = <int*> (<np.ndarray> d["feature"]).data
        cdef double* threshold = <double*> (<np.ndarray> d["threshold"]).data
        cdef double* value = <double*> (<np.ndarray> d["value"]).data
        cdef double* best_error = <double*> (<np.ndarray> d["best_error"]).data
        cdef double* init_error = <double*> (<np.ndarray> d["init_error"]).data
        cdef int* n_samples = <int*> (<np.ndarray> d["n_samples"]).data

        memcpy(self.children_left, children_left, self.capacity * sizeof(int))
        memcpy(self.children_right, children_right, self.capacity * sizeof(int))
        memcpy(self.feature, feature, self.capacity * sizeof(int))
        memcpy(self.threshold, threshold, self.capacity * sizeof(double))
        memcpy(self.value, value, self.capacity * self.value_stride * sizeof(double))
        memcpy(self.best_error, best_error, self.capacity * sizeof(double))
        memcpy(self.init_error, init_error, self.capacity * sizeof(double))
        memcpy(self.n_samples, n_samples, self.capacity * sizeof(int))

    cdef void resize(self, int capacity=-1):
        """Resize all inner arrays to `capacity`, if < 0 double capacity."""
        if capacity == self.capacity:
            return

        if capacity < 0:
            if self.capacity <= 0:
                capacity = 3 # default initial value
            else:
                capacity = 2 * self.capacity

        self.capacity = capacity

        cdef int* tmp_children_left = <int*> realloc(self.children_left, capacity * sizeof(int))
        if tmp_children_left != NULL: self.children_left = tmp_children_left

        cdef int* tmp_children_right = <int*> realloc(self.children_right, capacity * sizeof(int))
        if tmp_children_right != NULL: self.children_right = tmp_children_right

        cdef int* tmp_feature = <int*> realloc(self.feature, capacity * sizeof(int))
        if tmp_feature != NULL: self.feature = tmp_feature

        cdef double* tmp_threshold = <double*> realloc(self.threshold, capacity * sizeof(double))
        if tmp_threshold != NULL: self.threshold = tmp_threshold

        cdef double* tmp_value = <double*> realloc(self.value, capacity * self.value_stride * sizeof(double))
        if tmp_value != NULL: self.value = tmp_value

        cdef double* tmp_best_error = <double*> realloc(self.best_error, capacity * sizeof(double))
        if tmp_best_error != NULL: self.best_error = tmp_best_error

        cdef double* tmp_init_error = <double*> realloc(self.init_error, capacity * sizeof(double))
        if tmp_init_error != NULL: self.init_error = tmp_init_error

        cdef int* tmp_n_samples = <int*> realloc(self.n_samples, capacity * sizeof(int))
        if tmp_n_samples != NULL: self.n_samples = tmp_n_samples

        if tmp_children_left == NULL or \
           tmp_children_right == NULL or \
           tmp_feature == NULL or \
           tmp_threshold == NULL or \
           tmp_value == NULL or \
           tmp_best_error == NULL or \
           tmp_init_error == NULL or \
           tmp_n_samples == NULL:
            raise MemoryError()

        # if capacity smaller than node_count, adjust the counter
        if capacity < self.node_count:
            self.node_count = capacity

    cpdef build(self, np.ndarray X, np.ndarray y,
                np.ndarray sample_mask=None,
                np.ndarray X_argsorted=None,
                np.ndarray sample_weight=None):
        """Build a decision tree from the training set (X, y).

        Parameters
        ----------
        X : ndarray of shape [n_samples, n_features]
            The training input samples.

        y : ndarray of shape [n_samples, n_outputs]
            The target values.
        """
        # Check input before recursive partitioning
        if X.dtype != DTYPE or not np.isfortran(X):
            X = np.asarray(X, dtype=DTYPE, order="F")

        if y.dtype != DOUBLE or not y.flags.contiguous:
            y = np.asarray(y, dtype=DOUBLE, order="C")

        if sample_weight is not None:
            if sample_weight.dtype != DOUBLE or not sample_weight.flags.contiguous:
                sample_weight = np.asarray(
                        sample_weight, dtype=DOUBLE, order="C")

        if sample_mask is None:
            sample_mask = np_ones((X.shape[0],), dtype=np_bool)

        if X_argsorted is None:
            X_argsorted = np.asfortranarray(
                np.argsort(X.T, axis=1).astype(np.int32).T)

        # Pre-allocate some space
        cdef int init_capacity
        cdef int n_node_samples
        cdef double weighted_n_node_samples

        if self.max_depth <= 10:
            init_capacity = (2 ** (<int>(self.max_depth) + 1)) - 1
        else:
            init_capacity = 2047

        self.resize(init_capacity)
        cdef double* buffer_value = <double*> malloc(self.value_stride * sizeof(double))

        n_node_samples = np.sum(sample_mask)
        if sample_weight is not None:
            weighted_n_node_samples = np.sum(sample_weight[sample_mask])
        else:
            weighted_n_node_samples = n_node_samples

        # Build the tree by recursive partitioning
        self.recursive_partition(X,
                                 X_argsorted,
                                 y,
                                 sample_weight,
                                 sample_mask,
                                 n_node_samples,
                                 weighted_n_node_samples,
                                 0,
                                 -1,
                                 False,
                                 buffer_value)

        # Compactify
        self.resize(self.node_count)
        free(buffer_value)

    cdef void recursive_partition(self,
                                  np.ndarray[DTYPE_t, ndim=2, mode="fortran"] X,
                                  np.ndarray[np.int32_t, ndim=2, mode="fortran"] X_argsorted,
                                  np.ndarray[DOUBLE_t, ndim=2, mode="c"] y,
                                  np.ndarray[DOUBLE_t, ndim=1, mode="c"] sample_weight,
                                  np.ndarray sample_mask,
                                  int n_node_samples,
                                  double weighted_n_node_samples,
                                  int depth,
                                  int parent,
                                  int is_left_child,
                                  double* buffer_value) except *:
        """Recursive partition algorithm for the tree construction."""
        # Variables
        cdef Criterion criterion = self.criterion

        cdef DTYPE_t* X_ptr = <DTYPE_t*> X.data
        cdef int* X_argsorted_ptr = <int*> X_argsorted.data
        cdef DOUBLE_t* y_ptr = <DOUBLE_t*> y.data
        cdef BOOL_t* sample_mask_ptr = <BOOL_t*> sample_mask.data

        cdef DOUBLE_t* sample_weight_ptr = NULL
        if sample_weight is not None:
            sample_weight_ptr = <DOUBLE_t*> sample_weight.data
        cdef DOUBLE_t w = 1.0

        cdef Py_ssize_t X_stride = <Py_ssize_t> X.strides[1] / <int> X.itemsize
        cdef Py_ssize_t X_argsorted_stride = <Py_ssize_t> X_argsorted.strides[1] / <int> X_argsorted.itemsize
        cdef Py_ssize_t y_stride = <Py_ssize_t> y.strides[0] / <int> y.itemsize

        cdef int n_total_samples = y.shape[0]
        cdef int feature
        cdef double threshold
        cdef double best_error
        cdef double init_error

        cdef int i
        cdef np.ndarray sample_mask_left
        cdef np.ndarray sample_mask_right
        cdef BOOL_t* sample_mask_left_ptr = NULL
        cdef BOOL_t* sample_mask_right_ptr = NULL
        cdef int n_node_samples_left = 0
        cdef int n_node_samples_right = 0
        cdef double weighted_n_node_samples_left = 0.0
        cdef double weighted_n_node_samples_right = 0.0

        # Count samples
        if n_node_samples == 0:
            raise ValueError("Attempting to find a split "
                             "with an empty sample_mask.")

        if weighted_n_node_samples < 0.0:
            raise ValueError("Attempting to find a split with a negative "
                             "weighted number of samples.")

        # Split samples
        if depth < self.max_depth and \
           n_node_samples >= self.min_samples_split and \
           n_node_samples >= 2 * self.min_samples_leaf:
            self.find_split(X_ptr, X_stride,
                            X_argsorted_ptr, X_argsorted_stride,
                            y_ptr, y_stride,
                            sample_weight_ptr,
                            sample_mask_ptr,
                            n_node_samples,
                            weighted_n_node_samples,
                            n_total_samples,
                            &feature, &threshold, &best_error, &init_error)

        else:
            feature = -1
            criterion.init(y_ptr, y_stride,
                           sample_weight_ptr,
                           sample_mask_ptr,
                           n_node_samples,
                           weighted_n_node_samples,
                           n_total_samples)
            init_error = criterion.eval()

        criterion.init_value(buffer_value)

        # Current node is leaf
        if feature == -1:
            self.add_leaf(parent, is_left_child, buffer_value, init_error, n_node_samples)

        # Current node is internal node (= split node)
        else:
            # Sample mask is too sparse?
            if 1. * n_node_samples / n_total_samples <= self.min_density:
                X = X[sample_mask]
                X_argsorted = np.asfortranarray(np.argsort(X.T, axis=1).astype(np.int32).T)
                y = y[sample_mask]
                if sample_weight is not None:
                    sample_weight = sample_weight[sample_mask]
                    sample_weight_ptr = <DOUBLE_t*> sample_weight.data
                sample_mask = np_ones((n_node_samples, ), dtype=np_bool)

                n_total_samples = n_node_samples

                X_ptr = <DTYPE_t*> X.data
                X_stride = <Py_ssize_t> X.strides[1] / <int> X.itemsize
                sample_mask_ptr = <BOOL_t*> sample_mask.data

                # !! No need to update the other variables
                # X_argsorted_ptr = <int*> X_argsorted.data
                # y_ptr = <DOUBLE_t*> y.data
                # X_argsorted_stride = <Py_ssize_t> X_argsorted.strides[1] / <int> X_argsorted.itemsize
                # y_stride = <Py_ssize_t> y.strides[0] / <int> y.itemsize

            # Split
            X_ptr = X_ptr + feature * X_stride

            sample_mask_left = np_zeros((n_total_samples, ), dtype=np_bool)
            sample_mask_right = np_zeros((n_total_samples, ), dtype=np_bool)
            sample_mask_left_ptr = <BOOL_t*> sample_mask_left.data
            sample_mask_right_ptr = <BOOL_t*> sample_mask_right.data

            n_node_samples_left = 0
            n_node_samples_right = 0
            weighted_n_node_samples_left = 0.0
            weighted_n_node_samples_right = 0.0

            for i from 0 <= i < n_total_samples:
                if sample_mask_ptr[i]:
                    if sample_weight_ptr != NULL:
                        w = sample_weight_ptr[i]
                    if X_ptr[i] <= threshold:
                        sample_mask_left_ptr[i] = 1
                        n_node_samples_left += 1
                        weighted_n_node_samples_left += w
                    else:
                        sample_mask_right_ptr[i] = 1
                        n_node_samples_right += 1
                        weighted_n_node_samples_right += w

            # Make current node a leaf if no valid split was found
            if (weighted_n_node_samples_left <= 0 or
                weighted_n_node_samples_right <= 0):

                self.add_leaf(parent, is_left_child, buffer_value, init_error, n_node_samples)
                return

            node_id = self.add_split_node(parent, is_left_child, feature,
                                          threshold, buffer_value, best_error,
                                          init_error, n_node_samples)

            # Left child recursion
            self.recursive_partition(X, X_argsorted,
                                     y, sample_weight,
                                     sample_mask_left,
                                     n_node_samples_left,
                                     weighted_n_node_samples_left,
                                     depth + 1, node_id,
                                     True, buffer_value)

            # Right child recursion
            self.recursive_partition(X, X_argsorted,
                                     y, sample_weight,
                                     sample_mask_right,
                                     n_node_samples_right,
                                     weighted_n_node_samples_right,
                                     depth + 1, node_id,
                                     False, buffer_value)

    cdef int add_split_node(self, int parent, int is_left_child, int feature,
                                  double threshold, double* value,
                                  double best_error, double init_error,
                                  int n_samples):
        """Add a splitting node to the tree. The new node registers itself as
           the child of its parent. """
        cdef int node_id = self.node_count

        if node_id >= self.capacity:
            self.resize()

        self.feature[node_id] = feature
        self.threshold[node_id] = threshold

        cdef int offset_node = node_id * self.value_stride
        memcpy(self.value + offset_node, value, self.value_stride * sizeof(double))

        self.init_error[node_id] = init_error
        self.best_error[node_id] = best_error
        self.n_samples[node_id] = n_samples

        # set as left or right child of parent
        if parent > _TREE_LEAF:
            if is_left_child:
                self.children_left[parent] = node_id
            else:
                self.children_right[parent] = node_id

        self.node_count += 1

        return node_id

    cdef int add_leaf(self, int parent, int is_left_child, double* value,
                      double error, int n_samples):
        """Add a leaf to the tree. The new node registers itself as the
           child of its parent. """
        cdef int node_id = self.node_count

        if node_id >= self.capacity:
            self.resize()

        cdef int offset_node = node_id * self.n_outputs * self.max_n_classes
        memcpy(self.value + offset_node, value, self.value_stride * sizeof(double))

        self.init_error[node_id] = error
        self.best_error[node_id] = error
        self.n_samples[node_id] = n_samples

        if parent >= 0:
            if is_left_child:
                self.children_left[parent] = node_id
            else:
                self.children_right[parent] = node_id

        self.children_left[node_id] = _TREE_LEAF
        self.children_right[node_id] = _TREE_LEAF

        self.node_count += 1

        return node_id

    cdef void find_split(self, DTYPE_t* X_ptr, Py_ssize_t X_stride,
                         int* X_argsorted_ptr, Py_ssize_t X_argsorted_stride,
                         DOUBLE_t* y_ptr, Py_ssize_t y_stride,
                         DOUBLE_t* sample_weight_ptr,
                         BOOL_t* sample_mask_ptr,
                         int n_node_samples,
                         double weighted_n_node_samples,
                         int n_total_samples,
                         int* _best_i,
                         double* _best_t, double* _best_error,
                         double* _initial_error):
        """Find the best dimension and threshold that minimises the error."""
        if self.find_split_algorithm == _TREE_SPLIT_BEST:
            self.find_best_split(X_ptr, X_stride,
                                 X_argsorted_ptr, X_argsorted_stride,
                                 y_ptr, y_stride,
                                 sample_weight_ptr,
                                 sample_mask_ptr,
                                 n_node_samples,
                                 weighted_n_node_samples,
                                 n_total_samples, _best_i, _best_t,
                                 _best_error, _initial_error)

        elif self.find_split_algorithm == _TREE_SPLIT_RANDOM:
            self.find_random_split(X_ptr, X_stride,
                                   X_argsorted_ptr, X_argsorted_stride,
                                   y_ptr, y_stride,
                                   sample_weight_ptr,
                                   sample_mask_ptr,
                                   n_node_samples,
                                   weighted_n_node_samples,
                                   n_total_samples, _best_i, _best_t,
                                   _best_error, _initial_error)

    cdef void find_best_split(self, DTYPE_t* X_ptr, Py_ssize_t X_stride,
                              int* X_argsorted_ptr, Py_ssize_t X_argsorted_stride,
                              DOUBLE_t* y_ptr, Py_ssize_t y_stride,
                              DOUBLE_t* sample_weight_ptr,
                              BOOL_t* sample_mask_ptr,
                              int n_node_samples,
                              double weighted_n_node_samples,
                              int n_total_samples, int* _best_i,
                              double* _best_t, double* _best_error,
                              double* _initial_error):
        """Implementation of `find_split` that looks for the best threshold."""
        # Variables declarations
        cdef Criterion criterion = self.criterion
        cdef int n_features = self.n_features
        cdef int max_features = self.max_features
        cdef int visited_features = 0
        cdef int min_samples_leaf = self.min_samples_leaf
        cdef object random_state = self.random_state

        cdef int i, a, b, best_i = -1
        cdef np.int32_t feature_idx = -1
        cdef int n_left = 0

        cdef double t, initial_error, error
        cdef double best_error = INFINITY, best_t = INFINITY

        cdef DTYPE_t* X_i = NULL
        cdef int* X_argsorted_i = NULL
        cdef DTYPE_t X_a, X_b

        cdef np.ndarray[np.int32_t, ndim=1, mode="c"] features = self.features

        # Compute the initial criterion value in the node
        criterion.init(y_ptr, y_stride,
                       sample_weight_ptr,
                       sample_mask_ptr,
                       n_node_samples,
                       weighted_n_node_samples,
                       n_total_samples)
        initial_error = criterion.eval()

        if initial_error == 0:  # break early if the node is pure
            _best_i[0] = best_i
            _best_t[0] = best_t
            _best_error[0] = initial_error
            _initial_error[0] = initial_error

            return

        # Features to consider
        if max_features < 0 or max_features >= n_features:
            max_features = n_features
        else:
            random_state.shuffle(features)

        # Look for the best split
        for feature_idx from 0 <= feature_idx < n_features:
            i = features[feature_idx]

            # Get i-th col of X and X_sorted
            X_i = X_ptr + X_stride * i
            X_argsorted_i = X_argsorted_ptr + X_argsorted_stride * i

            # Reset the criterion for this feature
            criterion.reset()

            # Index of smallest sample in X_argsorted_i that is in the sample mask
            a = 0

            while sample_mask_ptr[X_argsorted_i[a]] == 0:
                a = a + 1

            # Check that the feature is not constant
            b = _smallest_sample_larger_than(a, X_i, X_argsorted_i,
                                             sample_mask_ptr, n_total_samples)

            if b == -1:
                continue # Skip that feature and don't count it as visited

            # Consider splits between two consecutive samples
            while True:
                # Find the following larger sample
                b = _smallest_sample_larger_than(a, X_i, X_argsorted_i,
                                                 sample_mask_ptr, n_total_samples)
                if b == -1:
                    break

                # Better split than the best so far?
                if not criterion.update(a, b,
                                        y_ptr, y_stride,
                                        X_argsorted_i,
                                        sample_weight_ptr,
                                        sample_mask_ptr):
                    a = b
                    continue

                # Only consider splits that respect min_leaf
                n_left = criterion.n_left
                if (n_left < min_samples_leaf or
                    (n_node_samples - n_left) < min_samples_leaf):
                    a = b
                    continue

                error = criterion.eval()

                if error < best_error:
                    X_a = X_i[X_argsorted_i[a]]
                    X_b = X_i[X_argsorted_i[b]]

                    t = X_a + (X_b - X_a) / 2.0
                    if t == X_b:
                        t = X_a

                    best_i = i
                    best_t = t
                    best_error = error

                # Proceed to the next interval
                a = b

            # Count one more visited feature
            visited_features += 1

            if visited_features >= max_features:
                break

        _best_i[0] = best_i
        _best_t[0] = best_t
        _best_error[0] = best_error
        _initial_error[0] = initial_error

    cdef void find_random_split(self, DTYPE_t* X_ptr, Py_ssize_t X_stride,
                                int* X_argsorted_ptr, Py_ssize_t X_argsorted_stride,
                                DOUBLE_t* y_ptr, Py_ssize_t y_stride,
                                DOUBLE_t* sample_weight_ptr,
                                BOOL_t* sample_mask_ptr,
                                int n_node_samples,
                                double weighted_n_node_samples,
                                int n_total_samples, int* _best_i,
                                double* _best_t, double* _best_error,
                                double* _initial_error):
        """Implementation of `find_split` that looks for the best threshold
           among randomly drawn thresholds at each feature."""
        # Variables declarations
        cdef Criterion criterion = self.criterion
        cdef int n_features = self.n_features
        cdef int max_features = self.max_features
        cdef int visited_features = 0
        cdef int min_samples_leaf = self.min_samples_leaf
        cdef object random_state = self.random_state

        cdef int i, a, b, c, best_i = -1
        cdef np.int32_t feature_idx = -1
        cdef int n_left = 0
        cdef double random

        cdef double t, initial_error, error
        cdef double best_error = INFINITY, best_t = INFINITY

        cdef DTYPE_t* X_i = NULL
        cdef int* X_argsorted_i = NULL
        cdef DTYPE_t X_a, X_b

        cdef np.ndarray[np.int32_t, ndim=1, mode="c"] features = self.features

        # Compute the initial criterion value in the node
        criterion.init(y_ptr, y_stride,
                       sample_weight_ptr,
                       sample_mask_ptr,
                       n_node_samples,
                       weighted_n_node_samples,
                       n_total_samples)
        initial_error = criterion.eval()

        if initial_error == 0:  # break early if the node is pure
            _best_i[0] = best_i
            _best_t[0] = best_t
            _best_error[0] = initial_error
            _initial_error[0] = initial_error

            return

        # Features to consider
        if max_features < 0 or max_features >= n_features:
            max_features = n_features
        else:
            random_state.shuffle(features)

        # Look for the best split
        for feature_idx from 0 <= feature_idx < n_features:
            i = features[feature_idx]

            # Get i-th col of X and X_sorted
            X_i = X_ptr + X_stride * i
            X_argsorted_i = X_argsorted_ptr + X_argsorted_stride * i

            # Reset the criterion for this feature
            criterion.reset()

            # Find min and max
            a = 0
            while sample_mask_ptr[X_argsorted_i[a]] == 0:
                a = a + 1
            X_a = X_i[X_argsorted_i[a]]

            b = n_total_samples - 1
            while sample_mask_ptr[X_argsorted_i[b]] == 0:
                b = b - 1
            X_b = X_i[X_argsorted_i[b]]

            if b <= a or X_a == X_b:
                continue # Skip that feature and don't count it as visited

            # Draw a random threshold in [a, b)
            random = random_state.rand()
            t = X_a + (random * (X_b - X_a))
            if t == X_b:
                t = X_a

            # Find the sample just greater than t
            c = a + 1

            while True:
                if sample_mask_ptr[X_argsorted_i[c]] != 0:
                    # FIXME why is t cast to DTYPE_t?
                    if X_i[X_argsorted_i[c]] > (<DTYPE_t> t) or c == b:
                        break

                c += 1

            # Better than the best so far?
            if not criterion.update(0, c,
                                    y_ptr, y_stride,
                                    X_argsorted_i,
                                    sample_weight_ptr,
                                    sample_mask_ptr):
                continue

            n_left = criterion.n_left

            if (n_left < min_samples_leaf or
                (n_node_samples - n_left) < min_samples_leaf):
                continue

            error = criterion.eval()

            if error < best_error:
                best_i = i
                best_t = t
                best_error = error

            # Count one more visited feature
            visited_features += 1

            if visited_features >= max_features:
                break

        _best_i[0] = best_i
        _best_t[0] = best_t
        _best_error[0] = best_error
        _initial_error[0] = initial_error

    cpdef predict(self, np.ndarray[DTYPE_t, ndim=2] X):
        """Predict target for X."""
        cdef int i, k, c
        cdef int n_samples = X.shape[0]
        cdef int node_id = 0
        cdef int offset_node
        cdef int offset_output

        cdef np.ndarray[np.float64_t, ndim=3] out
        out = np_zeros((n_samples, self.n_outputs, self.max_n_classes), dtype=np.float64)

        for i from 0 <= i < n_samples:
            node_id = 0

            # While node_id not a leaf
            while self.children_left[node_id] != _TREE_LEAF: # and self.children_right[node_id] != _TREE_LEAF:
                if X[i, self.feature[node_id]] <= self.threshold[node_id]:
                    node_id = self.children_left[node_id]
                else:
                    node_id = self.children_right[node_id]

            offset_node = node_id * self.value_stride

            for k from 0 <= k < self.n_outputs:
                offset_output = k * self.max_n_classes

                for c from 0 <= c < self.n_classes[k]:
                    out[i, k, c] = self.value[offset_node + offset_output + c]

        return out

    cpdef apply(self, np.ndarray[DTYPE_t, ndim=2] X):
        """Finds the terminal region (=leaf node) for each sample in X."""
        cdef int i = 0
        cdef int n_samples = X.shape[0]
        cdef int node_id = 0

        cdef np.ndarray[np.int32_t, ndim=1] out
        out = np_zeros((n_samples, ), dtype=np.int32)

        for i from 0 <= i < n_samples:
            node_id = 0

            # While node_id not a leaf
            while self.children_left[node_id] != _TREE_LEAF: # and self.children_right[node_id] != _TREE_LEAF:
                if X[i, self.feature[node_id]] <= self.threshold[node_id]:
                    node_id = self.children_left[node_id]
                else:
                    node_id = self.children_right[node_id]

            out[i] = node_id

        return out

    cpdef compute_feature_importances(self):
        """Computes the importance of each feature (aka variable)."""
        cdef int node
        cdef np.ndarray[np.float64_t, ndim=1] importances
        importances = np_zeros((self.n_features,), dtype=np.float64)

        for node from 0 <= node < self.node_count:
            if self.children_left[node] != _TREE_LEAF: # and self.children_right[node] != _TREE_LEAF:
                importances[self.feature[node]] += \
                    self.n_samples[node] * (self.init_error[node] - self.best_error[node])

        cdef double normalizer = np.sum(importances)

        if normalizer > 0.0:
            # Avoid dividing by zero (e.g., when root is pure)
            importances /= normalizer

        return importances


# =============================================================================
# Criterion
# =============================================================================

cdef class Criterion:
    """Interface for splitting criteria (regression and classification)."""

    cdef void init(self, DOUBLE_t* y, Py_ssize_t y_stride,
                         DOUBLE_t* sample_weight,
                         BOOL_t* sample_mask,
                         int n_samples,
                         double weighted_n_samples,
                         int n_total_samples):
        """Initialise the criterion."""
        pass

    cdef void reset(self):
        """Reset the criterion for a new feature index."""
        pass

    cdef bool update(self, int a, int b,
                           DOUBLE_t* y, Py_ssize_t y_stride,
                           int* X_argsorted_i,
                           DOUBLE_t* sample_weight,
                           BOOL_t* sample_mask):
        """Update the criteria for each value in interval [a,b) (where a and b
           are indices in `X_argsorted_i`)."""
        pass

    cdef double eval(self):
        """Evaluate the criteria (aka the split error)."""
        pass

    cdef void init_value(self, double* buffer_value):
        """Get the initial value of the criterion (`init` must be called
           before)."""
        pass


cdef class ClassificationCriterion(Criterion):
    """Abstract criterion for classification.

    Attributes
    ----------
    n_outputs : int
        The number of outputs.

    n_classes : int*
        n_classes[k] is the number of classes for output k.

    n_samples : int
        The number of samples.

    weighted_n_samples : double
        The weighted number of samples.

    label_count_stride : Py_ssize_t
        The stride between outputs in label_count_* arrays.

    label_count_left : double*
        label_count_left[k * label_count_stride + c] is the number of samples
        of class c left of splitting point for output k.

    label_count_right : double*
        label_count_rightt[k * label_count_stride + c] is the number of samples
        of class c right of splitting point for output k.

    label_count_init : double*
        label_count_init[k * label_count_stride + c] is the initial number of
        samples of class c for output k. Used to reset `label_count_right` for
        each feature.

    n_left : int
        The number of samples left of splitting point.

    n_right : int
        The number of samples right of splitting point.

    weighted_n_left : double
        The weighted number of samples left of splitting point.

    weighted_n_right : double
        The weighted number of samples right of splitting point.

    References
    ----------

    [1] Hastie et al. "Elements of Statistical Learning", 2009.
    """
    cdef int* n_classes

    cdef Py_ssize_t label_count_stride
    cdef double* label_count_left
    cdef double* label_count_right
    cdef double* label_count_init

    def __cinit__(self, int n_outputs, object n_classes):
        """Constructor."""
        cdef int k = 0

        self.n_outputs = n_outputs
        self.n_samples = 0
        self.weighted_n_samples = 0.0
        self.n_left = 0
        self.n_right = 0
        self.weighted_n_left = 0.0
        self.weighted_n_right = 0.0

        self.n_classes = <int*> malloc(n_outputs * sizeof(int))
        if self.n_classes == NULL:
            raise MemoryError()

        cdef Py_ssize_t label_count_stride = -1

        for k from 0 <= k < n_outputs:
            self.n_classes[k] = n_classes[k]

            if n_classes[k] > label_count_stride:
                label_count_stride = n_classes[k]

        self.label_count_stride = label_count_stride

        # Allocate
        self.label_count_left = <double*> calloc(n_outputs * label_count_stride, sizeof(double))
        self.label_count_right = <double*> calloc(n_outputs * label_count_stride, sizeof(double))
        self.label_count_init = <double*> calloc(n_outputs * label_count_stride, sizeof(double))

        # Check for allocation errors
        if self.label_count_left == NULL or \
           self.label_count_right == NULL or \
           self.label_count_init == NULL:
            free(self.n_classes)
            free(self.label_count_left)
            free(self.label_count_right)
            free(self.label_count_init)
            raise MemoryError()

    def __dealloc__(self):
        """Destructor."""
        free(self.n_classes)
        free(self.label_count_left)
        free(self.label_count_right)
        free(self.label_count_init)

    def __reduce__(self):
        return (ClassificationCriterion,
                (self.n_outputs, intp_to_ndarray(self.n_classes,
                                                 self.n_outputs)),
                self.__getstate__())

    def __getstate__(self):
        return {}

    def __setstate__(self, d):
        pass

    cdef void init(self, DOUBLE_t* y, Py_ssize_t y_stride,
                         DOUBLE_t* sample_weight,
                         BOOL_t* sample_mask,
                         int n_samples,
                         double weighted_n_samples,
                         int n_total_samples):
        """Initialise the criterion."""
        cdef int n_outputs = self.n_outputs
        cdef int* n_classes = self.n_classes
        cdef Py_ssize_t label_count_stride = self.label_count_stride
        cdef double* label_count_init = self.label_count_init

        cdef int k = 0
        cdef int c = 0
        cdef int j = 0
        cdef DOUBLE_t w = 1.0

        self.n_samples = n_samples
        self.weighted_n_samples = weighted_n_samples

        for k from 0 <= k < n_outputs:
            for c from 0 <= c < n_classes[k]:
                label_count_init[k * label_count_stride + c] = 0

        for j from 0 <= j < n_total_samples:
            if sample_mask[j] == 0:
                continue
            if sample_weight != NULL:
                w = sample_weight[j]

            for k from 0 <= k < n_outputs:
                c = <int>y[j * y_stride + k]
                label_count_init[k * label_count_stride + c] += w

        self.reset()

    cdef void reset(self):
        """Reset the criterion for a new feature index."""
        cdef int n_outputs = self.n_outputs
        cdef int* n_classes = self.n_classes
        cdef Py_ssize_t label_count_stride = self.label_count_stride
        cdef double* label_count_init = self.label_count_init
        cdef double* label_count_left = self.label_count_left
        cdef double* label_count_right = self.label_count_right

        cdef int k = 0
        cdef int c = 0
        self.n_left = 0
        self.n_right = self.n_samples
        self.weighted_n_left = 0.0
        self.weighted_n_right = self.weighted_n_samples

        for k from 0 <= k < n_outputs:
            for c from 0 <= c < n_classes[k]:
                # Reset left label counts to 0
                label_count_left[k * label_count_stride + c] = 0

                # Reset right label counts to the initial counts
                label_count_right[k * label_count_stride + c] = label_count_init[k * label_count_stride + c]

    cdef bool update(self, int a, int b,
                           DOUBLE_t* y, Py_ssize_t y_stride,
                           int* X_argsorted_i,
                           DOUBLE_t* sample_weight,
                           BOOL_t* sample_mask):
        """Update the criteria for each value in interval [a,b) (where a and b
           are indices in `X_argsorted_i`)."""
        cdef int n_outputs = self.n_outputs
        cdef int* n_classes = self.n_classes
        cdef Py_ssize_t label_count_stride = self.label_count_stride
        cdef double* label_count_left = self.label_count_left
        cdef double* label_count_right = self.label_count_right
        cdef int n_left = self.n_left
        cdef int n_right = self.n_right
        cdef double weighted_n_left = self.weighted_n_left
        cdef double weighted_n_right = self.weighted_n_right

        cdef int idx, k, c, s
        cdef DOUBLE_t w = 1.

        # post condition: all samples from [0:b) are on the left side
        for idx from a <= idx < b:
            s = X_argsorted_i[idx]

            if sample_mask[s] == 0:
                continue
            if sample_weight != NULL:
                w = sample_weight[s]

            for k from 0 <= k < n_outputs:
                c = <int>y[s * y_stride + k]
                label_count_left[k * label_count_stride + c] += w
                label_count_right[k * label_count_stride + c] -= w

            n_left += 1
            n_right -= 1
            weighted_n_left += w
            weighted_n_right -= w

        self.n_left = n_left
        self.n_right = n_right
        self.weighted_n_left = weighted_n_left
        self.weighted_n_right = weighted_n_right

        # Skip splits that result in nodes with net 0 or negative weight
        if (weighted_n_left <= 0 or
            (self.weighted_n_samples - weighted_n_left) <= 0):
            return False

        # Prevent any single class from having a net negative weight
        for k from 0 <= k < n_outputs:
            for c from 0 <= c < n_classes[k]:
                if (label_count_left[k * label_count_stride + c] < 0 or
                    label_count_right[k * label_count_stride + c] < 0):
                    return False

        return True

    cdef double eval(self):
        """Evaluate the criteria (aka the split error)."""
        pass

    cdef void init_value(self, double* buffer_value):
        """Get the initial value of the criterion (`init` must be called
           before)."""
        cdef int n_outputs = self.n_outputs
        cdef int* n_classes = self.n_classes
        cdef Py_ssize_t label_count_stride = self.label_count_stride
        cdef double* label_count_init = self.label_count_init

        cdef int k, c

        for k from 0 <= k < n_outputs:
            for c from 0 <= c < n_classes[k]:
                buffer_value[k * label_count_stride + c] = (
                    label_count_init[k * label_count_stride + c])


cdef class Gini(ClassificationCriterion):
    """Gini Index splitting criteria.

    Let the target be a classification outcome taking values in 0, 1, ..., K-1.
    If node m represents a region Rm with Nm observations, then let

        pmk = 1/ Nm \sum_{x_i in Rm} I(yi = k)

    be the proportion of class k observations in node m.

    The Gini Index is then defined as:

        index = \sum_{k=0}^{K-1} pmk (1 - pmk)
              = 1 - \sum_{k=0}^{K-1} pmk ** 2
    """

    cdef double eval(self):
        """Returns Gini index of left branch + Gini index of right branch."""
        cdef double n_samples = self.weighted_n_samples
        cdef int n_outputs = self.n_outputs
        cdef int* n_classes = self.n_classes
        cdef Py_ssize_t label_count_stride = self.label_count_stride
        cdef double* label_count_left = self.label_count_left
        cdef double* label_count_right = self.label_count_right
        cdef double n_left = self.weighted_n_left
        cdef double n_right = self.weighted_n_right

        cdef double total_left = 0.0
        cdef double total_right = 0.0
        cdef double H_left
        cdef double H_right
        cdef int k, c
        cdef double count_left, count_right

        for k from 0 <= k < n_outputs:
            H_left = n_left * n_left
            H_right = n_right * n_right

            for c from 0 <= c < n_classes[k]:
                count_left = label_count_left[k * label_count_stride + c]
                if count_left > 0:
                    H_left -= (count_left * count_left)

                count_right = label_count_right[k * label_count_stride + c]
                if count_right > 0:
                    H_right -= (count_right * count_right)

            if n_left == 0:
                H_left = 0
            else:
                H_left /= n_left

            if n_right == 0:
                H_right = 0
            else:
                H_right /= n_right

            total_left += H_left
            total_right += H_right

        return (total_left + total_right) / (n_samples * n_outputs)


cdef class Entropy(ClassificationCriterion):
    """Cross Entropy splitting criteria.

    Let the target be a classification outcome taking values in 0, 1, ..., K-1.
    If node m represents a region Rm with Nm observations, then let

        pmk = 1/ Nm \sum_{x_i in Rm} I(yi = k)

    be the proportion of class k observations in node m.

    The cross-entropy is then defined as

        cross-entropy = - \sum_{k=0}^{K-1} pmk log(pmk)
    """

    cdef double eval(self):
        """Returns Entropy of left branch + Entropy index of right branch. """
        cdef double n_samples = self.weighted_n_samples
        cdef int n_outputs = self.n_outputs
        cdef int* n_classes = self.n_classes
        cdef Py_ssize_t label_count_stride = self.label_count_stride
        cdef double* label_count_left = self.label_count_left
        cdef double* label_count_right = self.label_count_right
        cdef double n_left = self.weighted_n_left
        cdef double n_right = self.weighted_n_right

        cdef double total = 0.0
        cdef double H_left
        cdef double H_right
        cdef int k, c
        cdef double e1, e2

        for k from 0 <= k < n_outputs:
            H_left = 0.0
            H_right = 0.0

            for c from 0 <= c < n_classes[k]:
                if label_count_left[k * label_count_stride + c] > 0:
                    H_left -= ((label_count_left[k * label_count_stride + c] / n_left) * log(label_count_left[k * label_count_stride + c] / n_left))

                if self.label_count_right[k * label_count_stride + c] > 0:
                    H_right -= ((label_count_right[k * label_count_stride + c] / n_right) * log(label_count_right[k * label_count_stride + c] / n_right))

            e1 = (n_left / n_samples) * H_left
            e2 = (n_right / n_samples) * H_right

            total += e1 + e2

        return total / n_outputs


cdef class RegressionCriterion(Criterion):
    """Abstract criterion for regression.

    Computes variance of the target values left and right of the split point.
    Computation is linear in `n_samples` by using ::

        var = \sum_i^n (y_i - y_bar) ** 2
            = (\sum_i^n y_i ** 2) - n_samples y_bar ** 2

    Attributes
    ----------
    n_outputs : int
        The number of outputs.

    n_samples : int
        The number of samples

    weighted_n_samples : double
        The weighted number of samples.

    mean_left : double*
        mean_left[k] is the mean target value of the samples left of the split
        point for output k.

    mean_right : double*
        mean_right[k] is the mean target value of the samples right of the split
        point for output k.

    sq_sum_left : double*
        sq_sum_left[k] is the sum of squared target values left of the split
        point for output k.

    sq_sum_right : double*
        sq_sum_right[k] is the sum of squared target values right of the split
        point for output k.

    var_left : double*
        var_left[k] is the variance of the values left of the split point for
        output k.

    var_right : double*
        var_right[k] is the variance of the values riht of the split point for
        output k.

    n_left : int
        The number of samples left of split point.

    n_right : int
        The number of samples right of split point.

    weighted_n_left : double
        The weighted number of samples left of splitting point.

    weighted_n_right : double
        The weighted number of samples right of splitting point.
    """
    cdef double* mean_left
    cdef double* mean_right
    cdef double* mean_init
    cdef double* sq_sum_left
    cdef double* sq_sum_right
    cdef double* sq_sum_init
    cdef double* var_left
    cdef double* var_right

    def __cinit__(self, int n_outputs):
        """Constructor."""
        cdef int k = 0

        self.n_outputs = n_outputs

        self.n_samples = 0
        self.weighted_n_samples = 0.0
        self.n_left = 0
        self.n_right = 0
        self.weighted_n_left = 0.0
        self.weighted_n_right = 0.0

        # Allocate
        self.mean_left = <double*> calloc(n_outputs, sizeof(double))
        self.mean_right = <double*> calloc(n_outputs, sizeof(double))
        self.mean_init = <double*> calloc(n_outputs, sizeof(double))
        self.sq_sum_left = <double*> calloc(n_outputs, sizeof(double))
        self.sq_sum_right = <double*> calloc(n_outputs, sizeof(double))
        self.sq_sum_init = <double*> calloc(n_outputs, sizeof(double))
        self.var_left = <double*> calloc(n_outputs, sizeof(double))
        self.var_right = <double*> calloc(n_outputs, sizeof(double))

        # Check for allocation errors
        if self.mean_left == NULL or \
           self.mean_right == NULL or \
           self.mean_init == NULL or \
           self.sq_sum_left == NULL or \
           self.sq_sum_right == NULL or \
           self.sq_sum_init == NULL or \
           self.var_left == NULL or \
           self.var_right == NULL:
            free(self.mean_left)
            free(self.mean_right)
            free(self.mean_init)
            free(self.sq_sum_left)
            free(self.sq_sum_right)
            free(self.sq_sum_init)
            free(self.var_left)
            free(self.var_right)
            raise MemoryError()

    def __dealloc__(self):
        """Destructor."""
        free(self.mean_left)
        free(self.mean_right)
        free(self.mean_init)
        free(self.sq_sum_left)
        free(self.sq_sum_right)
        free(self.sq_sum_init)
        free(self.var_left)
        free(self.var_right)

    def __reduce__(self):
        return (RegressionCriterion,
                (self.n_outputs,),
                self.__getstate__())

    def __getstate__(self):
        return {}

    def __setstate__(self, d):
        pass

    cdef void init(self, DOUBLE_t* y, Py_ssize_t y_stride,
                         DOUBLE_t* sample_weight,
                         BOOL_t* sample_mask,
                         int n_samples,
                         double weighted_n_samples,
                         int n_total_samples):
        """Initialise the criterion class; assume all samples
           are in the right branch and store the mean and squared
           sum in `self.mean_init` and `self.sq_sum_init`. """
        cdef double* mean_left = self.mean_left
        cdef double* mean_right = self.mean_right
        cdef double* mean_init = self.mean_init
        cdef double* sq_sum_left = self.sq_sum_left
        cdef double* sq_sum_right = self.sq_sum_right
        cdef double* sq_sum_init = self.sq_sum_init
        cdef double* var_left = self.var_left
        cdef double* var_right = self.var_right
        cdef int n_outputs = self.n_outputs

        cdef int k = 0

        for k from 0 <= k < n_outputs:
            mean_left[k] = 0.0
            mean_right[k] = 0.0
            mean_init[k] = 0.0
            sq_sum_right[k] = 0.0
            sq_sum_left[k] = 0.0
            sq_sum_init[k] = 0.0
            var_left[k] = 0.0
            var_right[k] = 0.0

        self.n_samples = n_samples
        self.weighted_n_samples = weighted_n_samples

        cdef DOUBLE_t w = 1.0
        cdef DOUBLE_t y_jk = 0.0
        cdef int j = 0

        for j from 0 <= j < n_total_samples:
            if sample_mask[j] == 0:
                continue
            if sample_weight != NULL:
                w = sample_weight[j]

            for k from 0 <= k < n_outputs:
                y_jk = y[j * y_stride + k]
                sq_sum_init[k] += w * y_jk * y_jk
                mean_init[k] += w * y_jk

        for k from 0 <= k < n_outputs:
            mean_init[k] /= weighted_n_samples

        self.reset()

    cdef void reset(self):
        """Reset criterion for new feature.

        Assume all data in right branch and copy statistics of the
        whole dataset into the auxiliary variables of the
        right branch.
        """
        cdef double* mean_left = self.mean_left
        cdef double* mean_right = self.mean_right
        cdef double* mean_init = self.mean_init
        cdef double* sq_sum_left = self.sq_sum_left
        cdef double* sq_sum_right = self.sq_sum_right
        cdef double* sq_sum_init = self.sq_sum_init
        cdef double* var_left = self.var_left
        cdef double* var_right = self.var_right

        cdef double weighted_n_samples = self.weighted_n_samples
        cdef int n_outputs = self.n_outputs

        cdef int k = 0

        self.n_right = self.n_samples
        self.n_left = 0
        self.weighted_n_right = self.weighted_n_samples
        self.weighted_n_left = 0.0

        for k from 0 <= k < n_outputs:
            mean_right[k] = mean_init[k]
            mean_left[k] = 0.0
            sq_sum_right[k] = sq_sum_init[k]
            sq_sum_left[k] = 0.0
            var_left[k] = 0.0
            var_right[k] = (sq_sum_right[k] -
                weighted_n_samples * (mean_right[k] * mean_right[k]))

    cdef bool update(self, int a, int b,
                           DOUBLE_t* y, Py_ssize_t y_stride,
                           int* X_argsorted_i,
                           DOUBLE_t* sample_weight,
                           BOOL_t* sample_mask):
        """Update the criteria for each value in interval [a,b) (where a and b
           are indices in `X_argsorted_i`)."""
        cdef double* mean_left = self.mean_left
        cdef double* mean_right = self.mean_right
        cdef double* sq_sum_left = self.sq_sum_left
        cdef double* sq_sum_right = self.sq_sum_right
        cdef double* var_left = self.var_left
        cdef double* var_right = self.var_right

        cdef int n_samples = self.n_samples
        cdef double weighted_n_samples = self.weighted_n_samples
        cdef int n_outputs = self.n_outputs
        cdef int n_left = self.n_left
        cdef int n_right = self.n_right
        cdef double weighted_n_left = self.weighted_n_left
        cdef double weighted_n_right = self.weighted_n_right

        cdef DOUBLE_t w = 1.0
        cdef DOUBLE_t y_idx = 0.0
        cdef int idx, j, k

        # post condition: all samples from [0:b) are on the left side
        for idx from a <= idx < b:
            j = X_argsorted_i[idx]

            if sample_mask[j] == 0:
                continue
            if sample_weight != NULL:
                w = sample_weight[j]

            for k from 0 <= k < n_outputs:
                y_idx = y[j * y_stride + k]
                sq_sum_left[k] += w * (y_idx * y_idx)
                sq_sum_right[k] -= w * (y_idx * y_idx)

                mean_left[k] = ((weighted_n_left * mean_left[k] + w * y_idx) /
                                (weighted_n_left + w))
                mean_right[k] = (((weighted_n_samples - weighted_n_left) *
                                      mean_right[k] - w * y_idx) /
                                 (weighted_n_samples - weighted_n_left - w))

            n_left += 1
            self.n_left = n_left
            n_right -= 1
            self.n_right = n_right
            weighted_n_left += w
            self.weighted_n_left = weighted_n_left
            weighted_n_right -= w
            self.weighted_n_right = weighted_n_right

            for k from 0 <= k < n_outputs:
                var_left[k] = sq_sum_left[k] - weighted_n_left * (mean_left[k] * mean_left[k])
                var_right[k] = sq_sum_right[k] - weighted_n_right * (mean_right[k] * mean_right[k])

        # Skip splits that result in nodes with net 0 or negative weight
        if (weighted_n_left <= 0 or
            (self.weighted_n_samples - weighted_n_left) <= 0):
            return False

        return True

    cdef double eval(self):
        """Evaluate the criteria (aka the split error)."""
        pass

    cdef void init_value(self, double* buffer_value):
        """Get the initial value of the criterion (`init` must be called
           before)."""
        cdef int n_outputs = self.n_outputs
        cdef double* mean_init = self.mean_init

        cdef int k

        for k from 0 <= k < n_outputs:
            buffer_value[k] = mean_init[k]


cdef class MSE(RegressionCriterion):
    """Mean squared error impurity criterion.

        MSE = var_left + var_right
    """

    cdef double eval(self):
        cdef double* var_left = self.var_left
        cdef double* var_right = self.var_right

        cdef int n_outputs = self.n_outputs

        cdef int k
        cdef double total = 0.0

        for k from 0 <= k < n_outputs:
            total += var_left[k]
            total += var_right[k]

        return total / n_outputs


# =============================================================================
# Utils
# =============================================================================

cdef inline np.ndarray intp_to_ndarray(int* data, int size):
    """Encapsulate data into a 1D numpy array of int's."""
    cdef np.npy_intp shape[1]
    shape[0] = <np.npy_intp> size
    return np.PyArray_SimpleNewFromData(1, shape, np.NPY_INT, data)

cdef inline np.ndarray doublep_to_ndarray(double* data, int size):
    """Encapsulate data into a 1D numpy array of double's."""
    cdef np.npy_intp shape[1]
    shape[0] = <np.npy_intp> size
    return np.PyArray_SimpleNewFromData(1, shape, np.NPY_DOUBLE, data)

cdef inline int _smallest_sample_larger_than(int sample_idx,
                                             DTYPE_t* X_i,
                                             int* X_argsorted_i,
                                             BOOL_t* sample_mask,
                                             int n_total_samples):
    """Find the largest next sample.

    Find the index in the `X_i` array for sample who's feature
    `i` value is just about greater than those of the sample
    `X_argsorted_i[sample_idx]`.

    Returns
    -------
    next_sample_idx : int
        The index of the next smallest sample in `X_argsorted`
        with different feature value than `sample_idx` .
        I.e. `X_argsorted_i[sample_idx] < X_argsorted_i[next_sample_idx]`
        -1 if no such element exists.
    """
    cdef int idx = 0, j
    cdef DTYPE_t threshold = -DBL_MAX

    if sample_idx > -1:
        threshold = X_i[X_argsorted_i[sample_idx]]

    for idx from sample_idx < idx < n_total_samples:
        j = X_argsorted_i[idx]

        if sample_mask[j] == 0:
            continue

        if X_i[j] > threshold + 1.e-7:
            return idx

    return -1

def _random_sample_mask(int n_total_samples, int n_total_in_bag, random_state):
    """Create a random sample mask where ``n_total_in_bag`` elements are set.

    Parameters
    ----------
    n_total_samples : int
        The length of the resulting mask.

    n_total_in_bag : int
        The number of elements in the sample mask which are set to 1.

    random_state : np.RandomState
        A numpy ``RandomState`` object.

    Returns
    -------
    sample_mask : np.ndarray, shape=[n_total_samples]
        An ndarray where ``n_total_in_bag`` elements are set to ``True``
        the others are ``False``.
    """
    cdef np.ndarray[np.float64_t, ndim=1, mode="c"] rand = \
         random_state.rand(n_total_samples)
    cdef np.ndarray[BOOL_t, ndim=1, mode="c"] sample_mask = \
         np_zeros((n_total_samples,), dtype=np.int8)

    cdef int n_bagged = 0
    cdef int i = 0

    for i from 0 <= i < n_total_samples:
        if rand[i] * (n_total_samples - i) < (n_total_in_bag - n_bagged):
            sample_mask[i] = 1
            n_bagged += 1

    return sample_mask.astype(np_bool)
