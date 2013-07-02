# encoding: utf-8
# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False

# Author: Peter Prettenhofer, Brian Holt, Gilles Louppe
# Licence: BSD 3 clause

from libc.stdlib cimport calloc, free, malloc, realloc
from libc.string cimport memcpy

import numpy as np
cimport numpy as np

from cpython cimport bool


# =============================================================================
# Criterion
# =============================================================================

cdef class Criterion:
    """Interface for impurity criteria."""

    cdef void init(self, DOUBLE_t* y,
                         SIZE_t y_stride,
                         DOUBLE_t* sample_weight,
                         SIZE_t* samples,
                         SIZE_t start,
                         SIZE_t i,
                         SIZE_t end):
        """Initialize the criterion at node samples[start:end] and
           children samples[start:i] and samples[i:end]."""
        pass

    cdef double node_impurity(self):
        """Evaluate the impurity of the current node, i.e. the impurity of
           samples[start:end]."""
        pass

    cdef double children_impurity(self):
        """Evaluate the impurity in children nodes, i.e. the impurity of
           samples[start:i] + the impurity of samples[i:end]."""
        pass

    cdef bool update(self, SIZE_t j):
        """Update the collected statistics by moving samples[i:j] from the
           right children to the left children."""
        pass


cdef class ClassificationCriterion(Criterion):
    """Abstract criterion for classification."""
    cdef SIZE_t* n_classes
    cdef SIZE_t label_count_stride
    cdef double* label_count_left
    cdef double* label_count_right
    cdef double* label_count_init

    def __cinit__(self, SIZE_t n_outputs, object n_classes):
        cdef SIZE_t k = 0

        self.n_outputs = n_outputs
        self.n_node_samples = 0
        self.weighted_n_node_samples = 0.0
        self.n_left = 0
        self.n_right = 0
        self.weighted_n_left = 0.0
        self.weighted_n_right = 0.0

        self.n_classes = <SIZE_t*> malloc(n_outputs * sizeof(SIZE_t))
        if self.n_classes == NULL:
            raise MemoryError()

        cdef SIZE_t label_count_stride = 0

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
                (self.n_outputs, sizet_ptr_to_ndarray(self.n_classes, self.n_outputs)),
                self.__getstate__())

    def __getstate__(self):
        return {}

    def __setstate__(self, d):
        pass


# =============================================================================
# Utils
# =============================================================================

cdef inline np.ndarray int_ptr_to_ndarray(int* data, SIZE_t size):
    """Encapsulate data into a 1D numpy array of int's."""
    cdef np.npy_intp shape[1]
    shape[0] = <np.npy_intp> size
    return np.PyArray_SimpleNewFromData(1, shape, np.NPY_INT, data)

cdef inline np.ndarray sizet_ptr_to_ndarray(SIZE_t* data, SIZE_t size):
    """Encapsulate data into a 1D numpy array of intp's."""
    cdef np.npy_intp shape[1]
    shape[0] = <np.npy_intp> size
    return np.PyArray_SimpleNewFromData(1, shape, np.NPY_INTP, data)

cdef inline np.ndarray double_ptr_to_ndarray(double* data, SIZE_t size):
    """Encapsulate data into a 1D numpy array of double's."""
    cdef np.npy_intp shape[1]
    shape[0] = <np.npy_intp> size
    return np.PyArray_SimpleNewFromData(1, shape, np.NPY_DOUBLE, data)
