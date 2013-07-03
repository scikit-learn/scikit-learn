# encoding: utf-8
# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False

# Author: Peter Prettenhofer, Brian Holt, Gilles Louppe
# Licence: BSD 3 clause

from libc.stdlib cimport calloc, free, malloc, realloc
from libc.string cimport memcpy
from libc.math cimport log

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
                         SIZE_t pos,
                         SIZE_t end):
        """Initialize the criterion at node samples[start:end] and
           children samples[start:pos] and samples[pos:end]."""
        pass

    cdef void reset(self):
        """Reset the criterion at pos=0."""
        pass

    cdef void update(self, SIZE_t new_pos):
        """Update the collected statistics by moving samples[pos:new_pos] from the
           right child to the left child."""
        pass

    cdef double node_impurity(self):
        """Evaluate the impurity of the current node, i.e. the impurity of
           samples[start:end]."""
        pass

    cdef double children_impurity(self):
        """Evaluate the impurity in children nodes, i.e. the impurity of
           samples[start:pos] + the impurity of samples[pos:end]."""
        pass


cdef class ClassificationCriterion(Criterion):
    """Abstract criterion for classification."""
    cdef SIZE_t* n_classes
    cdef SIZE_t label_count_stride
    cdef double* label_count_left
    cdef double* label_count_right
    cdef double* label_count_total

    def __cinit__(self, SIZE_t n_outputs, object n_classes):
        # Default values
        self.y = NULL
        self.y_stride = 0
        self.sample_weight = NULL

        self.samples = NULL
        self.start = 0
        self.pos = 0
        self.end = 0

        self.n_outputs = n_outputs
        self.n_node_samples = 0
        self.weighted_n_node_samples = 0.0
        self.weighted_n_left = 0.0
        self.weighted_n_right = 0.0

        # Count labels for each output
        self.n_classes = <SIZE_t*> malloc(n_outputs * sizeof(SIZE_t))
        if self.n_classes == NULL:
            raise MemoryError()

        cdef SIZE_t k = 0
        cdef SIZE_t label_count_stride = 0

        for k from 0 <= k < n_outputs:
            self.n_classes[k] = n_classes[k]

            if n_classes[k] > label_count_stride:
                label_count_stride = n_classes[k]

        self.label_count_stride = label_count_stride

        # Allocate counters
        self.label_count_left = <double*> calloc(n_outputs * label_count_stride, sizeof(double))
        self.label_count_right = <double*> calloc(n_outputs * label_count_stride, sizeof(double))
        self.label_count_total = <double*> calloc(n_outputs * label_count_stride, sizeof(double))

        # Check for allocation errors
        if self.label_count_left == NULL or \
           self.label_count_right == NULL or \
           self.label_count_total == NULL:
            free(self.n_classes)
            free(self.label_count_left)
            free(self.label_count_right)
            free(self.label_count_total)
            raise MemoryError()

    def __dealloc__(self):
        """Destructor."""
        free(self.n_classes)
        free(self.label_count_left)
        free(self.label_count_right)
        free(self.label_count_total)

    def __reduce__(self):
        return (ClassificationCriterion,
                (self.n_outputs, sizet_ptr_to_ndarray(self.n_classes, self.n_outputs)),
                self.__getstate__())

    def __getstate__(self):
        return {}

    def __setstate__(self, d):
        pass

    cdef void init(self, DOUBLE_t* y,
                         SIZE_t y_stride,
                         DOUBLE_t* sample_weight,
                         SIZE_t* samples,
                         SIZE_t start,
                         SIZE_t pos,
                         SIZE_t end):
        """Initialize the criterion at node samples[start:end] and
           children samples[start:pos] and samples[pos:end]."""
        # Initialize fields
        self.y = y
        self.y_stride = y_stride
        self.sample_weight = sample_weight
        self.samples = samples
        self.start = start
        self.end = end
        self.n_node_samples = end - start
        self.weighted_n_node_samples = 0.

        # Initialize label_count_total and weighted_n_node_samples
        cdef SIZE_t n_outputs = self.n_outputs
        cdef SIZE_t* n_classes = self.n_classes
        cdef SIZE_t label_count_stride = self.label_count_stride
        cdef double* label_count_total = self.label_count_total

        cdef SIZE_t i = 0
        cdef SIZE_t p = 0
        cdef SIZE_t k = 0
        cdef SIZE_t c = 0
        cdef DOUBLE_t w = 1.0

        for k from 0 <= k < n_outputs:
            for c from 0 <= c < n_classes[k]:
                label_count_total[k * label_count_stride + c] = 0

        for p from start <= p < end:
            i = samples[p]

            if sample_weight != NULL:
                w = sample_weight[i]

            for k from 0 <= k < n_outputs:
                c = <SIZE_t> y[i * y_stride + k]
                label_count_total[k * label_count_stride + c] += w
                self.weighted_n_node_samples += w

        # Move cursor to i-th sample
        self.reset()
        self.update(pos)

    cdef void reset(self):
        """Reset the criterion at pos=0."""
        self.pos = 0

        self.weighted_n_left = 0.0
        self.weighted_n_right = self.weighted_n_node_samples

        cdef SIZE_t n_outputs = self.n_outputs
        cdef SIZE_t* n_classes = self.n_classes
        cdef SIZE_t label_count_stride = self.label_count_stride
        cdef double* label_count_total = self.label_count_total
        cdef double* label_count_left = self.label_count_left
        cdef double* label_count_right = self.label_count_right

        cdef SIZE_t k = 0
        cdef SIZE_t c = 0

        for k from 0 <= k < n_outputs:
            for c from 0 <= c < n_classes[k]:
                # Reset left label counts to 0
                label_count_left[k * label_count_stride + c] = 0

                # Reset right label counts to the initial counts
                label_count_right[k * label_count_stride + c] = label_count_total[k * label_count_stride + c]

    cdef void update(self, SIZE_t new_pos):
        """Update the collected statistics by moving samples[pos:new_pos] from the
           right child to the left child."""
        cdef DOUBLE_t* y = self.y
        cdef SIZE_t y_stride = self.y_stride
        cdef DOUBLE_t* sample_weight = self.sample_weight

        cdef SIZE_t* samples = self.samples
        cdef SIZE_t pos = self.pos

        cdef SIZE_t n_outputs = self.n_outputs
        cdef SIZE_t* n_classes = self.n_classes
        cdef SIZE_t label_count_stride = self.label_count_stride
        cdef double* label_count_total = self.label_count_total
        cdef double* label_count_left = self.label_count_left
        cdef double* label_count_right = self.label_count_right

        cdef double weighted_n_left = self.weighted_n_left
        cdef double weighted_n_right = self.weighted_n_right

        cdef SIZE_t i
        cdef SIZE_t p
        cdef SIZE_t k
        cdef SIZE_t c
        cdef DOUBLE_t w = 1.0

        # Note: We assume start <= pos < new_pos <= end

        for p from pos <= p < new_pos:
            i = samples[p]

            if sample_weight != NULL:
                w  = sample_weight[i]

            for k from 0 <= k < n_outputs:
                c = <SIZE_t> y[i * y_stride + k]
                label_count_left[k * label_count_stride + c] += w
                label_count_right[k * label_count_stride + c] -= w

            weighted_n_left += w
            weighted_n_right -= w

        self.weighted_n_left = weighted_n_left
        self.weighted_n_right = weighted_n_right

        self.pos = new_pos

        # TODO
        # # Skip splits that result in nodes with net 0 or negative weight
        # if (weighted_n_left <= 0 or
        #     (self.weighted_n_samples - weighted_n_left) <= 0):
        #     return False

        # # Prevent any single class from having a net negative weight
        # for k from 0 <= k < n_outputs:
        #     for c from 0 <= c < n_classes[k]:
        #         if (label_count_left[k * label_count_stride + c] < 0 or
        #             label_count_right[k * label_count_stride + c] < 0):
        #             return False

    cdef double node_impurity(self):
        pass

    cdef double children_impurity(self):
        pass


cdef class Entropy(ClassificationCriterion):
    """Cross Entropy impurity criteria.

    Let the target be a classification outcome taking values in 0, 1, ..., K-1.
    If node m represents a region Rm with Nm observations, then let

        pmk = 1/ Nm \sum_{x_i in Rm} I(yi = k)

    be the proportion of class k observations in node m.

    The cross-entropy is then defined as

        cross-entropy = - \sum_{k=0}^{K-1} pmk log(pmk)
    """
    cdef double node_impurity(self):
        """Evaluate the impurity of the current node, i.e. the impurity of
           samples[start:end]."""
        cdef double weighted_n_node_samples = self.weighted_n_node_samples

        cdef SIZE_t n_outputs = self.n_outputs
        cdef SIZE_t* n_classes = self.n_classes
        cdef SIZE_t label_count_stride = self.label_count_stride
        cdef double* label_count_total = self.label_count_total

        cdef double entropy = 0.0
        cdef double total = 0.0
        cdef double tmp
        cdef SIZE_t k
        cdef SIZE_t c

        for k from 0 <= k < n_outputs:
            entropy = 0.0

            for c from 0 <= c < n_classes[k]:
                tmp = label_count_total[k * label_count_stride + c]
                if tmp > 0.0:
                    tmp /= weighted_n_node_samples
                    entropy -= tmp * log(tmp)

            total += entropy

        return total / n_outputs

    cdef double children_impurity(self):
        """Evaluate the impurity in children nodes, i.e. the impurity of
           samples[start:pos] + the impurity of samples[pos:end]."""
        cdef double weighted_n_node_samples = self.weighted_n_node_samples
        cdef double weighted_n_left = self.weighted_n_left
        cdef double weighted_n_right = self.weighted_n_right

        cdef SIZE_t n_outputs = self.n_outputs
        cdef SIZE_t* n_classes = self.n_classes
        cdef SIZE_t label_count_stride = self.label_count_stride
        cdef double* label_count_left = self.label_count_left
        cdef double* label_count_right = self.label_count_right

        cdef double entropy_left = 0.0
        cdef double entropy_right = 0.0
        cdef double total = 0.0
        cdef double tmp
        cdef SIZE_t k
        cdef SIZE_t c

        for k from 0 <= k < n_outputs:
            entropy_left = 0.0
            entropy_right = 0.0

            for c from 0 <= c < n_classes[k]:
                tmp = label_count_left[k * label_count_stride + c]
                if tmp > 0.0:
                    tmp /= weighted_n_left
                    entropy_left -= tmp * log(tmp)

                tmp = label_count_right[k * label_count_stride + c]
                if tmp > 0.0:
                    tmp /= weighted_n_right
                    entropy_right -= tmp * log(tmp)

            total += weighted_n_left * entropy_left
            total += weighted_n_right * entropy_right

        return total / (weighted_n_node_samples * n_outputs)


cdef class Gini(ClassificationCriterion):
    """Gini Index impurity criteria.

    Let the target be a classification outcome taking values in 0, 1, ..., K-1.
    If node m represents a region Rm with Nm observations, then let

        pmk = 1/ Nm \sum_{x_i in Rm} I(yi = k)

    be the proportion of class k observations in node m.

    The Gini Index is then defined as:

        index = \sum_{k=0}^{K-1} pmk (1 - pmk)
              = 1 - \sum_{k=0}^{K-1} pmk ** 2
    """
    cdef double node_impurity(self):
        """Evaluate the impurity of the current node, i.e. the impurity of
           samples[start:end]."""
        cdef double weighted_n_node_samples = self.weighted_n_node_samples

        cdef SIZE_t n_outputs = self.n_outputs
        cdef SIZE_t* n_classes = self.n_classes
        cdef SIZE_t label_count_stride = self.label_count_stride
        cdef double* label_count_total = self.label_count_total

        cdef double gini = 0.0
        cdef double total = 0.0
        cdef double tmp
        cdef SIZE_t k
        cdef SIZE_t c

        for k from 0 <= k < n_outputs:
            gini = 0.0

            for c from 0 <= c < n_classes[k]:
                tmp = label_count_total[k * label_count_stride + c]
                gini += tmp * tmp

            if weighted_n_node_samples <= 0.0:
                gini = 0.0
            else:
                gini = 1.0 - gini / (weighted_n_node_samples * weighted_n_node_samples)

            total += gini

        return total / n_outputs

    cdef double children_impurity(self):
        """Evaluate the impurity in children nodes, i.e. the impurity of
           samples[start:pos] + the impurity of samples[pos:end]."""
        cdef double weighted_n_node_samples = self.weighted_n_node_samples
        cdef double weighted_n_left = self.weighted_n_left
        cdef double weighted_n_right = self.weighted_n_right

        cdef SIZE_t n_outputs = self.n_outputs
        cdef SIZE_t* n_classes = self.n_classes
        cdef SIZE_t label_count_stride = self.label_count_stride
        cdef double* label_count_left = self.label_count_left
        cdef double* label_count_right = self.label_count_right

        cdef double gini_left = 0.0
        cdef double gini_right = 0.0
        cdef double total = 0.0
        cdef double tmp
        cdef SIZE_t k
        cdef SIZE_t c

        for k from 0 <= k < n_outputs:
            gini_left = 0.0
            gini_right = 0.0

            for c from 0 <= c < n_classes[k]:
                tmp = label_count_left[k * label_count_stride + c]
                gini_left += tmp * tmp
                tmp = label_count_right[k * label_count_stride + c]
                gini_right += tmp * tmp

            if weighted_n_left <= 0.0:
                gini_left = 0.0
            else:
                gini_left = 1.0 - gini_left / (weighted_n_left * weighted_n_left)

            if weighted_n_right <= 0.0:
                gini_right = 0.0
            else:
                gini_right = 1.0 - gini_right / (weighted_n_right * weighted_n_right)

            total += weighted_n_left * gini_left
            total += weighted_n_right * gini_right

        return total / (weighted_n_node_samples * n_outputs)


cdef class RegressionCriterion(Criterion):
    """Abstract criterion for regression.

    Computes variance of the target values left and right of the split point.
    Computation is linear in `n_samples` by using ::

        var = \sum_i^n (y_i - y_bar) ** 2
            = (\sum_i^n y_i ** 2) - n_samples y_bar ** 2
    """
    cdef double* mean_left
    cdef double* mean_right
    cdef double* mean_total
    cdef double* sq_sum_left
    cdef double* sq_sum_right
    cdef double* sq_sum_total
    cdef double* var_left
    cdef double* var_right

    def __cinit__(self, SIZE_t n_outputs):
        # Default values
        self.y = NULL
        self.y_stride = 0
        self.sample_weight = NULL

        self.samples = NULL
        self.start = 0
        self.pos = 0
        self.end = 0

        self.n_outputs = n_outputs
        self.n_node_samples = 0
        self.weighted_n_node_samples = 0.0
        self.weighted_n_left = 0.0
        self.weighted_n_right = 0.0

        # Allocate accumulators
        self.mean_left = <double*> calloc(n_outputs, sizeof(double))
        self.mean_right = <double*> calloc(n_outputs, sizeof(double))
        self.mean_total = <double*> calloc(n_outputs, sizeof(double))
        self.sq_sum_left = <double*> calloc(n_outputs, sizeof(double))
        self.sq_sum_right = <double*> calloc(n_outputs, sizeof(double))
        self.sq_sum_total = <double*> calloc(n_outputs, sizeof(double))
        self.var_left = <double*> calloc(n_outputs, sizeof(double))
        self.var_right = <double*> calloc(n_outputs, sizeof(double))

        # Check for allocation errors
        if self.mean_left == NULL or \
           self.mean_right == NULL or \
           self.mean_total == NULL or \
           self.sq_sum_left == NULL or \
           self.sq_sum_right == NULL or \
           self.sq_sum_total == NULL or \
           self.var_left == NULL or \
           self.var_right == NULL:
            free(self.mean_left)
            free(self.mean_right)
            free(self.mean_total)
            free(self.sq_sum_left)
            free(self.sq_sum_right)
            free(self.sq_sum_total)
            free(self.var_left)
            free(self.var_right)
            raise MemoryError()

    def __dealloc__(self):
        """Destructor."""
        free(self.mean_left)
        free(self.mean_right)
        free(self.mean_total)
        free(self.sq_sum_left)
        free(self.sq_sum_right)
        free(self.sq_sum_total)
        free(self.var_left)
        free(self.var_right)

    def __reduce__(self):
        return (RegressionCriterion, (self.n_outputs,), self.__getstate__())

    def __getstate__(self):
        return {}

    def __setstate__(self, d):
        pass

    cdef void init(self, DOUBLE_t* y,
                         SIZE_t y_stride,
                         DOUBLE_t* sample_weight,
                         SIZE_t* samples,
                         SIZE_t start,
                         SIZE_t pos,
                         SIZE_t end):
        """Initialize the criterion at node samples[start:end] and
           children samples[start:pos] and samples[pos:end]."""
        # Initialize fields
        self.y = y
        self.y_stride = y_stride
        self.sample_weight = sample_weight
        self.samples = samples
        self.start = start
        self.end = end
        self.n_node_samples = end - start
        self.weighted_n_node_samples = 0.

        # Initialize accumulators
        cdef SIZE_t n_outputs = self.n_outputs
        cdef double* mean_left = self.mean_left
        cdef double* mean_right = self.mean_right
        cdef double* mean_total = self.mean_total
        cdef double* sq_sum_left = self.sq_sum_left
        cdef double* sq_sum_right = self.sq_sum_right
        cdef double* sq_sum_total = self.sq_sum_total
        cdef double* var_left = self.var_left
        cdef double* var_right = self.var_right

        cdef SIZE_t i = 0
        cdef SIZE_t p = 0
        cdef SIZE_t k = 0
        cdef DOUBLE_t y_ik = 0.0
        cdef DOUBLE_t w = 1.0

        for k from 0 <= k < n_outputs:
            mean_left[k] = 0.0
            mean_right[k] = 0.0
            mean_total[k] = 0.0
            sq_sum_right[k] = 0.0
            sq_sum_left[k] = 0.0
            sq_sum_total[k] = 0.0
            var_left[k] = 0.0
            var_right[k] = 0.0

        for p from start <= p < end:
            i = samples[p]

            if sample_weight != NULL:
                w = sample_weight[i]

            for k from 0 <= k < n_outputs:
                y_ik = y[i * y_stride + k]
                sq_sum_total[k] += w * y_ik * y_ik
                mean_total[k] += w * y_ik
                self.weighted_n_node_samples += w

        for k from 0 <= k < n_outputs:
            mean_total[k] /= self.weighted_n_node_samples

        # Move cursor to i-th sample
        self.reset()
        self.update(pos)

    cdef void reset(self):
        """Reset the criterion at pos=0."""
        self.pos = 0

        self.weighted_n_left = 0.0
        self.weighted_n_right = self.weighted_n_node_samples

        cdef SIZE_t n_outputs = self.n_outputs
        cdef double* mean_left = self.mean_left
        cdef double* mean_right = self.mean_right
        cdef double* mean_total = self.mean_total
        cdef double* sq_sum_left = self.sq_sum_left
        cdef double* sq_sum_right = self.sq_sum_right
        cdef double* sq_sum_total = self.sq_sum_total
        cdef double* var_left = self.var_left
        cdef double* var_right = self.var_right

        cdef SIZE_t k = 0

        for k from 0 <= k < n_outputs:
            mean_right[k] = mean_total[k]
            mean_left[k] = 0.0
            sq_sum_right[k] = sq_sum_total[k]
            sq_sum_left[k] = 0.0
            var_left[k] = 0.0
            var_right[k] = (sq_sum_right[k] - self.weighted_n_node_samples * (mean_right[k] * mean_right[k]))

    cdef void update(self, SIZE_t new_pos):
        """Update the collected statistics by moving samples[pos:new_pos] from the
           right child to the left child."""
        cdef DOUBLE_t* y = self.y
        cdef SIZE_t y_stride = self.y_stride
        cdef DOUBLE_t* sample_weight = self.sample_weight

        cdef SIZE_t* samples = self.samples
        cdef SIZE_t pos = self.pos

        cdef SIZE_t n_outputs = self.n_outputs
        cdef double* mean_left = self.mean_left
        cdef double* mean_right = self.mean_right
        cdef double* sq_sum_left = self.sq_sum_left
        cdef double* sq_sum_right = self.sq_sum_right
        cdef double* var_left = self.var_left
        cdef double* var_right = self.var_right

        cdef double weighted_n_left = self.weighted_n_left
        cdef double weighted_n_right = self.weighted_n_right

        cdef SIZE_t i
        cdef SIZE_t p
        cdef SIZE_t k
        cdef DOUBLE_t w = 1.0
        cdef DOUBLE_t y_ik = 0.0

        # Note: We assume start <= pos < new_pos <= end

        for p from pos <= p < new_pos:
            i = samples[p]

            if sample_weight != NULL:
                w  = sample_weight[i]

            for k from 0 <= k < n_outputs:
                y_ik = y[i * y_stride + k]
                sq_sum_left[k] += w * (y_ik * y_ik)
                sq_sum_right[k] -= w * (y_ik * y_ik)

                mean_left[k] = ((weighted_n_left * mean_left[k] + w * y_ik) / (weighted_n_left + w))
                mean_right[k] = ((weighted_n_right * mean_right[k] - w * y_ik) / (weighted_n_right - w))

            weighted_n_left += w
            weighted_n_right -= w

        for k from 0 <= k < n_outputs:
            var_left[k] = sq_sum_left[k] - weighted_n_left * (mean_left[k] * mean_left[k])
            var_right[k] = sq_sum_right[k] - weighted_n_right * (mean_right[k] * mean_right[k])

        self.weighted_n_left = weighted_n_left
        self.weighted_n_right = weighted_n_right

        self.pos = new_pos

        # TODO
        # # Skip splits that result in nodes with net 0 or negative weight
        # if (weighted_n_left <= 0 or
        #     (self.weighted_n_samples - weighted_n_left) <= 0):
        #     return False

cdef class MSE(RegressionCriterion):
    """Mean squared error impurity criterion.

        MSE = var_left + var_right
    """
    cdef double node_impurity(self):
        """Evaluate the impurity of the current node, i.e. the impurity of
           samples[start:end]."""
        cdef SIZE_t n_outputs = self.n_outputs
        cdef double* sq_sum_total = self.sq_sum_total
        cdef double* mean_total = self.mean_total
        cdef double weighted_n_node_samples = self.weighted_n_node_samples
        cdef double total = 0.0
        cdef SIZE_t k

        for k from 0 <= k < n_outputs:
            total += sq_sum_total[k] - weighted_n_node_samples * (mean_total[k] * mean_total[k])

        return total / n_outputs

    cdef double children_impurity(self):
        """Evaluate the impurity in children nodes, i.e. the impurity of
           samples[start:pos] + the impurity of samples[pos:end]."""
        cdef SIZE_t n_outputs = self.n_outputs
        cdef double* var_left = self.var_left
        cdef double* var_right = self.var_right
        cdef double total = 0.0
        cdef SIZE_t k

        for k from 0 <= k < n_outputs:
            total += var_left[k]
            total += var_right[k]

        return total / n_outputs


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
