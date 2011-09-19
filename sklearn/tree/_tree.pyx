# encoding: utf-8
# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
#
# Author: Peter Prettenhofer and Brian Holt
#
# License: BSD Style.


import numpy as np
cimport numpy as np

cimport cython

# Define a datatype for the data array 
DTYPE = np.float32
ctypedef np.float32_t DTYPE_t
ctypedef np.int8_t BOOL_t

cdef extern from "math.h":
    cdef extern double log(double x)
    cdef extern double pow(double base, double exponent)

cdef extern from "float.h":
    cdef extern double DBL_MAX

# Classification entropy measures
#
#    From Hastie et al. Elements of Statistical Learning, 2009.
#
#    If a target is a classification outcome taking on values 0,1,...,K-1
#    In node m, representing a region Rm with Nm observations, let
#
#       pmk = 1/ Nm \sum_{x_i in Rm} I(yi = k)
#
#    be the proportion of class k observations in node m

cdef class Criterion:
    """Interface for splitting criteria (regression and classification)"""

    cdef void init(self, DTYPE_t *y, BOOL_t *sample_mask, int n_samples,
                   int n_total_samples):
        """Initialise the criterion class for new split point."""
        pass

    cdef void reset(self):
        """Reset the criterion for a new feature index."""
        pass

    cdef int update(self, int a, int b, DTYPE_t *y, int *X_argsorted_i,
                    BOOL_t *sample_mask):
        """Update the criteria for each value in interval [a,b)

        a and b are indices in `X_argsorted_i`.
        """
        pass

    cdef double eval(self):
        """Evaluate the criteria (aka the split error)."""
        pass


cdef class ClassificationCriterion(Criterion):
    """Abstract criterion for classification.

    Attributes
    ----------
    n_classes : int
        number of classes
    n_samples : int
        number of samples
    label_count_left : int*
        label counts for samples left of splitting point.
    label_count_right : int*
        label counts for samples right of splitting point.
    label_count_init : int*
        Initial label counts for samples right of splitting point.
        Used to reset `label_count_right` for each feature.
    n_left : int
        number of samples left of splitting point.
    n_right : int
        number of samples right of splitting point.
    """
    cdef int n_classes
    cdef int n_samples
    cdef int* label_count_left
    cdef int* label_count_right
    cdef int* label_count_init
    cdef int n_left
    cdef int n_right

    # need to store ref to arrays to prevent GC
    cdef ndarray_label_count_left
    cdef ndarray_label_count_right
    cdef ndarray_label_count_init

    def __init__(self, int n_classes):
        cdef np.ndarray[np.int32_t, ndim=1] ndarray_label_count_left \
            = np.zeros((n_classes,), dtype=np.int32, order='C')
        cdef np.ndarray[np.int32_t, ndim=1] ndarray_label_count_right \
            = np.zeros((n_classes,), dtype=np.int32, order='C')
        cdef np.ndarray[np.int32_t, ndim=1] ndarray_label_count_init \
            = np.zeros((n_classes,), dtype=np.int32, order='C')

        self.n_classes = n_classes
        self.n_samples = 0
        self.n_left = 0
        self.n_right = 0
        self.label_count_left = <int*>ndarray_label_count_left.data
        self.label_count_right = <int*>ndarray_label_count_right.data
        self.label_count_init = <int*>ndarray_label_count_init.data
        self.ndarray_label_count_left = ndarray_label_count_left
        self.ndarray_label_count_right = ndarray_label_count_right
        self.ndarray_label_count_init = ndarray_label_count_init

    cdef void init(self, DTYPE_t *y, BOOL_t *sample_mask, int n_samples,
                   int n_total_samples):
        """Initialise the criterion class."""
        cdef int c = 0
        cdef int j = 0

        self.n_samples = n_samples

        for c from 0 <= c < self.n_classes:
            self.label_count_init[c] = 0

        for j from 0 <= j < n_total_samples:
            if sample_mask[j] == 0:
                continue
            c = <int>(y[j])
            self.label_count_init[c] += 1

        self.reset()

        # print "ClassificationCriterion.init: "
        # print "    label_count_left [",
        # for c from 0 <= c < self.n_classes:
        #     print self.label_count_left[c],
        # print "] ", self.n_left
        # print "    label_count_right [",
        # for c from 0 <= c < self.n_classes:
        #     print self.label_count_right[c],
        # print "] ", self.n_right

    cdef void reset(self):
        """Reset label_counts by setting `label_count_left to zero
        and copying the init array into the right."""
        cdef int c = 0
        self.n_left = 0
        self.n_right = self.n_samples

        for c from 0 <= c < self.n_classes:
            self.label_count_left[c] = 0
            self.label_count_right[c] = self.label_count_init[c]

    cdef int update(self, int a, int b, DTYPE_t *y, int *X_argsorted_i,
                    BOOL_t *sample_mask):
        """Update the criteria for each value in interval [a,b)

        a and b are indices in `X_argsorted_i`.
        """
        cdef int c
        # post condition: all samples from [0:b) are on the left side
        for idx from a <= idx < b:
            s = X_argsorted_i[idx]
            if sample_mask[s] == 0:
                continue
            c = <int>(y[s])
            self.label_count_right[c] -= 1
            self.label_count_left[c] += 1
            self.n_right -= 1
            self.n_left += 1

        # print "ClassificationCriterion.update: a = ", a, " b = ", b
        # print "    sorted [",
        # for c from 0 <= c < self.n_samples:
        #     print X_argsorted_i[c],
        # print "] "
        # print "    y [",
        # for c from 0 <= c < self.n_samples:
        #     print <int>y[X_argsorted_i[c]],
        # print "] "
        # print "    label_count_left [",
        # for c from 0 <= c < self.n_classes:
        #     print self.label_count_left[c],
        # print "] ", self.n_left
        # print "    label_count_right [",
        # for c from 0 <= c < self.n_classes:
        #     print self.label_count_right[c],
        # print "] ", self.n_right

        return self.n_left

    cdef double eval(self):
        pass


cdef class Gini(ClassificationCriterion):
    """Gini Index splitting criteria.

    Gini index = \sum_{k=0}^{K-1} pmk (1 - pmk)
               = 1 - \sum_{k=0}^{K-1} pmk ** 2
    """

    cdef double eval(self):
        """Returns Gini index of left branch + Gini index of right branch. """
        cdef double H_left = 1.0
        cdef double H_right = 1.0
        cdef int k
        cdef double e1, e2
        cdef double n_left = <double> self.n_left
        cdef double n_right = <double> self.n_right

        for k from 0 <= k < self.n_classes:
            if self.label_count_left[k] > 0:
                H_left -= ((self.label_count_left[k] / n_left)
                           * (self.label_count_left[k] / n_left))
            if self.label_count_right[k] > 0:
                H_right -= ((self.label_count_right[k] / n_right)
                            * (self.label_count_right[k] / n_right))

        e1 = (n_left / self.n_samples) * H_left
        e2 = (n_right / self.n_samples) * H_right
        return e1 + e2


cdef class Entropy(ClassificationCriterion):
    """Entropy splitting criteria.

    Cross Entropy = - \sum_{k=0}^{K-1} pmk log(pmk)
    """

    cdef double eval(self):
        """Returns Entropy of left branch + Entropy index of right branch. """
        cdef double H_left = 0.0
        cdef double H_right = 0.0
        cdef int k
        cdef double e1, e2
        cdef double n_left = <double> self.n_left
        cdef double n_right = <double> self.n_right

        for k from 0 <= k < self.n_classes:
            if self.label_count_left[k] > 0:
                H_left -= ((self.label_count_left[k] / n_left)
                           * log(self.label_count_left[k] / n_left))
            if self.label_count_right[k] > 0:
                H_right -= ((self.label_count_right[k] / n_right)
                            * log(self.label_count_right[k] / n_right))

        e1 = (n_left / self.n_samples) * H_left
        e2 = (n_right / self.n_samples) * H_right
        return e1 + e2


cdef class RegressionCriterion(Criterion):
    """Abstract criterion for regression. Computes variance of the
    target values left and right of the split point.

    Computation is linear in `n_samples` by using ::

        var = \sum_i^n (y_i - y_bar) ** 2
            = (\sum_i^n y_i ** 2) - n_samples y_bar ** 2

    Attributes
    ----------
    n_samples : int
        The number of samples
    mean_left : double
        The mean target value of the samples left of the split point.
    mean_right : double
        The mean target value of the samples right of the split.
    sq_sum_left : double
        The sum of squared target values left of the split point.
    sq_sum_right : double
        The sum of squared target values right of the split point.
    var_left : double
        The variance of the target values left of the split point.
    var_right : double
        The variance of the target values left of the split point.
    n_left : int
        number of samples left of split point.
    n_right : int
        number of samples right of split point.
    """

    cdef int n_samples
    cdef int n_right
    cdef int n_left

    cdef double mean_left
    cdef double mean_right
    cdef double mean_init

    cdef double sq_sum_right
    cdef double sq_sum_left
    cdef double sq_sum_init

    cdef double var_left
    cdef double var_right

    def __init__(self):
        self.n_samples = 0
        self.n_left = 0
        self.n_right = 0
        self.mean_left = 0.0
        self.mean_right = 0.0
        self.mean_init = 0.0
        self.sq_sum_right = 0.0
        self.sq_sum_left = 0.0
        self.sq_sum_init = 0.0
        self.var_left = 0.0
        self.var_right = 0.0

    cdef void init(self, DTYPE_t *y, BOOL_t *sample_mask, int n_samples,
                   int n_total_samples):
        """Initialise the criterion class; assume all samples
        are in the right branch and store the mean and squared
        sum in `self.mean_init` and `self.sq_sum_init`. """
        self.mean_left = 0.0
        self.mean_right = 0.0
        self.mean_init = 0.0
        self.sq_sum_right = 0.0
        self.sq_sum_left = 0.0
        self.sq_sum_init = 0.0
        self.var_left = 0.0
        self.var_right = 0.0
        self.n_samples = n_samples

        cdef int j = 0
        for j from 0 <= j < n_total_samples:
            if sample_mask[j] == 0:
                continue
            self.sq_sum_init = self.sq_sum_init + (y[j] * y[j])
            self.mean_init = self.mean_init + y[j]

        self.mean_init = self.mean_init / self.n_samples

        self.reset()

    cdef void reset(self):
        """Reset criterion for new feature.

        Assume all data in right branch and copy statistics of the
        whole dataset into the auxiliary variables of the
        right branch.
        """
        self.n_right = self.n_samples
        self.n_left = 0
        self.mean_right = self.mean_init
        self.mean_left = 0.0
        self.sq_sum_right = self.sq_sum_init
        self.sq_sum_left = 0.0
        self.var_left = 0.0
        self.var_right = self.sq_sum_right - \
            self.n_samples * (self.mean_right * self.mean_right)

    cdef int update(self, int a, int b, DTYPE_t *y, int *X_argsorted_i,
                    BOOL_t *sample_mask):
        """Update the criteria for each value in interval [a,b)

        a and b are indices in `X_argsorted_i`.
        """
        cdef double y_idx = 0.0
        cdef int idx, j
        # post condition: all samples from [0:b) are on the left side
        for idx from a <= idx < b:
            j = X_argsorted_i[idx]
            if sample_mask[j] == 0:
                continue
            y_idx = y[j]
            self.sq_sum_left = self.sq_sum_left + (y_idx * y_idx)
            self.sq_sum_right = self.sq_sum_right - (y_idx * y_idx)

            self.mean_left = (self.n_left * self.mean_left + y_idx) / \
                <double>(self.n_left + 1)
            self.mean_right = ((self.n_samples - self.n_left) * \
                self.mean_right - y_idx) / \
                <double>(self.n_samples - self.n_left - 1)

            self.n_right -= 1
            self.n_left += 1

            self.var_left = self.sq_sum_left - \
                self.n_left * (self.mean_left * self.mean_left)
            self.var_right = self.sq_sum_right - \
                self.n_right * (self.mean_right * self.mean_right)

        return self.n_left

    cdef double eval(self):
        pass


cdef class MSE(RegressionCriterion):
    """Mean squared error impurity criterion.

    MSE = var_left + var_right
    """

    cdef double eval(self):
        assert (self.n_left + self.n_right) == self.n_samples
        return (self.n_left / <double>self.n_samples) * self.var_left + \
               (self.n_right / <double>self.n_samples) * self.var_right


cdef int smallest_sample_larger_than(int sample_idx, DTYPE_t *X_i,
                                     int *X_argsorted_i, BOOL_t *sample_mask,
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
        if X_i[j] > threshold:
            return idx
    return -1


def _find_best_split(np.ndarray[DTYPE_t, ndim=2, mode="fortran"] X,
                     np.ndarray[DTYPE_t, ndim=1, mode="c"] y,
                     np.ndarray[np.int32_t, ndim=2, mode="fortran"] X_argsorted,
                     np.ndarray sample_mask,
                     np.ndarray[np.int32_t, ndim=1, mode="c"] feature_mask,
                     Criterion criterion,
                     int n_samples):
    """Find the best dimension and threshold that minimises the error.

    Parameters
    ----------
    X : ndarray, shape (n_total_samples, n_features), dtype=DTYPE_t
        The feature values.
    y : ndarray, shape (n_total_samples,), dtype=float
        The label to predict for each sample.
    X_argsorted : ndarray, shape (n_samples, n_features)
        Argsort of cols of `X`. `X_argsorted[0,j]` gives the example
        index of the smallest value of feature `j`.
    sample_mask : ndarray, shape (n_samples,), dtype=np.bool
        A mask for the samples to be considered. Only samples `j` for which
        sample_mask[j] != 0 are considered.
    feature_mask : ndarray, shape (n_samples,), dtype=int32
        A feature mask indicating active features.
    criterion : Criterion
        The criterion function to be minimized.
    n_samples : int
        The number of samples in the current sample_mask (i.e. `sample_mask.sum()`).

    Returns
    -------
    best_i : int
        The split feature or -1 if criterion not smaller than
        `parent_split_error`.
    best_t : DTYPE_t
        The split threshold
    initial_error : DTYPE_t
        The initial error contained in the node.
    """

    cdef int n_total_samples = X.shape[0]
    cdef int n_features = X.shape[1]
    cdef int i, a, b, best_i = -1
    cdef DTYPE_t t, initial_error, error
    cdef DTYPE_t best_error = np.inf, best_t = np.inf

    # Pointer access to ndarray data
    cdef DTYPE_t *y_ptr = <DTYPE_t *>y.data
    cdef DTYPE_t *X_i = NULL
    
    cdef int *X_argsorted_i = NULL

    # sample mask data pointer
    cdef BOOL_t *sample_mask_ptr = <BOOL_t *>sample_mask.data

    # Compute the column strides (increment in pointer elements to get
    # from column i to i + 1) for `X` and `X_argsorted`
    cdef int X_elem_stride = X.strides[0]
    cdef int X_col_stride = X.strides[1]
    cdef int X_stride = X_col_stride / X_elem_stride
    cdef int X_argsorted_elem_stride = X_argsorted.strides[0]
    cdef int X_argsorted_col_stride = X_argsorted.strides[1]
    cdef int X_argsorted_stride = X_argsorted_col_stride / X_argsorted_elem_stride

    # Compute the initial entropy in the node
    X_argsorted_i = <int *>X_argsorted.data
    criterion.init(y_ptr, sample_mask_ptr, n_samples, n_total_samples)
    initial_error = criterion.eval()
    if initial_error == 0: # break early if the node is pure
        return best_i, best_t, initial_error
    best_error = initial_error
    # print 'at init, best error = ', best_error

    for i from 0 <= i < n_features:
        if feature_mask[i] == 0:
            continue

        # get i-th col of X and X_sorted
        X_i = (<DTYPE_t *>X.data) + X_stride * i
        X_argsorted_i = (<int *>X_argsorted.data) + X_argsorted_stride * i

        # reset the criterion for this feature
        criterion.reset()

        # index of smallest sample in X_argsorted_i
        a = 0

        while True:
            b = smallest_sample_larger_than(a, X_i, X_argsorted_i, sample_mask_ptr,
                                            n_total_samples)

            # if -1 there's none and we are finished
            if b == -1:
                break

            criterion.update(a, b, y_ptr, X_argsorted_i, sample_mask_ptr)

            # get criterion value
            error = criterion.eval()
            # print 'error = ', error

            # check if current error is smaller than previous best
            # if this is never true best_i is -1.
            if error < best_error:
                t = (X_i[X_argsorted_i[a]] + X_i[X_argsorted_i[b]]) / 2.0
                best_i = i
                best_t = t
                best_error = error

            a = b

    return best_i, best_t, initial_error
