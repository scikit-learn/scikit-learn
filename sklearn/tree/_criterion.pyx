# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False

# Authors: Gilles Louppe <g.louppe@gmail.com>
#          Peter Prettenhofer <peter.prettenhofer@gmail.com>
#          Brian Holt <bdholt1@gmail.com>
#          Noel Dawe <noel@dawe.me>
#          Satrajit Gosh <satrajit.ghosh@gmail.com>
#          Lars Buitinck
#          Arnaud Joly <arnaud.v.joly@gmail.com>
#          Joel Nothman <joel.nothman@gmail.com>
#          Fares Hedayati <fares.hedayati@gmail.com>
#          Jacob Schreiber <jmschreiber91@gmail.com>
#          Nelson Liu <nelson@nelsonliu.me>
#
# License: BSD 3 clause

from libc.stdlib cimport calloc
from libc.stdlib cimport free
from libc.string cimport memcpy
from libc.string cimport memset
from libc.math cimport fabs

import numpy as np
cimport numpy as np
np.import_array()

from ._utils cimport log
from ._utils cimport safe_realloc
from ._utils cimport sizet_ptr_to_ndarray
from ._utils cimport WeightedMedianCalculator

cdef class Criterion:
    """Interface for impurity criteria.

    This object stores methods on how to calculate how good a split is using
    different metrics.
    """

    def __dealloc__(self):
        """Destructor."""

        free(self.sum_total)
        free(self.sum_left)
        free(self.sum_right)

    def __getstate__(self):
        return {}

    def __setstate__(self, d):
        pass

    cdef void init(self, DOUBLE_t* y, SIZE_t y_stride, DOUBLE_t* sample_weight,
                   double weighted_n_samples, SIZE_t* samples, SIZE_t start,
                   SIZE_t end) nogil:
        """Placeholder for a method which will initialize the criterion.

        Parameters
        ----------
        y: array-like, dtype=DOUBLE_t
            y is a buffer that can store values for n_outputs target variables
        y_stride: SIZE_t
            y_stride is used to index the kth output value as follows:
            y[i, k] = y[i * y_stride + k]
        sample_weight: array-like, dtype=DOUBLE_t
            The weight of each sample
        weighted_n_samples: DOUBLE_t
            The total weight of the samples being considered
        samples: array-like, dtype=DOUBLE_t
            Indices of the samples in X and y, where samples[start:end]
            correspond to the samples in this node
        start: SIZE_t
            The first sample to be used on this node
        end: SIZE_t
            The last sample used on this node

        """

        pass

    cdef void reset(self) nogil:
        """Reset the criterion at pos=start.

        This method must be implemented by the subclass.
        """

        pass

    cdef void reverse_reset(self) nogil:
        """Reset the criterion at pos=end.

        This method must be implemented by the subclass.
        """
        pass

    cdef void update(self, SIZE_t new_pos) nogil:
        """Updated statistics by moving samples[pos:new_pos] to the left child.

        This updates the collected statistics by moving samples[pos:new_pos]
        from the right child to the left child. It must be implemented by
        the subclass.

        Parameters
        ----------
        new_pos: SIZE_t
            New starting index position of the samples in the right child
        """

        pass

    cdef double node_impurity(self) nogil:
        """Placeholder for calculating the impurity of the node.

        Placeholder for a method which will evaluate the impurity of
        the current node, i.e. the impurity of samples[start:end]. This is the
        primary function of the criterion class.
        """

        pass

    cdef void children_impurity(self, double* impurity_left,
                                double* impurity_right) nogil:
        """Placeholder for calculating the impurity of children.

        Placeholder for a method which evaluates the impurity in
        children nodes, i.e. the impurity of samples[start:pos] + the impurity
        of samples[pos:end].

        Parameters
        ----------
        impurity_left: double pointer
            The memory address where the impurity of the left child should be
            stored.
        impurity_right: double pointer
            The memory address where the impurity of the right child should be
            stored
        """

        pass

    cdef void node_value(self, double* dest) nogil:
        """Placeholder for storing the node value.

        Placeholder for a method which will compute the node value
        of samples[start:end] and save the value into dest.

        Parameters
        ----------
        dest: double pointer
            The memory address where the node value should be stored.
        """

        pass

    cdef double proxy_impurity_improvement(self) nogil:
        """Compute a proxy of the impurity reduction

        This method is used to speed up the search for the best split.
        It is a proxy quantity such that the split that maximizes this value
        also maximizes the impurity improvement. It neglects all constant terms
        of the impurity decrease for a given split.

        The absolute impurity improvement is only computed by the
        impurity_improvement method once the best split has been found.
        """
        cdef double impurity_left
        cdef double impurity_right
        self.children_impurity(&impurity_left, &impurity_right)

        return (- self.weighted_n_right * impurity_right
                - self.weighted_n_left * impurity_left)

    cdef double impurity_improvement(self, double impurity) nogil:
        """Placeholder for improvement in impurity after a split.

        Placeholder for a method which computes the improvement
        in impurity when a split occurs. The weighted impurity improvement
        equation is the following:

            N_t / N * (impurity - N_t_R / N_t * right_impurity
                                - N_t_L / N_t * left_impurity)

        where N is the total number of samples, N_t is the number of samples
        at the current node, N_t_L is the number of samples in the left child,
        and N_t_R is the number of samples in the right child,

        Parameters
        ----------
        impurity: double
            The initial impurity of the node before the split

        Return
        ------
        double: improvement in impurity after the split occurs
        """

        cdef double impurity_left
        cdef double impurity_right

        self.children_impurity(&impurity_left, &impurity_right)

        return ((self.weighted_n_node_samples / self.weighted_n_samples) *
                (impurity - (self.weighted_n_right / 
                             self.weighted_n_node_samples * impurity_right)
                          - (self.weighted_n_left / 
                             self.weighted_n_node_samples * impurity_left)))


cdef class ClassificationCriterion(Criterion):
    """Abstract criterion for classification."""

    cdef SIZE_t* n_classes
    cdef SIZE_t sum_stride

    def __cinit__(self, SIZE_t n_outputs,
                  np.ndarray[SIZE_t, ndim=1] n_classes):
        """Initialize attributes for this criterion.

        Parameters
        ----------
        n_outputs: SIZE_t
            The number of targets, the dimensionality of the prediction
        n_classes: numpy.ndarray, dtype=SIZE_t
            The number of unique classes in each target
        """

        self.y = NULL
        self.y_stride = 0
        self.sample_weight = NULL

        self.samples = NULL
        self.start = 0
        self.pos = 0
        self.end = 0

        self.n_outputs = n_outputs
        self.n_samples = 0
        self.n_node_samples = 0
        self.weighted_n_node_samples = 0.0
        self.weighted_n_left = 0.0
        self.weighted_n_right = 0.0

        # Count labels for each output
        self.sum_total = NULL
        self.sum_left = NULL
        self.sum_right = NULL
        self.n_classes = NULL

        safe_realloc(&self.n_classes, n_outputs)

        cdef SIZE_t k = 0
        cdef SIZE_t sum_stride = 0

        # For each target, set the number of unique classes in that target,
        # and also compute the maximal stride of all targets
        for k in range(n_outputs):
            self.n_classes[k] = n_classes[k]

            if n_classes[k] > sum_stride:
                sum_stride = n_classes[k]

        self.sum_stride = sum_stride

        cdef SIZE_t n_elements = n_outputs * sum_stride
        self.sum_total = <double*> calloc(n_elements, sizeof(double))
        self.sum_left = <double*> calloc(n_elements, sizeof(double))
        self.sum_right = <double*> calloc(n_elements, sizeof(double))

        if (self.sum_total == NULL or 
                self.sum_left == NULL or
                self.sum_right == NULL):
            raise MemoryError()

    def __dealloc__(self):
        """Destructor."""
        free(self.n_classes)

    def __reduce__(self):
        return (type(self),
                (self.n_outputs,
                 sizet_ptr_to_ndarray(self.n_classes, self.n_outputs)),
                self.__getstate__())

    cdef void init(self, DOUBLE_t* y, SIZE_t y_stride,
                   DOUBLE_t* sample_weight, double weighted_n_samples,
                   SIZE_t* samples, SIZE_t start, SIZE_t end) nogil:
        """Initialize the criterion at node samples[start:end] and
        children samples[start:start] and samples[start:end].

        Parameters
        ----------
        y: array-like, dtype=DOUBLE_t
            The target stored as a buffer for memory efficiency
        y_stride: SIZE_t
            The stride between elements in the buffer, important if there
            are multiple targets (multi-output)
        sample_weight: array-like, dtype=DTYPE_t
            The weight of each sample
        weighted_n_samples: SIZE_t
            The total weight of all samples
        samples: array-like, dtype=SIZE_t
            A mask on the samples, showing which ones we want to use
        start: SIZE_t
            The first sample to use in the mask
        end: SIZE_t
            The last sample to use in the mask
        """

        self.y = y
        self.y_stride = y_stride
        self.sample_weight = sample_weight
        self.samples = samples
        self.start = start
        self.end = end
        self.n_node_samples = end - start
        self.weighted_n_samples = weighted_n_samples
        self.weighted_n_node_samples = 0.0

        cdef SIZE_t* n_classes = self.n_classes
        cdef double* sum_total = self.sum_total

        cdef SIZE_t i
        cdef SIZE_t p
        cdef SIZE_t k
        cdef SIZE_t c
        cdef DOUBLE_t w = 1.0
        cdef SIZE_t offset = 0

        for k in range(self.n_outputs):
            memset(sum_total + offset, 0, n_classes[k] * sizeof(double))
            offset += self.sum_stride

        for p in range(start, end):
            i = samples[p]

            # w is originally set to be 1.0, meaning that if no sample weights
            # are given, the default weight of each sample is 1.0
            if sample_weight != NULL:
                w = sample_weight[i]

            # Count weighted class frequency for each target
            for k in range(self.n_outputs):
                c = <SIZE_t> y[i * y_stride + k]
                sum_total[k * self.sum_stride + c] += w

            self.weighted_n_node_samples += w

        # Reset to pos=start
        self.reset()

    cdef void reset(self) nogil:
        """Reset the criterion at pos=start."""

        self.pos = self.start

        self.weighted_n_left = 0.0
        self.weighted_n_right = self.weighted_n_node_samples

        cdef double* sum_total = self.sum_total
        cdef double* sum_left = self.sum_left
        cdef double* sum_right = self.sum_right

        cdef SIZE_t* n_classes = self.n_classes
        cdef SIZE_t k

        for k in range(self.n_outputs):
            memset(sum_left, 0, n_classes[k] * sizeof(double))
            memcpy(sum_right, sum_total, n_classes[k] * sizeof(double))

            sum_total += self.sum_stride
            sum_left += self.sum_stride
            sum_right += self.sum_stride

    cdef void reverse_reset(self) nogil:
        """Reset the criterion at pos=end."""
        self.pos = self.end

        self.weighted_n_left = self.weighted_n_node_samples
        self.weighted_n_right = 0.0

        cdef double* sum_total = self.sum_total
        cdef double* sum_left = self.sum_left
        cdef double* sum_right = self.sum_right

        cdef SIZE_t* n_classes = self.n_classes
        cdef SIZE_t k

        for k in range(self.n_outputs):
            memset(sum_right, 0, n_classes[k] * sizeof(double))
            memcpy(sum_left, sum_total, n_classes[k] * sizeof(double))

            sum_total += self.sum_stride
            sum_left += self.sum_stride
            sum_right += self.sum_stride

    cdef void update(self, SIZE_t new_pos) nogil:
        """Updated statistics by moving samples[pos:new_pos] to the left child.

        Parameters
        ----------
        new_pos: SIZE_t
            The new ending position for which to move samples from the right
            child to the left child.
        """
        cdef DOUBLE_t* y = self.y
        cdef SIZE_t pos = self.pos
        cdef SIZE_t end = self.end

        cdef double* sum_left = self.sum_left
        cdef double* sum_right = self.sum_right
        cdef double* sum_total = self.sum_total

        cdef SIZE_t* n_classes = self.n_classes
        cdef SIZE_t* samples = self.samples
        cdef DOUBLE_t* sample_weight = self.sample_weight

        cdef SIZE_t i
        cdef SIZE_t p
        cdef SIZE_t k
        cdef SIZE_t c
        cdef SIZE_t label_index
        cdef DOUBLE_t w = 1.0

        # Update statistics up to new_pos
        #
        # Given that
        #   sum_left[x] +  sum_right[x] = sum_total[x]
        # and that sum_total is known, we are going to update
        # sum_left from the direction that require the least amount
        # of computations, i.e. from pos to new_pos or from end to new_po.

        if (new_pos - pos) <= (end - new_pos):
            for p in range(pos, new_pos):
                i = samples[p]

                if sample_weight != NULL:
                    w = sample_weight[i]

                for k in range(self.n_outputs):
                    label_index = (k * self.sum_stride +
                                   <SIZE_t> y[i * self.y_stride + k])
                    sum_left[label_index] += w

                self.weighted_n_left += w

        else:
            self.reverse_reset()

            for p in range(end - 1, new_pos - 1, -1):
                i = samples[p]

                if sample_weight != NULL:
                    w = sample_weight[i]

                for k in range(self.n_outputs):
                    label_index = (k * self.sum_stride +
                                   <SIZE_t> y[i * self.y_stride + k])
                    sum_left[label_index] -= w

                self.weighted_n_left -= w

        # Update right part statistics
        self.weighted_n_right = self.weighted_n_node_samples - self.weighted_n_left
        for k in range(self.n_outputs):
            for c in range(n_classes[k]):
                sum_right[c] = sum_total[c] - sum_left[c]

            sum_right += self.sum_stride
            sum_left += self.sum_stride
            sum_total += self.sum_stride

        self.pos = new_pos

    cdef double node_impurity(self) nogil:
        pass

    cdef void children_impurity(self, double* impurity_left,
                                double* impurity_right) nogil:
        pass

    cdef void node_value(self, double* dest) nogil:
        """Compute the node value of samples[start:end] and save it into dest.

        Parameters
        ----------
        dest: double pointer
            The memory address which we will save the node value into.
        """

        cdef double* sum_total = self.sum_total
        cdef SIZE_t* n_classes = self.n_classes
        cdef SIZE_t k

        for k in range(self.n_outputs):
            memcpy(dest, sum_total, n_classes[k] * sizeof(double))
            dest += self.sum_stride
            sum_total += self.sum_stride


cdef class Entropy(ClassificationCriterion):
    """Cross Entropy impurity criterion.

    This handles cases where the target is a classification taking values
    0, 1, ... K-2, K-1. If node m represents a region Rm with Nm observations,
    then let

        count_k = 1 / Nm \sum_{x_i in Rm} I(yi = k)

    be the proportion of class k observations in node m.

    The cross-entropy is then defined as

        cross-entropy = -\sum_{k=0}^{K-1} count_k log(count_k)
    """

    cdef double node_impurity(self) nogil:
        """Evaluate the impurity of the current node, i.e. the impurity of
        samples[start:end], using the cross-entropy criterion."""

        cdef SIZE_t* n_classes = self.n_classes
        cdef double* sum_total = self.sum_total
        cdef double entropy = 0.0
        cdef double count_k
        cdef SIZE_t k
        cdef SIZE_t c

        for k in range(self.n_outputs):
            for c in range(n_classes[k]):
                count_k = sum_total[c]
                if count_k > 0.0:
                    count_k /= self.weighted_n_node_samples
                    entropy -= count_k * log(count_k)

            sum_total += self.sum_stride

        return entropy / self.n_outputs

    cdef void children_impurity(self, double* impurity_left,
                                double* impurity_right) nogil:
        """Evaluate the impurity in children nodes

        i.e. the impurity of the left child (samples[start:pos]) and the
        impurity the right child (samples[pos:end]).

        Parameters
        ----------
        impurity_left: double pointer
            The memory address to save the impurity of the left node
        impurity_right: double pointer
            The memory address to save the impurity of the right node
        """

        cdef SIZE_t* n_classes = self.n_classes
        cdef double* sum_left = self.sum_left
        cdef double* sum_right = self.sum_right
        cdef double entropy_left = 0.0
        cdef double entropy_right = 0.0
        cdef double count_k
        cdef SIZE_t k
        cdef SIZE_t c

        for k in range(self.n_outputs):
            for c in range(n_classes[k]):
                count_k = sum_left[c]
                if count_k > 0.0:
                    count_k /= self.weighted_n_left
                    entropy_left -= count_k * log(count_k)

                count_k = sum_right[c]
                if count_k > 0.0:
                    count_k /= self.weighted_n_right
                    entropy_right -= count_k * log(count_k)

            sum_left += self.sum_stride
            sum_right += self.sum_stride

        impurity_left[0] = entropy_left / self.n_outputs
        impurity_right[0] = entropy_right / self.n_outputs


cdef class Gini(ClassificationCriterion):
    """Gini Index impurity criterion.

    This handles cases where the target is a classification taking values
    0, 1, ... K-2, K-1. If node m represents a region Rm with Nm observations,
    then let

        count_k = 1/ Nm \sum_{x_i in Rm} I(yi = k)

    be the proportion of class k observations in node m.

    The Gini Index is then defined as:

        index = \sum_{k=0}^{K-1} count_k (1 - count_k)
              = 1 - \sum_{k=0}^{K-1} count_k ** 2
    """

    cdef double node_impurity(self) nogil:
        """Evaluate the impurity of the current node, i.e. the impurity of
        samples[start:end] using the Gini criterion."""


        cdef SIZE_t* n_classes = self.n_classes
        cdef double* sum_total = self.sum_total
        cdef double gini = 0.0
        cdef double sq_count
        cdef double count_k
        cdef SIZE_t k
        cdef SIZE_t c

        for k in range(self.n_outputs):
            sq_count = 0.0

            for c in range(n_classes[k]):
                count_k = sum_total[c]
                sq_count += count_k * count_k

            gini += 1.0 - sq_count / (self.weighted_n_node_samples *
                                      self.weighted_n_node_samples)

            sum_total += self.sum_stride

        return gini / self.n_outputs

    cdef void children_impurity(self, double* impurity_left,
                                double* impurity_right) nogil:
        """Evaluate the impurity in children nodes

        i.e. the impurity of the left child (samples[start:pos]) and the
        impurity the right child (samples[pos:end]) using the Gini index.

        Parameters
        ----------
        impurity_left: DTYPE_t
            The memory address to save the impurity of the left node to
        impurity_right: DTYPE_t
            The memory address to save the impurity of the right node to
        """

        cdef SIZE_t* n_classes = self.n_classes
        cdef double* sum_left = self.sum_left
        cdef double* sum_right = self.sum_right
        cdef double gini_left = 0.0
        cdef double gini_right = 0.0
        cdef double sq_count_left
        cdef double sq_count_right
        cdef double count_k
        cdef SIZE_t k
        cdef SIZE_t c

        for k in range(self.n_outputs):
            sq_count_left = 0.0
            sq_count_right = 0.0

            for c in range(n_classes[k]):
                count_k = sum_left[c]
                sq_count_left += count_k * count_k

                count_k = sum_right[c]
                sq_count_right += count_k * count_k

            gini_left += 1.0 - sq_count_left / (self.weighted_n_left *
                                                self.weighted_n_left)

            gini_right += 1.0 - sq_count_right / (self.weighted_n_right *
                                                  self.weighted_n_right)

            sum_left += self.sum_stride
            sum_right += self.sum_stride

        impurity_left[0] = gini_left / self.n_outputs
        impurity_right[0] = gini_right / self.n_outputs


cdef class RegressionCriterion(Criterion):
    """Abstract regression criterion.

    This handles cases where the target is a continuous value, and is
    evaluated by computing the variance of the target values left and right
    of the split point. The computation takes linear time with `n_samples`
    by using ::

        var = \sum_i^n (y_i - y_bar) ** 2
            = (\sum_i^n y_i ** 2) - n_samples * y_bar ** 2
    """

    cdef double sq_sum_total

    def __cinit__(self, SIZE_t n_outputs, SIZE_t n_samples):
        """Initialize parameters for this criterion.

        Parameters
        ----------
        n_outputs: SIZE_t
            The number of targets to be predicted

        n_samples: SIZE_t
            The total number of samples to fit on
        """

        # Default values
        self.y = NULL
        self.y_stride = 0
        self.sample_weight = NULL

        self.samples = NULL
        self.start = 0
        self.pos = 0
        self.end = 0

        self.n_outputs = n_outputs
        self.n_samples = n_samples
        self.n_node_samples = 0
        self.weighted_n_node_samples = 0.0
        self.weighted_n_left = 0.0
        self.weighted_n_right = 0.0

        self.sq_sum_total = 0.0

        # Allocate accumulators. Make sure they are NULL, not uninitialized,
        # before an exception can be raised (which triggers __dealloc__).
        self.sum_total = NULL
        self.sum_left = NULL
        self.sum_right = NULL

        # Allocate memory for the accumulators
        self.sum_total = <double*> calloc(n_outputs, sizeof(double))
        self.sum_left = <double*> calloc(n_outputs, sizeof(double))
        self.sum_right = <double*> calloc(n_outputs, sizeof(double))

        if (self.sum_total == NULL or 
                self.sum_left == NULL or
                self.sum_right == NULL):
            raise MemoryError()

    def __reduce__(self):
        return (type(self), (self.n_outputs, self.n_samples), self.__getstate__())

    cdef void init(self, DOUBLE_t* y, SIZE_t y_stride, DOUBLE_t* sample_weight,
                   double weighted_n_samples, SIZE_t* samples, SIZE_t start,
                   SIZE_t end) nogil:
        """Initialize the criterion at node samples[start:end] and
           children samples[start:start] and samples[start:end]."""
        # Initialize fields
        self.y = y
        self.y_stride = y_stride
        self.sample_weight = sample_weight
        self.samples = samples
        self.start = start
        self.end = end
        self.n_node_samples = end - start
        self.weighted_n_samples = weighted_n_samples
        self.weighted_n_node_samples = 0.

        cdef SIZE_t i
        cdef SIZE_t p
        cdef SIZE_t k
        cdef DOUBLE_t y_ik
        cdef DOUBLE_t w_y_ik
        cdef DOUBLE_t w = 1.0

        self.sq_sum_total = 0.0
        memset(self.sum_total, 0, self.n_outputs * sizeof(double))

        for p in range(start, end):
            i = samples[p]

            if sample_weight != NULL:
                w = sample_weight[i]

            for k in range(self.n_outputs):
                y_ik = y[i * y_stride + k]
                w_y_ik = w * y_ik
                self.sum_total[k] += w_y_ik
                self.sq_sum_total += w_y_ik * y_ik

            self.weighted_n_node_samples += w

        # Reset to pos=start
        self.reset()

    cdef void reset(self) nogil:
        """Reset the criterion at pos=start."""
        cdef SIZE_t n_bytes = self.n_outputs * sizeof(double)
        memset(self.sum_left, 0, n_bytes)
        memcpy(self.sum_right, self.sum_total, n_bytes)

        self.weighted_n_left = 0.0
        self.weighted_n_right = self.weighted_n_node_samples
        self.pos = self.start

    cdef void reverse_reset(self) nogil:
        """Reset the criterion at pos=end."""
        cdef SIZE_t n_bytes = self.n_outputs * sizeof(double)
        memset(self.sum_right, 0, n_bytes)
        memcpy(self.sum_left, self.sum_total, n_bytes)

        self.weighted_n_right = 0.0
        self.weighted_n_left = self.weighted_n_node_samples
        self.pos = self.end

    cdef void update(self, SIZE_t new_pos) nogil:
        """Updated statistics by moving samples[pos:new_pos] to the left."""

        cdef double* sum_left = self.sum_left
        cdef double* sum_right = self.sum_right
        cdef double* sum_total = self.sum_total

        cdef double* sample_weight = self.sample_weight
        cdef SIZE_t* samples = self.samples

        cdef DOUBLE_t* y = self.y
        cdef SIZE_t pos = self.pos
        cdef SIZE_t end = self.end
        cdef SIZE_t i
        cdef SIZE_t p
        cdef SIZE_t k
        cdef DOUBLE_t w = 1.0
        cdef DOUBLE_t y_ik

        # Update statistics up to new_pos
        #
        # Given that
        #           sum_left[x] +  sum_right[x] = sum_total[x]
        # and that sum_total is known, we are going to update
        # sum_left from the direction that require the least amount
        # of computations, i.e. from pos to new_pos or from end to new_po.

        if (new_pos - pos) <= (end - new_pos):
            for p in range(pos, new_pos):
                i = samples[p]

                if sample_weight != NULL:
                    w = sample_weight[i]

                for k in range(self.n_outputs):
                    y_ik = y[i * self.y_stride + k]
                    sum_left[k] += w * y_ik

                self.weighted_n_left += w
        else:
            self.reverse_reset()

            for p in range(end - 1, new_pos - 1, -1):
                i = samples[p]

                if sample_weight != NULL:
                    w = sample_weight[i]

                for k in range(self.n_outputs):
                    y_ik = y[i * self.y_stride + k]
                    sum_left[k] -= w * y_ik

                self.weighted_n_left -= w

        self.weighted_n_right = (self.weighted_n_node_samples - 
                                 self.weighted_n_left)
        for k in range(self.n_outputs):
            sum_right[k] = sum_total[k] - sum_left[k]

        self.pos = new_pos

    cdef double node_impurity(self) nogil:
        pass

    cdef void children_impurity(self, double* impurity_left,
                                double* impurity_right) nogil:
        pass

    cdef void node_value(self, double* dest) nogil:
        """Compute the node value of samples[start:end] into dest."""

        cdef SIZE_t k

        for k in range(self.n_outputs):
            dest[k] = self.sum_total[k] / self.weighted_n_node_samples


cdef class MSE(RegressionCriterion):
    """Mean squared error impurity criterion.

        MSE = var_left + var_right
    """

    cdef double node_impurity(self) nogil:
        """Evaluate the impurity of the current node, i.e. the impurity of
           samples[start:end]."""

        cdef double* sum_total = self.sum_total
        cdef double impurity
        cdef SIZE_t k

        impurity = self.sq_sum_total / self.weighted_n_node_samples
        for k in range(self.n_outputs):
            impurity -= (sum_total[k] / self.weighted_n_node_samples)**2.0

        return impurity / self.n_outputs

    cdef double proxy_impurity_improvement(self) nogil:
        """Compute a proxy of the impurity reduction

        This method is used to speed up the search for the best split.
        It is a proxy quantity such that the split that maximizes this value
        also maximizes the impurity improvement. It neglects all constant terms
        of the impurity decrease for a given split.

        The absolute impurity improvement is only computed by the
        impurity_improvement method once the best split has been found.
        """

        cdef double* sum_left = self.sum_left
        cdef double* sum_right = self.sum_right

        cdef SIZE_t k
        cdef double proxy_impurity_left = 0.0
        cdef double proxy_impurity_right = 0.0

        for k in range(self.n_outputs):
            proxy_impurity_left += sum_left[k] * sum_left[k]
            proxy_impurity_right += sum_right[k] * sum_right[k]

        return (proxy_impurity_left / self.weighted_n_left +
                proxy_impurity_right / self.weighted_n_right)

    cdef void children_impurity(self, double* impurity_left,
                                double* impurity_right) nogil:
        """Evaluate the impurity in children nodes, i.e. the impurity of the
           left child (samples[start:pos]) and the impurity the right child
           (samples[pos:end])."""


        cdef DOUBLE_t* y = self.y
        cdef DOUBLE_t* sample_weight = self.sample_weight
        cdef SIZE_t* samples = self.samples
        cdef SIZE_t pos = self.pos
        cdef SIZE_t start = self.start

        cdef double* sum_left = self.sum_left
        cdef double* sum_right = self.sum_right

        cdef double sq_sum_left = 0.0
        cdef double sq_sum_right

        cdef SIZE_t i
        cdef SIZE_t p
        cdef SIZE_t k
        cdef DOUBLE_t w = 1.0
        cdef DOUBLE_t y_ik

        for p in range(start, pos):
            i = samples[p]

            if sample_weight != NULL:
                w = sample_weight[i]

            for k in range(self.n_outputs):
                y_ik = y[i * self.y_stride + k]
                sq_sum_left += w * y_ik * y_ik

        sq_sum_right = self.sq_sum_total - sq_sum_left

        impurity_left[0] = sq_sum_left / self.weighted_n_left
        impurity_right[0] = sq_sum_right / self.weighted_n_right

        for k in range(self.n_outputs):
            impurity_left[0] -= (sum_left[k] / self.weighted_n_left) ** 2.0
            impurity_right[0] -= (sum_right[k] / self.weighted_n_right) ** 2.0 

        impurity_left[0] /= self.n_outputs
        impurity_right[0] /= self.n_outputs

cdef class MAE(RegressionCriterion):
    """Mean absolute error impurity criterion

       MAE = (1 / n)*(\sum_i |y_i - f_i|), where y_i is the true
       value and f_i is the predicted value."""
    def __dealloc__(self):
        """Destructor."""
        free(self.node_medians)

    cdef np.ndarray left_child
    cdef np.ndarray right_child
    cdef DOUBLE_t* node_medians

    def __cinit__(self, SIZE_t n_outputs, SIZE_t n_samples):
        """Initialize parameters for this criterion.

        Parameters
        ----------
        n_outputs: SIZE_t
            The number of targets to be predicted

        n_samples: SIZE_t
            The total number of samples to fit on
        """

        # Default values
        self.y = NULL
        self.y_stride = 0
        self.sample_weight = NULL

        self.samples = NULL
        self.start = 0
        self.pos = 0
        self.end = 0

        self.n_outputs = n_outputs
        self.n_samples = n_samples
        self.n_node_samples = 0
        self.weighted_n_node_samples = 0.0
        self.weighted_n_left = 0.0
        self.weighted_n_right = 0.0

        # Allocate accumulators. Make sure they are NULL, not uninitialized,
        # before an exception can be raised (which triggers __dealloc__).
        self.node_medians = NULL

        # Allocate memory for the accumulators
        safe_realloc(&self.node_medians, n_outputs)

        if (self.node_medians == NULL):
            raise MemoryError()

        self.left_child = np.empty(n_outputs, dtype='object')
        self.right_child = np.empty(n_outputs, dtype='object')
        # initialize WeightedMedianCalculators
        for k in range(n_outputs):
            self.left_child[k] = WeightedMedianCalculator(n_samples)
            self.right_child[k] = WeightedMedianCalculator(n_samples)

    cdef void init(self, DOUBLE_t* y, SIZE_t y_stride, DOUBLE_t* sample_weight,
                   double weighted_n_samples, SIZE_t* samples, SIZE_t start,
                   SIZE_t end) nogil:
        """Initialize the criterion at node samples[start:end] and
           children samples[start:start] and samples[start:end]."""

        cdef SIZE_t i, p, k
        cdef DOUBLE_t y_ik
        cdef DOUBLE_t w = 1.0

        # Initialize fields
        self.y = y
        self.y_stride = y_stride
        self.sample_weight = sample_weight
        self.samples = samples
        self.start = start
        self.end = end
        self.n_node_samples = end - start
        self.weighted_n_samples = weighted_n_samples
        self.weighted_n_node_samples = 0.

        cdef void** left_child
        cdef void** right_child

        left_child = <void**> self.left_child.data
        right_child = <void**> self.right_child.data

        for k in range(self.n_outputs):
            (<WeightedMedianCalculator> left_child[k]).reset()
            (<WeightedMedianCalculator> right_child[k]).reset()

        for p in range(start, end):
            i = samples[p]

            if sample_weight != NULL:
                w = sample_weight[i]

            for k in range(self.n_outputs):
                y_ik = y[i * y_stride + k]

                # push all values to the right side,
                # since pos = start initially anyway
                (<WeightedMedianCalculator> right_child[k]).push(y_ik, w)

            self.weighted_n_node_samples += w
        # calculate the node medians
        for k in range(self.n_outputs):
            self.node_medians[k] = (<WeightedMedianCalculator> right_child[k]).get_median()

        # Reset to pos=start
        self.reset()

    cdef void reset(self) nogil:
        """Reset the criterion at pos=start."""

        cdef SIZE_t i, k
        cdef DOUBLE_t value
        cdef DOUBLE_t weight

        cdef void** left_child = <void**> self.left_child.data
        cdef void** right_child = <void**> self.right_child.data

        self.weighted_n_left = 0.0
        self.weighted_n_right = self.weighted_n_node_samples
        self.pos = self.start

        # reset the WeightedMedianCalculators, left should have no
        # elements and right should have all elements.

        for k in range(self.n_outputs):
            # if left has no elements, it's already reset
            for i in range((<WeightedMedianCalculator> left_child[k]).size()):
                # remove everything from left and put it into right
                (<WeightedMedianCalculator> left_child[k]).pop(&value,
                                                               &weight)
                (<WeightedMedianCalculator> right_child[k]).push(value,
                                                                 weight)

    cdef void reverse_reset(self) nogil:
        """Reset the criterion at pos=end."""

        self.weighted_n_right = 0.0
        self.weighted_n_left = self.weighted_n_node_samples
        self.pos = self.end

        cdef DOUBLE_t value
        cdef DOUBLE_t weight
        cdef void** left_child = <void**> self.left_child.data
        cdef void** right_child = <void**> self.right_child.data

        # reverse reset the WeightedMedianCalculators, right should have no
        # elements and left should have all elements.
        for k in range(self.n_outputs):
            # if right has no elements, it's already reset
            for i in range((<WeightedMedianCalculator> right_child[k]).size()):
                # remove everything from right and put it into left
                (<WeightedMedianCalculator> right_child[k]).pop(&value,
                                                                &weight)
                (<WeightedMedianCalculator> left_child[k]).push(value,
                                                                weight)

    cdef void update(self, SIZE_t new_pos) nogil:
        """Updated statistics by moving samples[pos:new_pos] to the left."""

        cdef DOUBLE_t* sample_weight = self.sample_weight
        cdef SIZE_t* samples = self.samples

        cdef void** left_child = <void**> self.left_child.data
        cdef void** right_child = <void**> self.right_child.data

        cdef DOUBLE_t* y = self.y
        cdef SIZE_t pos = self.pos
        cdef SIZE_t end = self.end
        cdef SIZE_t i, p, k
        cdef DOUBLE_t w = 1.0
        cdef DOUBLE_t y_ik

        # Update statistics up to new_pos
        #
        # We are going to update right_child and left_child
        # from the direction that require the least amount of
        # computations, i.e. from pos to new_pos or from end to new_pos.

        if (new_pos - pos) <= (end - new_pos):
            for p in range(pos, new_pos):
                i = samples[p]

                if sample_weight != NULL:
                    w = sample_weight[i]

                for k in range(self.n_outputs):
                    y_ik = y[i * self.y_stride + k]
                    # remove y_ik and its weight w from right and add to left
                    (<WeightedMedianCalculator> right_child[k]).remove(y_ik, w)
                    (<WeightedMedianCalculator> left_child[k]).push(y_ik, w)

                self.weighted_n_left += w
        else:
            self.reverse_reset()

            for p in range(end - 1, new_pos - 1, -1):
                i = samples[p]

                if sample_weight != NULL:
                    w = sample_weight[i]

                for k in range(self.n_outputs):
                    y_ik = y[i * self.y_stride + k]
                    # remove y_ik and its weight w from left and add to right
                    (<WeightedMedianCalculator> left_child[k]).remove(y_ik, w)
                    (<WeightedMedianCalculator> right_child[k]).push(y_ik, w)

                self.weighted_n_left -= w

        self.weighted_n_right = (self.weighted_n_node_samples -
                                 self.weighted_n_left)
        self.pos = new_pos

    cdef void node_value(self, double* dest) nogil:
        """Computes the node value of samples[start:end] into dest."""

        cdef SIZE_t k
        for k in range(self.n_outputs):
            dest[k] = <double> self.node_medians[k]

    cdef double node_impurity(self) nogil:
        """Evaluate the impurity of the current node, i.e. the impurity of
           samples[start:end]"""

        cdef DOUBLE_t* y = self.y
        cdef DOUBLE_t* sample_weight = self.sample_weight
        cdef SIZE_t* samples = self.samples
        cdef SIZE_t i, p, k
        cdef DOUBLE_t y_ik
        cdef DOUBLE_t w_y_ik

        cdef double impurity = 0.0

        for k in range(self.n_outputs):
            for p in range(self.start, self.end):
                i = samples[p]

                y_ik = y[i * self.y_stride + k]

                impurity += <double> fabs((<double> y_ik) - <double> self.node_medians[k])
        return impurity / (self.weighted_n_node_samples * self.n_outputs)

    cdef void children_impurity(self, double* impurity_left,
                                double* impurity_right) nogil:
        """Evaluate the impurity in children nodes, i.e. the impurity of the
           left child (samples[start:pos]) and the impurity the right child
           (samples[pos:end]).
        """

        cdef DOUBLE_t* y = self.y
        cdef DOUBLE_t* sample_weight = self.sample_weight
        cdef SIZE_t* samples = self.samples

        cdef SIZE_t start = self.start
        cdef SIZE_t pos = self.pos
        cdef SIZE_t end = self.end

        cdef SIZE_t i, p, k
        cdef DOUBLE_t y_ik
        cdef DOUBLE_t median

        cdef void** left_child = <void**> self.left_child.data
        cdef void** right_child = <void**> self.right_child.data

        impurity_left[0] = 0.0
        impurity_right[0] = 0.0

        for k in range(self.n_outputs):
            median = (<WeightedMedianCalculator> left_child[k]).get_median()
            for p in range(start, pos):
                i = samples[p]

                y_ik = y[i * self.y_stride + k]

                impurity_left[0] += <double>fabs((<double> y_ik) -
                                                 <double> median)
        impurity_left[0] /= <double>((self.weighted_n_left) * self.n_outputs)

        for k in range(self.n_outputs):
            median = (<WeightedMedianCalculator> right_child[k]).get_median()
            for p in range(pos, end):
                i = samples[p]

                y_ik = y[i * self.y_stride + k]

                impurity_right[0] += <double>fabs((<double> y_ik) -
                                                  <double> median)
        impurity_right[0] /= <double>((self.weighted_n_right) *
                                      self.n_outputs)


cdef class FriedmanMSE(MSE):
    """Mean squared error impurity criterion with improvement score by Friedman

    Uses the formula (35) in Friedmans original Gradient Boosting paper:

        diff = mean_left - mean_right
        improvement = n_left * n_right * diff^2 / (n_left + n_right)
    """

    cdef double proxy_impurity_improvement(self) nogil:
        """Compute a proxy of the impurity reduction

        This method is used to speed up the search for the best split.
        It is a proxy quantity such that the split that maximizes this value
        also maximizes the impurity improvement. It neglects all constant terms
        of the impurity decrease for a given split.

        The absolute impurity improvement is only computed by the
        impurity_improvement method once the best split has been found.
        """

        cdef double* sum_left = self.sum_left
        cdef double* sum_right = self.sum_right

        cdef double total_sum_left = 0.0
        cdef double total_sum_right = 0.0

        cdef SIZE_t k
        cdef double diff = 0.0

        for k in range(self.n_outputs):
            total_sum_left += sum_left[k]
            total_sum_right += sum_right[k]

        diff = (self.weighted_n_right * total_sum_left -
                self.weighted_n_left * total_sum_right)

        return diff * diff / (self.weighted_n_left * self.weighted_n_right)

    cdef double impurity_improvement(self, double impurity) nogil:
        cdef double* sum_left = self.sum_left
        cdef double* sum_right = self.sum_right

        cdef double total_sum_left = 0.0
        cdef double total_sum_right = 0.0

        cdef SIZE_t k
        cdef double diff = 0.0

        for k in range(self.n_outputs):
            total_sum_left += sum_left[k]
            total_sum_right += sum_right[k]

        diff = (self.weighted_n_right * total_sum_left -
                self.weighted_n_left * total_sum_right) / self.n_outputs

        return (diff * diff / (self.weighted_n_left * self.weighted_n_right * 
                               self.weighted_n_node_samples))
