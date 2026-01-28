# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from libc.string cimport memcpy
from libc.string cimport memset
from libc.math cimport INFINITY

import numpy as np
cimport numpy as cnp
cnp.import_array()

from scipy.special.cython_special cimport xlogy

from sklearn.tree._utils cimport log
from sklearn.tree._utils cimport WeightedFenwickTree
from sklearn.tree._partitioner cimport sort

# EPSILON is used in the Poisson criterion
cdef float64_t EPSILON = 10 * np.finfo('double').eps

cdef class Criterion:
    """Interface for impurity criteria.

    This object stores methods on how to calculate how good a split is using
    different metrics.
    """
    def __getstate__(self):
        return {}

    def __setstate__(self, d):
        pass

    cdef int init(
        self,
        const float64_t[:, ::1] y,
        const float64_t[:] sample_weight,
        float64_t weighted_n_samples,
        const intp_t[:] sample_indices,
        intp_t start,
        intp_t end,
    ) except -1 nogil:
        """Placeholder for a method which will initialize the criterion.

        Returns -1 in case of failure to allocate memory (and raise MemoryError)
        or 0 otherwise.

        Parameters
        ----------
        y : ndarray, dtype=float64_t
            y is a buffer that can store values for n_outputs target variables
            stored as a Cython memoryview.
        sample_weight : ndarray, dtype=float64_t
            The weight of each sample stored as a Cython memoryview.
        weighted_n_samples : float64_t
            The total weight of the samples being considered
        sample_indices : ndarray, dtype=intp_t
            A mask on the samples. Indices of the samples in X and y we want to use,
            where sample_indices[start:end] correspond to the samples in this node.
        start : intp_t
            The first sample to be used on this node
        end : intp_t
            The last sample used on this node

        """
        pass

    cdef void init_missing(self, intp_t n_missing) noexcept nogil:
        """Initialize sum_missing if there are missing values.

        This method assumes that caller placed the missing samples in
        self.sample_indices[-n_missing:]

        Parameters
        ----------
        n_missing: intp_t
            Number of missing values for specific feature.
        """
        pass

    cdef int reset(self) except -1 nogil:
        """Reset the criterion at pos=start.

        This method must be implemented by the subclass.
        """
        pass

    cdef int reverse_reset(self) except -1 nogil:
        """Reset the criterion at pos=end.

        This method must be implemented by the subclass.
        """
        pass

    cdef int update(self, intp_t new_pos) except -1 nogil:
        """Updated statistics by moving sample_indices[pos:new_pos] to the left child.

        This updates the collected statistics by moving sample_indices[pos:new_pos]
        from the right child to the left child. It must be implemented by
        the subclass.

        Parameters
        ----------
        new_pos : intp_t
            New starting index position of the sample_indices in the right child
        """
        pass

    cdef float64_t node_impurity(self) noexcept nogil:
        """Placeholder for calculating the impurity of the node.

        Placeholder for a method which will evaluate the impurity of
        the current node, i.e. the impurity of sample_indices[start:end]. This is the
        primary function of the criterion class. The smaller the impurity the
        better.
        """
        pass

    cdef void children_impurity(self, float64_t* impurity_left,
                                float64_t* impurity_right) noexcept nogil:
        """Placeholder for calculating the impurity of children.

        Placeholder for a method which evaluates the impurity in
        children nodes, i.e. the impurity of sample_indices[start:pos] + the impurity
        of sample_indices[pos:end].

        Parameters
        ----------
        impurity_left : float64_t pointer
            The memory address where the impurity of the left child should be
            stored.
        impurity_right : float64_t pointer
            The memory address where the impurity of the right child should be
            stored
        """
        pass

    cdef void node_value(self, float64_t* dest) noexcept nogil:
        """Placeholder for storing the node value.

        Placeholder for a method which will compute the node value
        of sample_indices[start:end] and save the value into dest.

        Parameters
        ----------
        dest : float64_t pointer
            The memory address where the node value should be stored.
        """
        pass

    cdef void clip_node_value(self, float64_t* dest, float64_t lower_bound, float64_t upper_bound) noexcept nogil:
        pass

    cdef float64_t middle_value(self) noexcept nogil:
        """Compute the middle value of a split for monotonicity constraints

        This method is implemented in ClassificationCriterion and RegressionCriterion.
        """
        pass

    cdef float64_t proxy_impurity_improvement(self) noexcept nogil:
        """Compute a proxy of the impurity reduction.

        This method is used to speed up the search for the best split.
        It is a proxy quantity such that the split that maximizes this value
        also maximizes the impurity improvement. It neglects all constant terms
        of the impurity decrease for a given split.

        The absolute impurity improvement is only computed by the
        impurity_improvement method once the best split has been found.
        """
        cdef float64_t impurity_left
        cdef float64_t impurity_right
        self.children_impurity(&impurity_left, &impurity_right)

        return (- self.weighted_n_right * impurity_right
                - self.weighted_n_left * impurity_left)

    cdef float64_t impurity_improvement(self, float64_t impurity_parent,
                                        float64_t impurity_left,
                                        float64_t impurity_right) noexcept nogil:
        """Compute the improvement in impurity.

        This method computes the improvement in impurity when a split occurs.
        The weighted impurity improvement equation is the following:

            N_t / N * (impurity - N_t_R / N_t * right_impurity
                                - N_t_L / N_t * left_impurity)

        where N is the total number of samples, N_t is the number of samples
        at the current node, N_t_L is the number of samples in the left child,
        and N_t_R is the number of samples in the right child,

        Parameters
        ----------
        impurity_parent : float64_t
            The initial impurity of the parent node before the split

        impurity_left : float64_t
            The impurity of the left child

        impurity_right : float64_t
            The impurity of the right child

        Return
        ------
        float64_t : improvement in impurity after the split occurs
        """
        return ((self.weighted_n_node_samples / self.weighted_n_samples) *
                (impurity_parent - (self.weighted_n_right /
                                    self.weighted_n_node_samples * impurity_right)
                                 - (self.weighted_n_left /
                                    self.weighted_n_node_samples * impurity_left)))

    cdef bint check_monotonicity(
        self,
        cnp.int8_t monotonic_cst,
        float64_t lower_bound,
        float64_t upper_bound,
    ) noexcept nogil:
        pass

    cdef inline bint _check_monotonicity(
        self,
        cnp.int8_t monotonic_cst,
        float64_t lower_bound,
        float64_t upper_bound,
        float64_t value_left,
        float64_t value_right,
    ) noexcept nogil:
        cdef:
            bint check_lower_bound = (
                (value_left >= lower_bound) &
                (value_right >= lower_bound)
            )
            bint check_upper_bound = (
                (value_left <= upper_bound) &
                (value_right <= upper_bound)
            )
            bint check_monotonic_cst = (
                (value_left - value_right) * monotonic_cst <= 0
            )
        return check_lower_bound & check_upper_bound & check_monotonic_cst

    cdef void init_sum_missing(self):
        """Init sum_missing to hold sums for missing values."""

cdef inline void _move_sums_classification(
    ClassificationCriterion criterion,
    float64_t[:, ::1] sum_1,
    float64_t[:, ::1] sum_2,
    float64_t* weighted_n_1,
    float64_t* weighted_n_2,
    bint put_missing_in_1,
) noexcept nogil:
    """Distribute sum_total and sum_missing into sum_1 and sum_2.

    If there are missing values and:
    - put_missing_in_1 is True, then missing values to go sum_1. Specifically:
        sum_1 = sum_missing
        sum_2 = sum_total - sum_missing

    - put_missing_in_1 is False, then missing values go to sum_2. Specifically:
        sum_1 = 0
        sum_2 = sum_total
    """
    cdef intp_t k, c, n_bytes
    if criterion.n_missing != 0 and put_missing_in_1:
        for k in range(criterion.n_outputs):
            n_bytes = criterion.n_classes[k] * sizeof(float64_t)
            memcpy(&sum_1[k, 0], &criterion.sum_missing[k, 0], n_bytes)

        for k in range(criterion.n_outputs):
            for c in range(criterion.n_classes[k]):
                sum_2[k, c] = criterion.sum_total[k, c] - criterion.sum_missing[k, c]

        weighted_n_1[0] = criterion.weighted_n_missing
        weighted_n_2[0] = criterion.weighted_n_node_samples - criterion.weighted_n_missing
    else:
        # Assigning sum_2 = sum_total for all outputs.
        for k in range(criterion.n_outputs):
            n_bytes = criterion.n_classes[k] * sizeof(float64_t)
            memset(&sum_1[k, 0], 0, n_bytes)
            memcpy(&sum_2[k, 0], &criterion.sum_total[k, 0], n_bytes)

        weighted_n_1[0] = 0.0
        weighted_n_2[0] = criterion.weighted_n_node_samples


cdef class ClassificationCriterion(Criterion):
    """Abstract criterion for classification."""

    def __cinit__(self, intp_t n_outputs,
                  cnp.ndarray[intp_t, ndim=1] n_classes):
        """Initialize attributes for this criterion.

        Parameters
        ----------
        n_outputs : intp_t
            The number of targets, the dimensionality of the prediction
        n_classes : numpy.ndarray, dtype=intp_t
            The number of unique classes in each target
        """
        self.start = 0
        self.pos = 0
        self.end = 0
        self.missing_go_to_left = 0

        self.n_outputs = n_outputs
        self.n_samples = 0
        self.n_node_samples = 0
        self.weighted_n_node_samples = 0.0
        self.weighted_n_left = 0.0
        self.weighted_n_right = 0.0
        self.weighted_n_missing = 0.0

        self.n_classes = np.empty(n_outputs, dtype=np.intp)

        cdef intp_t k = 0
        cdef intp_t max_n_classes = 0

        # For each target, set the number of unique classes in that target,
        # and also compute the maximal stride of all targets
        for k in range(n_outputs):
            self.n_classes[k] = n_classes[k]

            if n_classes[k] > max_n_classes:
                max_n_classes = n_classes[k]

        self.max_n_classes = max_n_classes

        # Count labels for each output
        self.sum_total = np.zeros((n_outputs, max_n_classes), dtype=np.float64)
        self.sum_left = np.zeros((n_outputs, max_n_classes), dtype=np.float64)
        self.sum_right = np.zeros((n_outputs, max_n_classes), dtype=np.float64)

    def __reduce__(self):
        return (type(self),
                (self.n_outputs, np.asarray(self.n_classes)), self.__getstate__())

    cdef int init(
        self,
        const float64_t[:, ::1] y,
        const float64_t[:] sample_weight,
        float64_t weighted_n_samples,
        const intp_t[:] sample_indices,
        intp_t start,
        intp_t end
    ) except -1 nogil:
        """Initialize the criterion.

        This initializes the criterion at node sample_indices[start:end] and children
        sample_indices[start:start] and sample_indices[start:end].

        Returns -1 in case of failure to allocate memory (and raise MemoryError)
        or 0 otherwise.

        Parameters
        ----------
        y : ndarray, dtype=float64_t
            The target stored as a buffer for memory efficiency.
        sample_weight : ndarray, dtype=float64_t
            The weight of each sample stored as a Cython memoryview.
        weighted_n_samples : float64_t
            The total weight of all samples
        sample_indices : ndarray, dtype=intp_t
            A mask on the samples. Indices of the samples in X and y we want to use,
            where sample_indices[start:end] correspond to the samples in this node.
        start : intp_t
            The first sample to use in the mask
        end : intp_t
            The last sample to use in the mask
        """
        self.y = y
        self.sample_weight = sample_weight
        self.sample_indices = sample_indices
        self.start = start
        self.end = end
        self.n_node_samples = end - start
        self.weighted_n_samples = weighted_n_samples
        self.weighted_n_node_samples = 0.0

        cdef intp_t i
        cdef intp_t p
        cdef intp_t k
        cdef intp_t c
        cdef float64_t w = 1.0

        for k in range(self.n_outputs):
            memset(&self.sum_total[k, 0], 0, self.n_classes[k] * sizeof(float64_t))

        for p in range(start, end):
            i = sample_indices[p]

            # w is originally set to be 1.0, meaning that if no sample weights
            # are given, the default weight of each sample is 1.0.
            if sample_weight is not None:
                w = sample_weight[i]

            # Count weighted class frequency for each target
            for k in range(self.n_outputs):
                c = <intp_t> self.y[i, k]
                self.sum_total[k, c] += w

            self.weighted_n_node_samples += w

        # Reset to pos=start
        self.reset()
        return 0

    cdef void init_sum_missing(self):
        """Init sum_missing to hold sums for missing values."""
        self.sum_missing = np.zeros((self.n_outputs, self.max_n_classes), dtype=np.float64)

    cdef void init_missing(self, intp_t n_missing) noexcept nogil:
        """Initialize sum_missing if there are missing values.

        This method assumes that caller placed the missing samples in
        self.sample_indices[-n_missing:]
        """
        cdef intp_t i, p, k, c
        cdef float64_t w = 1.0

        self.n_missing = n_missing
        if n_missing == 0:
            return

        memset(&self.sum_missing[0, 0], 0, self.max_n_classes * self.n_outputs * sizeof(float64_t))

        self.weighted_n_missing = 0.0

        # The missing samples are assumed to be in self.sample_indices[-n_missing:]
        for p in range(self.end - n_missing, self.end):
            i = self.sample_indices[p]
            if self.sample_weight is not None:
                w = self.sample_weight[i]

            for k in range(self.n_outputs):
                c = <intp_t> self.y[i, k]
                self.sum_missing[k, c] += w

            self.weighted_n_missing += w

    cdef int reset(self) except -1 nogil:
        """Reset the criterion at pos=start.

        Returns -1 in case of failure to allocate memory (and raise MemoryError)
        or 0 otherwise.
        """
        self.pos = self.start
        _move_sums_classification(
            self,
            self.sum_left,
            self.sum_right,
            &self.weighted_n_left,
            &self.weighted_n_right,
            self.missing_go_to_left,
        )
        return 0

    cdef int reverse_reset(self) except -1 nogil:
        """Reset the criterion at pos=end.

        Returns -1 in case of failure to allocate memory (and raise MemoryError)
        or 0 otherwise.
        """
        self.pos = self.end
        _move_sums_classification(
            self,
            self.sum_right,
            self.sum_left,
            &self.weighted_n_right,
            &self.weighted_n_left,
            not self.missing_go_to_left
        )
        return 0

    cdef int update(self, intp_t new_pos) except -1 nogil:
        """Updated statistics by moving sample_indices[pos:new_pos] to the left child.

        Returns -1 in case of failure to allocate memory (and raise MemoryError)
        or 0 otherwise.

        Parameters
        ----------
        new_pos : intp_t
            The new ending position for which to move sample_indices from the right
            child to the left child.
        """
        cdef intp_t pos = self.pos
        # The missing samples are assumed to be in
        # self.sample_indices[-self.n_missing:] that is
        # self.sample_indices[end_non_missing:self.end].
        cdef intp_t end_non_missing = self.end - self.n_missing
        cdef intp_t i
        cdef intp_t p
        cdef intp_t k
        cdef intp_t c
        cdef float64_t w = 1.0

        # Update statistics up to new_pos
        #
        # Given that
        #   sum_left[x] +  sum_right[x] = sum_total[x]
        # and that sum_total is known, we are going to update
        # sum_left from the direction that require the least amount
        # of computations, i.e. from pos to new_pos or from end to new_po.
        if (new_pos - pos) <= (end_non_missing - new_pos):
            for p in range(pos, new_pos):
                i = self.sample_indices[p]

                if self.sample_weight is not None:
                    w = self.sample_weight[i]

                for k in range(self.n_outputs):
                    self.sum_left[k, <intp_t> self.y[i, k]] += w

                self.weighted_n_left += w

        else:
            self.reverse_reset()

            for p in range(end_non_missing - 1, new_pos - 1, -1):
                i = self.sample_indices[p]

                if self.sample_weight is not None:
                    w = self.sample_weight[i]

                for k in range(self.n_outputs):
                    self.sum_left[k, <intp_t> self.y[i, k]] -= w

                self.weighted_n_left -= w

        # Update right part statistics
        self.weighted_n_right = self.weighted_n_node_samples - self.weighted_n_left
        for k in range(self.n_outputs):
            for c in range(self.n_classes[k]):
                self.sum_right[k, c] = self.sum_total[k, c] - self.sum_left[k, c]

        self.pos = new_pos
        return 0

    cdef float64_t node_impurity(self) noexcept nogil:
        pass

    cdef void children_impurity(self, float64_t* impurity_left,
                                float64_t* impurity_right) noexcept nogil:
        pass

    cdef void node_value(self, float64_t* dest) noexcept nogil:
        """Compute the node value of sample_indices[start:end] and save it into dest.

        Parameters
        ----------
        dest : float64_t pointer
            The memory address which we will save the node value into.
        """
        cdef intp_t k, c

        for k in range(self.n_outputs):
            for c in range(self.n_classes[k]):
                dest[c] = self.sum_total[k, c] / self.weighted_n_node_samples
            dest += self.max_n_classes

    cdef inline void clip_node_value(
        self, float64_t * dest, float64_t lower_bound, float64_t upper_bound
    ) noexcept nogil:
        """Clip the values in dest such that predicted probabilities stay between
        `lower_bound` and `upper_bound` when monotonic constraints are enforced.
        Note that monotonicity constraints are only supported for:
        - single-output trees and
        - binary classifications.
        """
        if dest[0] < lower_bound:
            dest[0] = lower_bound
        elif dest[0] > upper_bound:
            dest[0] = upper_bound

        # Values for binary classification must sum to 1.
        dest[1] = 1 - dest[0]

    cdef inline float64_t middle_value(self) noexcept nogil:
        """Compute the middle value of a split for monotonicity constraints as the simple average
        of the left and right children values.

        Note that monotonicity constraints are only supported for:
        - single-output trees and
        - binary classifications.
        """
        return (
            (self.sum_left[0, 0] / (2 * self.weighted_n_left)) +
            (self.sum_right[0, 0] / (2 * self.weighted_n_right))
        )

    cdef inline bint check_monotonicity(
        self,
        cnp.int8_t monotonic_cst,
        float64_t lower_bound,
        float64_t upper_bound,
    ) noexcept nogil:
        """Check monotonicity constraint is satisfied at the current classification split"""
        cdef:
            float64_t value_left = self.sum_left[0][0] / self.weighted_n_left
            float64_t value_right = self.sum_right[0][0] / self.weighted_n_right

        return self._check_monotonicity(monotonic_cst, lower_bound, upper_bound, value_left, value_right)


cdef class Entropy(ClassificationCriterion):
    r"""Cross Entropy impurity criterion.

    This handles cases where the target is a classification taking values
    0, 1, ... K-2, K-1. If node m represents a region Rm with Nm observations,
    then let

        count_k = 1 / Nm \sum_{x_i in Rm} I(yi = k)

    be the proportion of class k observations in node m.

    The cross-entropy is then defined as

        cross-entropy = -\sum_{k=0}^{K-1} count_k log(count_k)
    """

    cdef float64_t node_impurity(self) noexcept nogil:
        """Evaluate the impurity of the current node.

        Evaluate the cross-entropy criterion as impurity of the current node,
        i.e. the impurity of sample_indices[start:end]. The smaller the impurity the
        better.
        """
        cdef float64_t entropy = 0.0
        cdef float64_t count_k
        cdef intp_t k
        cdef intp_t c

        for k in range(self.n_outputs):
            for c in range(self.n_classes[k]):
                count_k = self.sum_total[k, c]
                if count_k > 0.0:
                    count_k /= self.weighted_n_node_samples
                    entropy -= count_k * log(count_k)

        return entropy / self.n_outputs

    cdef void children_impurity(self, float64_t* impurity_left,
                                float64_t* impurity_right) noexcept nogil:
        """Evaluate the impurity in children nodes.

        i.e. the impurity of the left child (sample_indices[start:pos]) and the
        impurity the right child (sample_indices[pos:end]).

        Parameters
        ----------
        impurity_left : float64_t pointer
            The memory address to save the impurity of the left node
        impurity_right : float64_t pointer
            The memory address to save the impurity of the right node
        """
        cdef float64_t entropy_left = 0.0
        cdef float64_t entropy_right = 0.0
        cdef float64_t count_k
        cdef intp_t k
        cdef intp_t c

        for k in range(self.n_outputs):
            for c in range(self.n_classes[k]):
                count_k = self.sum_left[k, c]
                if count_k > 0.0:
                    count_k /= self.weighted_n_left
                    entropy_left -= count_k * log(count_k)

                count_k = self.sum_right[k, c]
                if count_k > 0.0:
                    count_k /= self.weighted_n_right
                    entropy_right -= count_k * log(count_k)

        impurity_left[0] = entropy_left / self.n_outputs
        impurity_right[0] = entropy_right / self.n_outputs


cdef class Gini(ClassificationCriterion):
    r"""Gini Index impurity criterion.

    This handles cases where the target is a classification taking values
    0, 1, ... K-2, K-1. If node m represents a region Rm with Nm observations,
    then let

        count_k = 1/ Nm \sum_{x_i in Rm} I(yi = k)

    be the proportion of class k observations in node m.

    The Gini Index is then defined as:

        index = \sum_{k=0}^{K-1} count_k (1 - count_k)
              = 1 - \sum_{k=0}^{K-1} count_k ** 2
    """

    cdef float64_t node_impurity(self) noexcept nogil:
        """Evaluate the impurity of the current node.

        Evaluate the Gini criterion as impurity of the current node,
        i.e. the impurity of sample_indices[start:end]. The smaller the impurity the
        better.
        """
        cdef float64_t gini = 0.0
        cdef float64_t sq_count
        cdef float64_t count_k
        cdef intp_t k
        cdef intp_t c

        for k in range(self.n_outputs):
            sq_count = 0.0

            for c in range(self.n_classes[k]):
                count_k = self.sum_total[k, c]
                sq_count += count_k * count_k

            gini += 1.0 - sq_count / (self.weighted_n_node_samples *
                                      self.weighted_n_node_samples)

        return gini / self.n_outputs

    cdef void children_impurity(self, float64_t* impurity_left,
                                float64_t* impurity_right) noexcept nogil:
        """Evaluate the impurity in children nodes.

        i.e. the impurity of the left child (sample_indices[start:pos]) and the
        impurity the right child (sample_indices[pos:end]) using the Gini index.

        Parameters
        ----------
        impurity_left : float64_t pointer
            The memory address to save the impurity of the left node to
        impurity_right : float64_t pointer
            The memory address to save the impurity of the right node to
        """
        cdef float64_t gini_left = 0.0
        cdef float64_t gini_right = 0.0
        cdef float64_t sq_count_left
        cdef float64_t sq_count_right
        cdef float64_t count_k
        cdef intp_t k
        cdef intp_t c

        for k in range(self.n_outputs):
            sq_count_left = 0.0
            sq_count_right = 0.0

            for c in range(self.n_classes[k]):
                count_k = self.sum_left[k, c]
                sq_count_left += count_k * count_k

                count_k = self.sum_right[k, c]
                sq_count_right += count_k * count_k

            gini_left += 1.0 - sq_count_left / (self.weighted_n_left *
                                                self.weighted_n_left)

            gini_right += 1.0 - sq_count_right / (self.weighted_n_right *
                                                  self.weighted_n_right)

        impurity_left[0] = gini_left / self.n_outputs
        impurity_right[0] = gini_right / self.n_outputs


cdef inline void _move_sums_regression(
    RegressionCriterion criterion,
    float64_t[::1] sum_1,
    float64_t[::1] sum_2,
    float64_t* weighted_n_1,
    float64_t* weighted_n_2,
    bint put_missing_in_1,
) noexcept nogil:
    """Distribute sum_total and sum_missing into sum_1 and sum_2.

    If there are missing values and:
    - put_missing_in_1 is True, then missing values to go sum_1. Specifically:
        sum_1 = sum_missing
        sum_2 = sum_total - sum_missing

    - put_missing_in_1 is False, then missing values go to sum_2. Specifically:
        sum_1 = 0
        sum_2 = sum_total
    """
    cdef:
        intp_t i
        intp_t n_bytes = criterion.n_outputs * sizeof(float64_t)
        bint has_missing = criterion.n_missing != 0

    if has_missing and put_missing_in_1:
        memcpy(&sum_1[0], &criterion.sum_missing[0], n_bytes)
        for i in range(criterion.n_outputs):
            sum_2[i] = criterion.sum_total[i] - criterion.sum_missing[i]
        weighted_n_1[0] = criterion.weighted_n_missing
        weighted_n_2[0] = criterion.weighted_n_node_samples - criterion.weighted_n_missing
    else:
        memset(&sum_1[0], 0, n_bytes)
        # Assigning sum_2 = sum_total for all outputs.
        memcpy(&sum_2[0], &criterion.sum_total[0], n_bytes)
        weighted_n_1[0] = 0.0
        weighted_n_2[0] = criterion.weighted_n_node_samples


cdef class RegressionCriterion(Criterion):
    r"""Abstract regression criterion.

    This handles cases where the target is a continuous value, and is
    evaluated by computing the variance of the target values left and right
    of the split point. The computation takes linear time with `n_samples`
    by using ::

        var = \sum_i^n (y_i - y_bar) ** 2
            = (\sum_i^n y_i ** 2) - n_samples * y_bar ** 2
    """

    def __cinit__(self, intp_t n_outputs, intp_t n_samples):
        """Initialize parameters for this criterion.

        Parameters
        ----------
        n_outputs : intp_t
            The number of targets to be predicted

        n_samples : intp_t
            The total number of samples to fit on
        """
        # Default values
        self.start = 0
        self.pos = 0
        self.end = 0

        self.n_outputs = n_outputs
        self.n_samples = n_samples
        self.n_node_samples = 0
        self.weighted_n_node_samples = 0.0
        self.weighted_n_left = 0.0
        self.weighted_n_right = 0.0
        self.weighted_n_missing = 0.0

        self.sq_sum_total = 0.0

        self.sum_total = np.zeros(n_outputs, dtype=np.float64)
        self.sum_left = np.zeros(n_outputs, dtype=np.float64)
        self.sum_right = np.zeros(n_outputs, dtype=np.float64)

    def __reduce__(self):
        return (type(self), (self.n_outputs, self.n_samples), self.__getstate__())

    cdef int init(
        self,
        const float64_t[:, ::1] y,
        const float64_t[:] sample_weight,
        float64_t weighted_n_samples,
        const intp_t[:] sample_indices,
        intp_t start,
        intp_t end,
    ) except -1 nogil:
        """Initialize the criterion.

        This initializes the criterion at node sample_indices[start:end] and children
        sample_indices[start:start] and sample_indices[start:end].
        """
        # Initialize fields
        self.y = y
        self.sample_weight = sample_weight
        self.sample_indices = sample_indices
        self.start = start
        self.end = end
        self.n_node_samples = end - start
        self.weighted_n_samples = weighted_n_samples
        self.weighted_n_node_samples = 0.

        cdef intp_t i
        cdef intp_t p
        cdef intp_t k
        cdef float64_t y_ik
        cdef float64_t w_y_ik
        cdef float64_t w = 1.0
        self.sq_sum_total = 0.0
        memset(&self.sum_total[0], 0, self.n_outputs * sizeof(float64_t))

        for p in range(start, end):
            i = sample_indices[p]

            if sample_weight is not None:
                w = sample_weight[i]

            for k in range(self.n_outputs):
                y_ik = self.y[i, k]
                w_y_ik = w * y_ik
                self.sum_total[k] += w_y_ik
                self.sq_sum_total += w_y_ik * y_ik

            self.weighted_n_node_samples += w

        # Reset to pos=start
        self.reset()
        return 0

    cdef void init_sum_missing(self):
        """Init sum_missing to hold sums for missing values."""
        self.sum_missing = np.zeros(self.n_outputs, dtype=np.float64)

    cdef void init_missing(self, intp_t n_missing) noexcept nogil:
        """Initialize sum_missing if there are missing values.

        This method assumes that caller placed the missing samples in
        self.sample_indices[-n_missing:]
        """
        cdef intp_t i, p, k
        cdef float64_t y_ik
        cdef float64_t w_y_ik
        cdef float64_t w = 1.0

        self.n_missing = n_missing
        if n_missing == 0:
            return

        memset(&self.sum_missing[0], 0, self.n_outputs * sizeof(float64_t))

        self.weighted_n_missing = 0.0

        # The missing samples are assumed to be in self.sample_indices[-n_missing:]
        for p in range(self.end - n_missing, self.end):
            i = self.sample_indices[p]
            if self.sample_weight is not None:
                w = self.sample_weight[i]

            for k in range(self.n_outputs):
                y_ik = self.y[i, k]
                w_y_ik = w * y_ik
                self.sum_missing[k] += w_y_ik

            self.weighted_n_missing += w

    cdef int reset(self) except -1 nogil:
        """Reset the criterion at pos=start."""
        self.pos = self.start
        _move_sums_regression(
            self,
            self.sum_left,
            self.sum_right,
            &self.weighted_n_left,
            &self.weighted_n_right,
            self.missing_go_to_left
        )
        return 0

    cdef int reverse_reset(self) except -1 nogil:
        """Reset the criterion at pos=end."""
        self.pos = self.end
        _move_sums_regression(
            self,
            self.sum_right,
            self.sum_left,
            &self.weighted_n_right,
            &self.weighted_n_left,
            not self.missing_go_to_left
        )
        return 0

    cdef int update(self, intp_t new_pos) except -1 nogil:
        """Updated statistics by moving sample_indices[pos:new_pos] to the left."""
        cdef intp_t pos = self.pos

        # The missing samples are assumed to be in
        # self.sample_indices[-self.n_missing:] that is
        # self.sample_indices[end_non_missing:self.end].
        cdef intp_t end_non_missing = self.end - self.n_missing
        cdef intp_t i
        cdef intp_t p
        cdef intp_t k
        cdef float64_t w = 1.0

        # Update statistics up to new_pos
        #
        # Given that
        #           sum_left[x] +  sum_right[x] = sum_total[x]
        # and that sum_total is known, we are going to update
        # sum_left from the direction that require the least amount
        # of computations, i.e. from pos to new_pos or from end to new_pos.
        if (new_pos - pos) <= (end_non_missing - new_pos):
            for p in range(pos, new_pos):
                i = self.sample_indices[p]

                if self.sample_weight is not None:
                    w = self.sample_weight[i]

                for k in range(self.n_outputs):
                    self.sum_left[k] += w * self.y[i, k]

                self.weighted_n_left += w
        else:
            self.reverse_reset()

            for p in range(end_non_missing - 1, new_pos - 1, -1):
                i = self.sample_indices[p]

                if self.sample_weight is not None:
                    w = self.sample_weight[i]

                for k in range(self.n_outputs):
                    self.sum_left[k] -= w * self.y[i, k]

                self.weighted_n_left -= w

        self.weighted_n_right = (self.weighted_n_node_samples -
                                 self.weighted_n_left)
        for k in range(self.n_outputs):
            self.sum_right[k] = self.sum_total[k] - self.sum_left[k]

        self.pos = new_pos
        return 0

    cdef float64_t node_impurity(self) noexcept nogil:
        pass

    cdef void children_impurity(self, float64_t* impurity_left,
                                float64_t* impurity_right) noexcept nogil:
        pass

    cdef void node_value(self, float64_t* dest) noexcept nogil:
        """Compute the node value of sample_indices[start:end] into dest."""
        cdef intp_t k

        for k in range(self.n_outputs):
            dest[k] = self.sum_total[k] / self.weighted_n_node_samples

    cdef inline void clip_node_value(self, float64_t* dest, float64_t lower_bound, float64_t upper_bound) noexcept nogil:
        """Clip the value in dest between lower_bound and upper_bound for monotonic constraints."""
        if dest[0] < lower_bound:
            dest[0] = lower_bound
        elif dest[0] > upper_bound:
            dest[0] = upper_bound

    cdef float64_t middle_value(self) noexcept nogil:
        """Compute the middle value of a split for monotonicity constraints as the simple average
        of the left and right children values.

        Monotonicity constraints are only supported for single-output trees we can safely assume
        n_outputs == 1.
        """
        return (
            (self.sum_left[0] / (2 * self.weighted_n_left)) +
            (self.sum_right[0] / (2 * self.weighted_n_right))
        )

    cdef bint check_monotonicity(
        self,
        cnp.int8_t monotonic_cst,
        float64_t lower_bound,
        float64_t upper_bound,
    ) noexcept nogil:
        """Check monotonicity constraint is satisfied at the current regression split"""
        cdef:
            float64_t value_left = self.sum_left[0] / self.weighted_n_left
            float64_t value_right = self.sum_right[0] / self.weighted_n_right

        return self._check_monotonicity(monotonic_cst, lower_bound, upper_bound, value_left, value_right)


cdef class MSE(RegressionCriterion):
    """Mean squared error impurity criterion.

        MSE = var_left + var_right
    """

    cdef float64_t node_impurity(self) noexcept nogil:
        """Evaluate the impurity of the current node.

        Evaluate the MSE criterion as impurity of the current node,
        i.e. the impurity of sample_indices[start:end]. The smaller the impurity the
        better.
        """
        cdef float64_t impurity
        cdef intp_t k

        impurity = self.sq_sum_total / self.weighted_n_node_samples
        for k in range(self.n_outputs):
            impurity -= (self.sum_total[k] / self.weighted_n_node_samples)**2.0

        return impurity / self.n_outputs

    cdef float64_t proxy_impurity_improvement(self) noexcept nogil:
        """Compute a proxy of the impurity reduction.

        This method is used to speed up the search for the best split.
        It is a proxy quantity such that the split that maximizes this value
        also maximizes the impurity improvement. It neglects all constant terms
        of the impurity decrease for a given split.

        The absolute impurity improvement is only computed by the
        impurity_improvement method once the best split has been found.

        The MSE proxy is derived from

            sum_{i left}(y_i - y_pred_L)^2 + sum_{i right}(y_i - y_pred_R)^2
            = sum(y_i^2) - n_L * mean_{i left}(y_i)^2 - n_R * mean_{i right}(y_i)^2

        Neglecting constant terms, this gives:

            - 1/n_L * sum_{i left}(y_i)^2 - 1/n_R * sum_{i right}(y_i)^2
        """
        cdef intp_t k
        cdef float64_t proxy_impurity_left = 0.0
        cdef float64_t proxy_impurity_right = 0.0

        for k in range(self.n_outputs):
            proxy_impurity_left += self.sum_left[k] * self.sum_left[k]
            proxy_impurity_right += self.sum_right[k] * self.sum_right[k]

        return (proxy_impurity_left / self.weighted_n_left +
                proxy_impurity_right / self.weighted_n_right)

    cdef void children_impurity(self, float64_t* impurity_left,
                                float64_t* impurity_right) noexcept nogil:
        """Evaluate the impurity in children nodes.

        i.e. the impurity of the left child (sample_indices[start:pos]) and the
        impurity the right child (sample_indices[pos:end]).
        """
        cdef const float64_t[:] sample_weight = self.sample_weight
        cdef const intp_t[:] sample_indices = self.sample_indices
        cdef intp_t pos = self.pos
        cdef intp_t start = self.start

        cdef float64_t y_ik

        cdef float64_t sq_sum_left = 0.0
        cdef float64_t sq_sum_right

        cdef intp_t i
        cdef intp_t p
        cdef intp_t k
        cdef float64_t w = 1.0

        cdef intp_t end_non_missing

        for p in range(start, pos):
            i = sample_indices[p]

            if sample_weight is not None:
                w = sample_weight[i]

            for k in range(self.n_outputs):
                y_ik = self.y[i, k]
                sq_sum_left += w * y_ik * y_ik

        if self.missing_go_to_left:
            # add up the impact of these missing values on the left child
            # statistics.
            # Note: this only impacts the square sum as the sum
            # is modified elsewhere.
            end_non_missing = self.end - self.n_missing

            for p in range(end_non_missing, self.end):
                i = sample_indices[p]
                if sample_weight is not None:
                    w = sample_weight[i]

                for k in range(self.n_outputs):
                    y_ik = self.y[i, k]
                    sq_sum_left += w * y_ik * y_ik

        sq_sum_right = self.sq_sum_total - sq_sum_left

        impurity_left[0] = sq_sum_left / self.weighted_n_left
        impurity_right[0] = sq_sum_right / self.weighted_n_right

        for k in range(self.n_outputs):
            impurity_left[0] -= (self.sum_left[k] / self.weighted_n_left) ** 2.0
            impurity_right[0] -= (self.sum_right[k] / self.weighted_n_right) ** 2.0

        impurity_left[0] /= self.n_outputs
        impurity_right[0] /= self.n_outputs


# Helper for Pinball criterion:

cdef void precompute_pinball_losses(
    const float64_t[::1] sorted_y,
    const intp_t[::1] ranks,
    const float64_t[:] sample_weight,
    const intp_t[:] sample_indices,
    WeightedFenwickTree tree,
    intp_t start,
    intp_t end,
    float64_t alpha,
    float64_t[::1] pinball_losses,
    float64_t[::1] quantiles,
) noexcept nogil:
    """
    Fill `pinball_losses` and `quantiles`.

    If start < end:
        Forward pass: Computes the "prefix" losses/quantiles
        i.e the losses for each set of indices sample_indices[start:start + i]
        with i in {1, ..., n}, where n = end - start.
    Else:
        Backward pass: Computes the "suffix" losses/quantiles
        i.e the losses for each set of indices sample_indices[start - i:start]
        with i in {1, ..., n}, where n = start - end.

    Parameters
    ----------
    sorted_y : const float64_t[::1]
        Target values, sorted
    ranks : const intp_t[::1]
        Ranks of the node-local values of y for points in sample_indices such that:
        sorted_y[ranks[p]] == y[sample_indices[p]] for any p in [start, end) or
        (end, start].
    sample_weight : const float64_t[:]
    sample_indices : const intp_t[:]
        indices indicating which samples to use. Shape: (n_samples,)
    tree : WeightedFenwickTree
        pre-instanciated tree
    start : intp_t
        Start index in `sample_indices`
    end : intp_t
        End index (exclusive) in `sample_indices`
    alpha : float64_t
        Quantile level for the pinball loss (between 0 and 1)
    pinball_losses : float64_t[::1]
        array to store (increment) the computed pinball losses. Shape: (n,)
        with n := end - start
    quantiles : float64_t[::1]
        array to store (overwrite) the computed quantiles. Shape: (n,)

    Complexity: O(n log n)
    """
    cdef:
        intp_t p, i, step, n, rank, quantile_rank, quantile_prev_rank
        float64_t w = 1.
        float64_t target_weight, quantile
        float64_t w_right, w_left, wy_left, wy_right

    if start < end:
        step = 1
        n = end - start
    else:
        n = start - end
        step = -1

    tree.reset(n)

    p = start
    # We iterate exactly `n` samples starting at absolute index `start` and
    # move by `step` (+1 for the forward pass, -1 for the backward pass).
    for _ in range(n):
        i = sample_indices[p]
        if sample_weight is not None:
            w = sample_weight[i]
        # Activate sample i at its rank:
        rank = ranks[p]
        tree.add(rank, sorted_y[rank], w)

        # Weighted quantile by cumulative weight: the quantile is where the
        # cumulative weight crosses alpha of the total weight.
        target_weight = alpha * tree.total_w
        # find the smallest activated rank with cumulative weight > target_weight
        # while returning the prefix sums (`w_left` and `wy_left`)
        # up to (and excluding) that index:
        quantile_rank = tree.search(target_weight, &w_left, &wy_left, &quantile_prev_rank)

        if quantile_rank != quantile_prev_rank:
            # Exact match for target_weight fell between two consecutive ranks:
            # cumulative weight up to `quantile_rank` excluded is exactly target_weight.
            # In that case, `quantile_prev_rank` is the activated rank such that
            # the cumulative weight up to it included is exactly target_weight.
            # In this case we take the mid-point to match with
            # sklearn.utils.stats._weighted_percentile(..., average=True)
            quantile = (sorted_y[quantile_prev_rank] + sorted_y[quantile_rank]) / 2
        else:
            # if there are no exact match for target_weight in the cumulative weights
            # `quantile_rank == quantile_prev_rank` and the quantile is:
            quantile = sorted_y[quantile_rank]

        # Convert left prefix sums into right-hand complements.
        w_right = tree.total_w - w_left
        wy_right = tree.total_wy - wy_left

        quantiles[p] = quantile
        # Pinball loss identity for the alpha-quantile at the current set:
        #   sum_{y_i >= q} w_i * alpha * (y_i - q) = alpha * (wy_right - q * w_right)
        #   sum_{y_i <  q} w_i * (1-alpha) * (q - y_i) = (1-alpha) * (q * w_left - wy_left)
        pinball_losses[p] += (
            alpha * (wy_right - quantile * w_right)
            + (1.0 - alpha) * (quantile * w_left - wy_left)
        )
        p += step


cdef inline void compute_ranks(
    float64_t* sorted_y,
    intp_t* sorted_indices,
    intp_t* ranks,
    intp_t n
) noexcept nogil:
    """Sort `sorted_y` inplace and fill `ranks` accordingly"""
    cdef intp_t i
    for i in range(n):
        sorted_indices[i] = i
    sort(sorted_y, sorted_indices, n)
    for i in range(n):
        ranks[sorted_indices[i]] = i


def _py_precompute_pinball_losses(
    const float64_t[:, ::1] ys,
    const float64_t[:] sample_weight,
    const intp_t[:] sample_indices,
    const intp_t start,
    const intp_t end,
    const intp_t n,
    const float64_t alpha,
):
    """Used for testing precompute_pinball_losses."""
    cdef:
        intp_t p, i
        intp_t s = start
        intp_t e = end
        WeightedFenwickTree tree = WeightedFenwickTree(n)
        float64_t[::1] sorted_y = np.empty(n, dtype=np.float64)
        intp_t[::1] sorted_indices = np.empty(n, dtype=np.intp)
        intp_t[::1] ranks = np.empty(n, dtype=np.intp)
        float64_t[::1] pinball_losses = np.zeros(n, dtype=np.float64)
        float64_t[::1] quantiles = np.empty(n, dtype=np.float64)

    if start > end:
        s = end + 1
        e = start + 1
    for p in range(s, e):
        i = sample_indices[p]
        sorted_y[p - s] = ys[i, 0]
    compute_ranks(&sorted_y[0], &sorted_indices[0], &ranks[s], n)

    precompute_pinball_losses(
        sorted_y, ranks, sample_weight, sample_indices, tree,
        start, end, alpha, pinball_losses, quantiles
    )
    return np.asarray(pinball_losses)[s:e], np.asarray(quantiles)[s:e]


cdef class Pinball(Criterion):
    r"""Pinball loss impurity criterion.

    This criterion generalizes the Mean Absolute Error (MAE) by using
    quantile regression with a specified quantile level (alpha).
    MAE corresponds to alpha=0.5 (median).

    Pinball loss = (1 / n)*(\sum_i rho_alpha(y_i - q_i)), where y_i is the true
    value, q_i is the predicted quantile, and rho_alpha is the pinball loss function:
        rho_alpha(u) = u * (alpha - I(u < 0))
    In a decision tree, that prediction is the (weighted) alpha-quantile
    of the targets in the node.

    How this implementation works
    -----------------------------
    This class precomputes in `reset`, for the current node,
    the pinball loss values and corresponding quantiles for all
    potential split positions: every p in [start, end).

    For that:
    - We first compute the rank of each samples node-local sorted order of target values.
      `self.ranks[p]` gives the rank of sample p.
    - While iterating the segment of indices (p in [start, end)), we
        * "activate" one sample at a time at its rank within a prefix sum tree,
          the `WeightedFenwickTree`: `tree.add(rank, y, weight)`
          The tree maintains cumulative sums of weights and of `weight * y`
        * search for the target weight alpha * total_weight in the tree:
          `tree.search(alpha * current_total_weight)`.
          This allows us to retrieve/compute:
            * the current weighted quantile value
            * the pinball loss contribution via the standard pinball-loss identity:
              PL = alpha * (wy_right - quantile * w_right) + (1-alpha) * (quantile * w_left - wy_left)
    - We perform two such passes:
        * one forward from `start` to `end - 1` to fill `left_pinball_losses[p]` and
          `left_quantiles[p]` for left children.
        * one backward from `end - 1` down to `start` to fill
          `right_pinball_losses[p]` and `right_quantiles[p]` for right children.

    Complexity: time complexity is O(n log n), indeed:
    - computing ranks is based on sorting: O(n log n)
    - add and search operations in the Fenwick tree are O(log n).
      => the forward and backward passes are O(n log n).

    How the other methods use the precomputations
    --------------------------------------------
    - `reset` performs the precomputation described above.
      It also stores the node weighted quantile per output in
      `node_quantiles` (prediction value of the node).

    - `update(new_pos)` only updates `weighted_n_left` and `weighted_n_right`;
      no recomputation of losses is needed.

    - `children_impurity` reads the precomputed pinball losses at
      `left_pinball_losses[pos - 1]` and `right_pinball_losses[pos]` and scales
      them by the corresponding child weights and `n_outputs` to report the
      impurity of each child.

    - `middle_value` and `check_monotonicity` use the precomputed
      `left_quantiles[pos - 1]` and `right_quantiles[pos]` to derive the
      mid-point value and to validate monotonic constraints when enabled.

    - Missing values are not supported for Pinball: `init_missing` raises.

    For a complementary, in-depth discussion of the mathematics and design
    choices, see the external report:
    https://github.com/cakedev0/fast-mae-split/blob/main/report.ipynb
    """
    cdef float64_t alpha
    cdef float64_t[::1] node_quantiles
    cdef float64_t[::1] left_pinball_losses
    cdef float64_t[::1] right_pinball_losses
    cdef float64_t[::1] left_quantiles
    cdef float64_t[::1] right_quantiles
    cdef float64_t[::1] sorted_y
    cdef intp_t [::1] sorted_indices
    cdef intp_t[::1] ranks
    cdef WeightedFenwickTree prefix_sum_tree

    def __cinit__(self, intp_t n_outputs, intp_t n_samples, float64_t alpha=0.5):
        """Initialize parameters for this criterion.

        Parameters
        ----------
        n_outputs : intp_t
            The number of targets to be predicted

        n_samples : intp_t
            The total number of samples to fit on
        """
        self.alpha = alpha

        # Default values
        self.start = 0
        self.pos = 0
        self.end = 0

        self.n_outputs = n_outputs
        self.n_samples = n_samples
        self.n_node_samples = 0
        self.weighted_n_node_samples = 0.0
        self.weighted_n_left = 0.0
        self.weighted_n_right = 0.0

        self.node_quantiles = np.zeros(n_outputs, dtype=np.float64)

        # Note: this criterion has a  n_samples x 64 bytes memory footprint, which is
        # fine as it's instantiated only once to build an entire tree
        self.left_pinball_losses = np.empty(n_samples, dtype=np.float64)
        self.right_pinball_losses = np.empty(n_samples, dtype=np.float64)
        self.left_quantiles = np.empty(n_samples, dtype=np.float64)
        self.right_quantiles = np.empty(n_samples, dtype=np.float64)
        self.ranks = np.empty(n_samples, dtype=np.intp)
        # Important: The arrays declared above are indexed with
        # the absolute position `p` in `sample_indices` (not with a 0-based offset).
        # The forward and backward passes in `reset` method ensure that
        # for any current split position `pos` we can read:
        # - left child precomputed values at `p = pos - 1`, and
        # - right child precomputed values at `p = pos`.

        self.prefix_sum_tree = WeightedFenwickTree(n_samples)
        # used memory: 2 float64 arrays of size n_samples + 1
        # we reuse a single `WeightedFenwickTree` instance to build prefix
        # and suffix aggregates over the node samples.

        # Work buffer arrays, used with 0-based offset:
        self.sorted_y = np.empty(n_samples, dtype=np.float64)
        self.sorted_indices = np.empty(n_samples, dtype=np.intp)

    cdef int init(
        self,
        const float64_t[:, ::1] y,
        const float64_t[:] sample_weight,
        float64_t weighted_n_samples,
        const intp_t[:] sample_indices,
        intp_t start,
        intp_t end,
    ) except -1 nogil:
        """Initialize the criterion.

        This initializes the criterion at node sample_indices[start:end] and children
        sample_indices[start:start] and sample_indices[start:end].

        WARNING: sample_indices will be modified in-place externally
        after this method is called
        """
        cdef:
            intp_t i, p
            intp_t n = end - start
            float64_t w = 1.0
        # Initialize fields
        self.y = y
        self.sample_weight = sample_weight
        self.sample_indices = sample_indices
        self.start = start
        self.end = end
        self.n_node_samples = n
        self.weighted_n_samples = weighted_n_samples
        self.weighted_n_node_samples = 0.

        for p in range(start, end):
            i = sample_indices[p]
            if sample_weight is not None:
                w = sample_weight[i]
            self.weighted_n_node_samples += w

        # Reset to pos=start
        self.reset()
        return 0

    cdef void init_missing(self, intp_t n_missing) noexcept nogil:
        """Raise error if n_missing != 0."""
        if n_missing == 0:
            return
        with gil:
            raise ValueError("missing values is not supported for Pinball criterion.")

    cdef int reset(self) except -1 nogil:
        """Reset the criterion at pos=start.

        Returns -1 in case of failure to allocate memory (and raise MemoryError)
        or 0 otherwise.

        Reset might be called after an external class has changed
        inplace self.sample_indices[start:end], hence re-computing
        the pinball loss is needed.
        """
        cdef intp_t k, p, i

        self.weighted_n_left = 0.0
        self.weighted_n_right = self.weighted_n_node_samples
        self.pos = self.start

        n_bytes = self.n_node_samples * sizeof(float64_t)
        memset(&self.left_pinball_losses[self.start],  0, n_bytes)
        memset(&self.right_pinball_losses[self.start], 0, n_bytes)

        # Multi-output handling:
        # pinball losses are accumulated across outputs by
        # incrementing `left_pinball_losses` and `right_pinball_losses` on each pass.
        # The per-output quantiles arrays are overwritten at each output iteration
        # as they are only used for monotonicity checks when `n_outputs == 1`.

        # Precompute pinball losses (summed over each output)
        # and quantiles (used only when n_outputs=1)
        # of the right and left child of all possible splits
        # for the current ordering of `sample_indices`
        # Precomputation is needed here and can't be done step-by-step in the update method
        # like for other criterions. Indeed, we don't have efficient ways to update right child
        # statistics when removing samples from it. So we compute right child losses/quantiles by
        # traversing from right to left (and hence only adding samples).
        for k in range(self.n_outputs):

            # 1) Node-local ordering:
            # for each output k, the values `y[sample_indices[p], k]` for p
            # in [start, end) are copied into self.sorted_y[0:n_node_samples]`
            # and ranked with `compute_ranks`.
            # The resulting `self.ranks[p]` gives the rank of sample p in the
            # node-local sorted order.
            for p in range(self.start, self.end):
                i = self.sample_indices[p]
                self.sorted_y[p - self.start] = self.y[i, k]

            compute_ranks(
                &self.sorted_y[0],
                &self.sorted_indices[0],
                &self.ranks[self.start],
                self.n_node_samples,
            )

            # 2) Forward pass
            # from `start` to `end - 1` to fill `left_pinball_losses[p]` and
            # `left_quantiles[p]` for left children.
            precompute_pinball_losses(
                self.sorted_y, self.ranks, self.sample_weight, self.sample_indices,
                self.prefix_sum_tree, self.start, self.end, self.alpha,
                # left_pinball_losses is incremented, left_quantiles is overwritten
                self.left_pinball_losses, self.left_quantiles
            )
            # 3) Backward pass
            # from `end - 1` down to `start` to fill `right_pinball_losses[p]`
            # and `right_quantiles[p]` for right children.
            precompute_pinball_losses(
                self.sorted_y, self.ranks, self.sample_weight, self.sample_indices,
                self.prefix_sum_tree, self.end - 1, self.start - 1, self.alpha,
                # right_pinball_losses is incremented, right_quantiles is overwritten
                self.right_pinball_losses, self.right_quantiles
            )

            # Store the quantile for the current node: when p == self.start all the
            # node's data points are sent to the right child, so the current node
            # quantile value and the right child quantile value would be equal.
            self.node_quantiles[k] = self.right_quantiles[self.start]

        return 0

    cdef int reverse_reset(self) except -1 nogil:
        """For this class, this method is never called."""
        raise NotImplementedError("This method is not implemented for this subclass")

    cdef int update(self, intp_t new_pos) except -1 nogil:
        """Updated statistics by moving sample_indices[pos:new_pos] to the left.
        new_pos is guaranteed to be greater than pos.

        Returns -1 in case of failure to allocate memory (and raise MemoryError)
        or 0 otherwise.

        Time complexity: O(new_pos - pos) (which usually is O(1), at least for dense data).
        """
        cdef intp_t pos = self.pos
        cdef intp_t i, p
        cdef float64_t w = 1.0

        # Update statistics up to new_pos
        for p in range(pos, new_pos):
            i = self.sample_indices[p]
            if self.sample_weight is not None:
                w = self.sample_weight[i]
            self.weighted_n_left += w

        self.weighted_n_right = (self.weighted_n_node_samples -
                                 self.weighted_n_left)
        self.pos = new_pos
        return 0

    cdef void node_value(self, float64_t* dest) noexcept nogil:
        """Computes the node value of sample_indices[start:end] into dest."""
        cdef intp_t k
        for k in range(self.n_outputs):
            dest[k] = <float64_t> self.node_quantiles[k]

    cdef inline float64_t middle_value(self) noexcept nogil:
        """Compute the middle value of a split for monotonicity constraints as the simple average
        of the left and right children values.

        Monotonicity constraints are only supported for single-output trees we can safely assume
        n_outputs == 1.
        """
        return (
            self.left_quantiles[self.pos - 1]
            + self.right_quantiles[self.pos]
        ) / 2

    cdef inline bint check_monotonicity(
        self,
        cnp.int8_t monotonic_cst,
        float64_t lower_bound,
        float64_t upper_bound,
    ) noexcept nogil:
        """Check monotonicity constraint is satisfied at the current regression split"""
        return self._check_monotonicity(
            monotonic_cst, lower_bound, upper_bound,
            self.left_quantiles[self.pos - 1], self.right_quantiles[self.pos])

    cdef float64_t node_impurity(self) noexcept nogil:
        """Evaluate the impurity of the current node.

        Evaluate the Pinball criterion as impurity of the current node,
        i.e. the impurity of sample_indices[start:end]. The smaller the impurity the
        better.

        Time complexity: O(1) (precomputed in `.reset()`)
        """
        return (
            self.right_pinball_losses[self.start]
            / (self.weighted_n_node_samples * self.n_outputs)
        )

    cdef void children_impurity(self, float64_t* p_impurity_left,
                                float64_t* p_impurity_right) noexcept nogil:
        """Evaluate the impurity in children nodes.

        i.e. the impurity of the left child (sample_indices[start:pos]) and the
        impurity the right child (sample_indices[pos:end]).

        Time complexity: O(1) (precomputed in `.reset()`)
        """
        cdef float64_t impurity_left = 0.0
        cdef float64_t impurity_right = 0.0

        # if pos == start, left child is empty, hence impurity is 0
        if self.pos > self.start:
            impurity_left += self.left_pinball_losses[self.pos - 1]
        p_impurity_left[0] = impurity_left / (self.weighted_n_left *
                                              self.n_outputs)

        # if pos == end, right child is empty, hence impurity is 0
        if self.pos < self.end:
            impurity_right += self.right_pinball_losses[self.pos]
        p_impurity_right[0] = impurity_right / (self.weighted_n_right *
                                                self.n_outputs)

    def __reduce__(self):
        return (type(self), (self.n_outputs, self.n_samples, self.alpha), self.__getstate__())

    # this method is copied from the RegressionCriterion abstract class:
    cdef inline void clip_node_value(self, float64_t* dest, float64_t lower_bound, float64_t upper_bound) noexcept nogil:
        """Clip the value in dest between lower_bound and upper_bound for monotonic constraints."""
        if dest[0] < lower_bound:
            dest[0] = lower_bound
        elif dest[0] > upper_bound:
            dest[0] = upper_bound


cdef class MAE(Pinball):
    """
    The median is just the quantile alpha=0.5
    And the absolute error is twice the pinball_loss (with alpha=0.5)
    """

    # XXX: Trust the instanciater to pass alpha=0.5 to the __cinit__...

    cdef float64_t node_impurity(self) noexcept nogil:
        return 2 * Pinball.node_impurity(self)

    cdef void children_impurity(self, float64_t* p_impurity_left,
                                float64_t* p_impurity_right) noexcept nogil:
        Pinball.children_impurity(self, p_impurity_left, p_impurity_right)
        p_impurity_left[0] *= 2
        p_impurity_right[0] *= 2

    def __reduce__(self):
        return (type(self), (self.n_outputs, self.n_samples), self.__getstate__())


cdef class Poisson(RegressionCriterion):
    """Half Poisson deviance as impurity criterion.

    Poisson deviance = 2/n * sum(y_true * log(y_true/y_pred) + y_pred - y_true)

    Note that the deviance is >= 0, and since we have `y_pred = mean(y_true)`
    at the leaves, one always has `sum(y_pred - y_true) = 0`. It remains the
    implemented impurity (factor 2 is skipped):
        1/n * sum(y_true * log(y_true/y_pred)
    """
    # FIXME in 1.0:
    # min_impurity_split with default = 0 forces us to use a non-negative
    # impurity like the Poisson deviance. Without this restriction, one could
    # throw away the 'constant' term sum(y_true * log(y_true)) and just use
    # Poisson loss = - 1/n * sum(y_true * log(y_pred))
    #              = - 1/n * sum(y_true * log(mean(y_true))
    #              = - mean(y_true) * log(mean(y_true))
    # With this trick (used in proxy_impurity_improvement()), as for MSE,
    # children_impurity would only need to go over left xor right split, not
    # both. This could be faster.

    cdef float64_t node_impurity(self) noexcept nogil:
        """Evaluate the impurity of the current node.

        Evaluate the Poisson criterion as impurity of the current node,
        i.e. the impurity of sample_indices[start:end]. The smaller the impurity the
        better.
        """
        return self.poisson_loss(self.start, self.end, self.sum_total,
                                 self.weighted_n_node_samples)

    cdef float64_t proxy_impurity_improvement(self) noexcept nogil:
        """Compute a proxy of the impurity reduction.

        This method is used to speed up the search for the best split.
        It is a proxy quantity such that the split that maximizes this value
        also maximizes the impurity improvement. It neglects all constant terms
        of the impurity decrease for a given split.

        The absolute impurity improvement is only computed by the
        impurity_improvement method once the best split has been found.

        The Poisson proxy is derived from:

              sum_{i left }(y_i * log(y_i / y_pred_L))
            + sum_{i right}(y_i * log(y_i / y_pred_R))
            = sum(y_i * log(y_i) - n_L * mean_{i left}(y_i) * log(mean_{i left}(y_i))
                                 - n_R * mean_{i right}(y_i) * log(mean_{i right}(y_i))

        Neglecting constant terms, this gives

            - sum{i left }(y_i) * log(mean{i left}(y_i))
            - sum{i right}(y_i) * log(mean{i right}(y_i))
        """
        cdef intp_t k
        cdef float64_t proxy_impurity_left = 0.0
        cdef float64_t proxy_impurity_right = 0.0
        cdef float64_t y_mean_left = 0.
        cdef float64_t y_mean_right = 0.

        for k in range(self.n_outputs):
            if (self.sum_left[k] <= EPSILON) or (self.sum_right[k] <= EPSILON):
                # Poisson loss does not allow non-positive predictions. We
                # therefore forbid splits that have child nodes with
                # sum(y_i) <= 0.
                # Since sum_right = sum_total - sum_left, it can lead to
                # floating point rounding error and will not give zero. Thus,
                # we relax the above comparison to sum(y_i) <= EPSILON.
                return -INFINITY
            else:
                y_mean_left = self.sum_left[k] / self.weighted_n_left
                y_mean_right = self.sum_right[k] / self.weighted_n_right
                proxy_impurity_left -= self.sum_left[k] * log(y_mean_left)
                proxy_impurity_right -= self.sum_right[k] * log(y_mean_right)

        return - proxy_impurity_left - proxy_impurity_right

    cdef void children_impurity(self, float64_t* impurity_left,
                                float64_t* impurity_right) noexcept nogil:
        """Evaluate the impurity in children nodes.

        i.e. the impurity of the left child (sample_indices[start:pos]) and the
        impurity of the right child (sample_indices[pos:end]) for Poisson.
        """
        cdef intp_t start = self.start
        cdef intp_t pos = self.pos
        cdef intp_t end = self.end

        impurity_left[0] = self.poisson_loss(start, pos, self.sum_left,
                                             self.weighted_n_left)

        impurity_right[0] = self.poisson_loss(pos, end, self.sum_right,
                                              self.weighted_n_right)

    cdef inline float64_t poisson_loss(
        self,
        intp_t start,
        intp_t end,
        const float64_t[::1] y_sum,
        float64_t weight_sum
    ) noexcept nogil:
        """Helper function to compute Poisson loss (~deviance) of a given node.
        """
        cdef const float64_t[:, ::1] y = self.y
        cdef const float64_t[:] sample_weight = self.sample_weight
        cdef const intp_t[:] sample_indices = self.sample_indices

        cdef float64_t y_mean = 0.
        cdef float64_t poisson_loss = 0.
        cdef float64_t w = 1.0
        cdef intp_t i, k, p
        cdef intp_t n_outputs = self.n_outputs

        for k in range(n_outputs):
            if y_sum[k] <= EPSILON:
                # y_sum could be computed from the subtraction
                # sum_right = sum_total - sum_left leading to a potential
                # floating point rounding error.
                # Thus, we relax the comparison y_sum <= 0 to
                # y_sum <= EPSILON.
                return INFINITY

            y_mean = y_sum[k] / weight_sum

            for p in range(start, end):
                i = sample_indices[p]

                if sample_weight is not None:
                    w = sample_weight[i]

                poisson_loss += w * xlogy(y[i, k], y[i, k] / y_mean)
        return poisson_loss / (weight_sum * n_outputs)
