# cython: language_level=3
# cython: boundscheck=False, wraparound=False, initializedcheck=False, cdivision=True

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
#          Adam Li <adam2392@gmail.com>
#          Jong Shin <jshinm@gmail.com>
#
# License: BSD 3 clause

from libc.math cimport INFINITY, fabs
from libc.string cimport memcpy, memset

import numpy as np

cimport numpy as cnp

cnp.import_array()

from scipy.special.cython_special cimport xlogy

from ._utils cimport WeightedMedianCalculator, log


# EPSILON is used in the Poisson criterion
cdef double EPSILON = 10 * np.finfo('double').eps

cdef class BaseCriterion:
    """This is an abstract interface for criterion.

    For example, a tree model could
    be either supervisedly, or unsupervisedly computing impurity on samples of
    covariates, or labels, or both. Although scikit-learn currently only contains
    supervised tree methods, this class enables 3rd party packages to leverage
    scikit-learn's Cython code for criteria.

    The downstream classes _must_ implement methods to compute the impurity
    in current node and in children nodes.

    This object stores methods on how to calculate how good a split is using
    a set API.

    Samples in the "current" node are stored in `samples[start:end]` which is
    partitioned around `pos` (an index in `start:end`) so that:
       - the samples of left child node are stored in `samples[start:pos]`
       - the samples of right child node are stored in `samples[pos:end]`
    """
    def __getstate__(self):
        return {}

    def __setstate__(self, d):
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

    cdef int update(self, SIZE_t new_pos) except -1 nogil:
        """Updated statistics by moving sample_indices[pos:new_pos] to the left child.

        This updates the collected statistics by moving sample_indices[pos:new_pos]
        from the right child to the left child. It must be implemented by
        the subclass.

        Parameters
        ----------
        new_pos : SIZE_t
            New starting index position of the sample_indices in the right child
        """
        pass

    cdef double node_impurity(self) noexcept nogil:
        """Placeholder for calculating the impurity of the node.

        Placeholder for a method which will evaluate the impurity of
        the current node, i.e. the impurity of sample_indices[start:end]. This is the
        primary function of the criterion class. The smaller the impurity the
        better.
        """
        pass

    cdef void children_impurity(self, double* impurity_left,
                                double* impurity_right) noexcept nogil:
        """Placeholder for calculating the impurity of children.

        Placeholder for a method which evaluates the impurity in
        children nodes, i.e. the impurity of sample_indices[start:pos] + the impurity
        of sample_indices[pos:end].

        Parameters
        ----------
        impurity_left : double pointer
            The memory address where the impurity of the left child should be
            stored.
        impurity_right : double pointer
            The memory address where the impurity of the right child should be
            stored
        """
        pass

    cdef void node_value(self, double* dest) noexcept nogil:
        """Placeholder for storing the node value.

        Placeholder for a method which will compute the node value
        of sample_indices[start:end] and save the value into dest.

        Parameters
        ----------
        dest : double pointer
            The memory address where the node value should be stored.
        """
        pass

    cdef double proxy_impurity_improvement(self) noexcept nogil:
        """Compute a proxy of the impurity reduction.

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

    cdef double impurity_improvement(self, double impurity_parent,
                                     double impurity_left,
                                     double impurity_right) noexcept nogil:
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
        impurity_parent : double
            The initial impurity of the parent node before the split

        impurity_left : double
            The impurity of the left child

        impurity_right : double
            The impurity of the right child

        Return
        ------
        double : improvement in impurity after the split occurs
        """
        return ((self.weighted_n_node_samples / self.weighted_n_samples) *
                (impurity_parent - (self.weighted_n_right /
                                    self.weighted_n_node_samples * impurity_right)
                                 - (self.weighted_n_left /
                                    self.weighted_n_node_samples * impurity_left)))

    cdef void set_sample_pointers(
        self,
        SIZE_t start,
        SIZE_t end
    ) noexcept nogil:
        """Abstract method which will set sample pointers in the criterion.

        The dataset array that we compute criteria on is assumed to consist of 'N'
        ordered samples or rows (i.e. sorted). Since we pass this by reference, we
        use sample pointers to move the start and end around to consider only a subset of data.
        This function should also update relevant statistics that the class uses to compute the final criterion.

        Parameters
        ----------
        start : SIZE_t
            The index of the first sample to be used on computation of criteria of the current node.
        end : SIZE_t
            The last sample used on this node
        """
        pass


cdef class Criterion(BaseCriterion):
    """Interface for impurity criteria.

    The supervised criterion computes the impurity of a node and the reduction of
    impurity of a split on that node using the distribution of labels in parent and
    children nodes. It also computes the output statistics such as the mean in regression
    and class probabilities in classification. Instances of this class are responsible
    for compute splits' impurity difference.

    Criterion is the base class for criteria used in supervised tree-based models
    with a homogeneous float64-dtyped y.
    """
    cdef int init(
        self,
        const DOUBLE_t[:, ::1] y,
        const DOUBLE_t[:] sample_weight,
        double weighted_n_samples,
        const SIZE_t[:] sample_indices
    ) except -1 nogil:
        """Placeholder for a method which will initialize the criterion.

        Returns -1 in case of failure to allocate memory (and raise MemoryError)
        or 0 otherwise.

        Parameters
        ----------
        y : ndarray, dtype=DOUBLE_t
            y is a buffer that can store values for n_outputs target variables
            stored as a Cython memoryview.
        sample_weight : ndarray, dtype=DOUBLE_t
            The weight of each sample stored as a Cython memoryview.
        weighted_n_samples : double
            The total weight of the samples being considered
        sample_indices : ndarray, dtype=SIZE_t
            A mask on the samples. Indices of the samples in X and y we want to use,
            where sample_indices[start:end] correspond to the samples in this node.
        """
        pass

    cdef void init_missing(self, SIZE_t n_missing) noexcept nogil:
        """Initialize sum_missing if there are missing values.

        This method assumes that caller placed the missing samples in
        self.sample_indices[-n_missing:]

        Parameters
        ----------
        n_missing: SIZE_t
            Number of missing values for specific feature.
        """
        pass

    cdef void clip_node_value(self, double* dest, double lower_bound, double upper_bound) noexcept nogil:
        pass

    cdef double middle_value(self) noexcept nogil:
        """Compute the middle value of a split for monotonicity constraints

        This method is implemented in ClassificationCriterion and RegressionCriterion.
        """
        pass

    cdef bint check_monotonicity(
        self,
        cnp.int8_t monotonic_cst,
        double lower_bound,
        double upper_bound,
    ) noexcept nogil:
        pass

    cdef inline bint _check_monotonicity(
        self,
        cnp.int8_t monotonic_cst,
        double lower_bound,
        double upper_bound,
        double value_left,
        double value_right,
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

    cdef void node_samples(
        self,
        vector[vector[DOUBLE_t]]& dest
    ) noexcept nogil:
        """Copy the samples of the current node into dest.

        Parameters
        ----------
        dest : reference vector[vector[DOUBLE_t]]
            The vector of vectors where the samples should be copied.
            This is passed by reference and modified in place.
        """
        cdef SIZE_t i, j, k

        # Resize the destination vector of vectors
        dest.resize(self.n_node_samples)

        # Loop over the samples
        for i in range(self.n_node_samples):
            # Get the index of the current sample
            j = self.sample_indices[self.start + i]

            # Get the sample values for each output
            for k in range(self.n_outputs):
                dest[i].push_back(self.y[j, k])


cdef inline void _move_sums_classification(
    ClassificationCriterion criterion,
    double[:, ::1] sum_1,
    double[:, ::1] sum_2,
    double* weighted_n_1,
    double* weighted_n_2,
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
    cdef SIZE_t k, c, n_bytes
    if criterion.n_missing != 0 and put_missing_in_1:
        for k in range(criterion.n_outputs):
            n_bytes = criterion.n_classes[k] * sizeof(double)
            memcpy(&sum_1[k, 0], &criterion.sum_missing[k, 0], n_bytes)

        for k in range(criterion.n_outputs):
            for c in range(criterion.n_classes[k]):
                sum_2[k, c] = criterion.sum_total[k, c] - criterion.sum_missing[k, c]

        weighted_n_1[0] = criterion.weighted_n_missing
        weighted_n_2[0] = criterion.weighted_n_node_samples - criterion.weighted_n_missing
    else:
        # Assigning sum_2 = sum_total for all outputs.
        for k in range(criterion.n_outputs):
            n_bytes = criterion.n_classes[k] * sizeof(double)
            memset(&sum_1[k, 0], 0, n_bytes)
            memcpy(&sum_2[k, 0], &criterion.sum_total[k, 0], n_bytes)

        weighted_n_1[0] = 0.0
        weighted_n_2[0] = criterion.weighted_n_node_samples


cdef class ClassificationCriterion(Criterion):
    """Abstract criterion for classification."""

    def __cinit__(self, SIZE_t n_outputs,
                  cnp.ndarray[SIZE_t, ndim=1] n_classes):
        """Initialize attributes for this criterion.

        Parameters
        ----------
        n_outputs : SIZE_t
            The number of targets, the dimensionality of the prediction
        n_classes : numpy.ndarray, dtype=SIZE_t
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

        cdef SIZE_t k = 0
        cdef SIZE_t max_n_classes = 0

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
        const DOUBLE_t[:, ::1] y,
        const DOUBLE_t[:] sample_weight,
        double weighted_n_samples,
        const SIZE_t[:] sample_indices
    ) except -1 nogil:
        """Initialize the criterion.

        Returns -1 in case of failure to allocate memory (and raise MemoryError)
        or 0 otherwise.

        Parameters
        ----------
        y : ndarray, dtype=DOUBLE_t
            The target stored as a buffer for memory efficiency.
        sample_weight : ndarray, dtype=DOUBLE_t
            The weight of each sample stored as a Cython memoryview.
        weighted_n_samples : double
            The total weight of all samples
        sample_indices : ndarray, dtype=SIZE_t
            A mask on the samples. Indices of the samples in X and y we want to use,
            where sample_indices[start:end] correspond to the samples in this node.
        """
        self.y = y
        self.sample_weight = sample_weight
        self.sample_indices = sample_indices
        self.weighted_n_samples = weighted_n_samples

        return 0

    cdef void set_sample_pointers(
        self,
        SIZE_t start,
        SIZE_t end
    ) noexcept nogil:
        """Set sample pointers in the criterion."""
        self.n_node_samples = end - start
        self.start = start
        self.end = end

        self.weighted_n_node_samples = 0.0

        cdef SIZE_t i
        cdef SIZE_t p
        cdef SIZE_t k
        cdef SIZE_t c
        cdef DOUBLE_t w = 1.0

        for k in range(self.n_outputs):
            memset(&self.sum_total[k, 0], 0, self.n_classes[k] * sizeof(double))

        for p in range(start, end):
            i = self.sample_indices[p]

            # w is originally set to be 1.0, meaning that if no sample weights
            # are given, the default weight of each sample is 1.0.
            if self.sample_weight is not None:
                w = self.sample_weight[i]

            # Count weighted class frequency for each target
            for k in range(self.n_outputs):
                c = <SIZE_t> self.y[i, k]
                self.sum_total[k, c] += w

            self.weighted_n_node_samples += w

        # Reset to pos=start
        self.reset()

    cdef void init_sum_missing(self):
        """Init sum_missing to hold sums for missing values."""
        self.sum_missing = np.zeros((self.n_outputs, self.max_n_classes), dtype=np.float64)

    cdef void init_missing(self, SIZE_t n_missing) noexcept nogil:
        """Initialize sum_missing if there are missing values.

        This method assumes that caller placed the missing samples in
        self.sample_indices[-n_missing:]
        """
        cdef SIZE_t i, p, k, c
        cdef DOUBLE_t w = 1.0

        self.n_missing = n_missing
        if n_missing == 0:
            return

        memset(&self.sum_missing[0, 0], 0, self.max_n_classes * self.n_outputs * sizeof(double))

        self.weighted_n_missing = 0.0

        # The missing samples are assumed to be in self.sample_indices[-n_missing:]
        for p in range(self.end - n_missing, self.end):
            i = self.sample_indices[p]
            if self.sample_weight is not None:
                w = self.sample_weight[i]

            for k in range(self.n_outputs):
                c = <SIZE_t> self.y[i, k]
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

    cdef int update(self, SIZE_t new_pos) except -1 nogil:
        """Updated statistics by moving sample_indices[pos:new_pos] to the left child.

        Returns -1 in case of failure to allocate memory (and raise MemoryError)
        or 0 otherwise.

        Parameters
        ----------
        new_pos : SIZE_t
            The new ending position for which to move sample_indices from the right
            child to the left child.
        """
        cdef SIZE_t pos = self.pos
        # The missing samples are assumed to be in
        # self.sample_indices[-self.n_missing:] that is
        # self.sample_indices[end_non_missing:self.end].
        cdef SIZE_t end_non_missing = self.end - self.n_missing

        cdef const SIZE_t[:] sample_indices = self.sample_indices
        cdef const DOUBLE_t[:] sample_weight = self.sample_weight

        cdef SIZE_t i
        cdef SIZE_t p
        cdef SIZE_t k
        cdef SIZE_t c
        cdef DOUBLE_t w = 1.0

        # Update statistics up to new_pos
        #
        # Given that
        #   sum_left[x] +  sum_right[x] = sum_total[x]
        # and that sum_total is known, we are going to update
        # sum_left from the direction that require the least amount
        # of computations, i.e. from pos to new_pos or from end to new_po.
        if (new_pos - pos) <= (end_non_missing - new_pos):
            for p in range(pos, new_pos):
                i = sample_indices[p]

                if sample_weight is not None:
                    w = sample_weight[i]

                for k in range(self.n_outputs):
                    self.sum_left[k, <SIZE_t> self.y[i, k]] += w

                self.weighted_n_left += w

        else:
            self.reverse_reset()

            for p in range(end_non_missing - 1, new_pos - 1, -1):
                i = sample_indices[p]

                if sample_weight is not None:
                    w = sample_weight[i]

                for k in range(self.n_outputs):
                    self.sum_left[k, <SIZE_t> self.y[i, k]] -= w

                self.weighted_n_left -= w

        # Update right part statistics
        self.weighted_n_right = self.weighted_n_node_samples - self.weighted_n_left
        for k in range(self.n_outputs):
            for c in range(self.n_classes[k]):
                self.sum_right[k, c] = self.sum_total[k, c] - self.sum_left[k, c]

        self.pos = new_pos
        return 0

    cdef double node_impurity(self) noexcept nogil:
        pass

    cdef void children_impurity(self, double* impurity_left,
                                double* impurity_right) noexcept nogil:
        pass

    cdef void node_value(self, double* dest) noexcept nogil:
        """Compute the node value of sample_indices[start:end] and save it into dest.

        Parameters
        ----------
        dest : double pointer
            The memory address which we will save the node value into.
        """
        cdef SIZE_t k

        for k in range(self.n_outputs):
            memcpy(dest, &self.sum_total[k, 0], self.n_classes[k] * sizeof(double))
            dest += self.max_n_classes

    cdef void clip_node_value(self, double * dest, double lower_bound, double upper_bound) noexcept nogil:
        """Clip the value in dest between lower_bound and upper_bound for monotonic constraints.

        Note that monotonicity constraints are only supported for:
        - single-output trees and
        - binary classifications.
        """
        if dest[0] < lower_bound:
            dest[0] = lower_bound
        elif dest[0] > upper_bound:
            dest[0] = upper_bound

        # Class proportions for binary classification must sum to 1.
        dest[1] = 1 - dest[0]

    cdef inline double middle_value(self) noexcept nogil:
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
        double lower_bound,
        double upper_bound,
    ) noexcept nogil:
        """Check monotonicity constraint is satisfied at the current classification split"""
        cdef:
            double value_left = self.sum_left[0][0] / self.weighted_n_left
            double value_right = self.sum_right[0][0] / self.weighted_n_right

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

    cdef double node_impurity(self) noexcept nogil:
        """Evaluate the impurity of the current node.

        Evaluate the cross-entropy criterion as impurity of the current node,
        i.e. the impurity of sample_indices[start:end]. The smaller the impurity the
        better.
        """
        cdef double entropy = 0.0
        cdef double count_k
        cdef SIZE_t k
        cdef SIZE_t c

        for k in range(self.n_outputs):
            for c in range(self.n_classes[k]):
                count_k = self.sum_total[k, c]
                if count_k > 0.0:
                    count_k /= self.weighted_n_node_samples
                    entropy -= count_k * log(count_k)

        return entropy / self.n_outputs

    cdef void children_impurity(self, double* impurity_left,
                                double* impurity_right) noexcept nogil:
        """Evaluate the impurity in children nodes.

        i.e. the impurity of the left child (sample_indices[start:pos]) and the
        impurity the right child (sample_indices[pos:end]).

        Parameters
        ----------
        impurity_left : double pointer
            The memory address to save the impurity of the left node
        impurity_right : double pointer
            The memory address to save the impurity of the right node
        """
        cdef double entropy_left = 0.0
        cdef double entropy_right = 0.0
        cdef double count_k
        cdef SIZE_t k
        cdef SIZE_t c

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

    cdef double node_impurity(self) noexcept nogil:
        """Evaluate the impurity of the current node.

        Evaluate the Gini criterion as impurity of the current node,
        i.e. the impurity of sample_indices[start:end]. The smaller the impurity the
        better.
        """
        cdef double gini = 0.0
        cdef double sq_count
        cdef double count_k
        cdef SIZE_t k
        cdef SIZE_t c

        for k in range(self.n_outputs):
            sq_count = 0.0

            for c in range(self.n_classes[k]):
                count_k = self.sum_total[k, c]
                sq_count += count_k * count_k

            gini += 1.0 - sq_count / (self.weighted_n_node_samples *
                                      self.weighted_n_node_samples)

        return gini / self.n_outputs

    cdef void children_impurity(self, double* impurity_left,
                                double* impurity_right) noexcept nogil:
        """Evaluate the impurity in children nodes.

        i.e. the impurity of the left child (sample_indices[start:pos]) and the
        impurity the right child (sample_indices[pos:end]) using the Gini index.

        Parameters
        ----------
        impurity_left : double pointer
            The memory address to save the impurity of the left node to
        impurity_right : double pointer
            The memory address to save the impurity of the right node to
        """
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
    double[::1] sum_1,
    double[::1] sum_2,
    double* weighted_n_1,
    double* weighted_n_2,
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
        SIZE_t i
        SIZE_t n_bytes = criterion.n_outputs * sizeof(double)
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

    def __cinit__(self, SIZE_t n_outputs, SIZE_t n_samples):
        """Initialize parameters for this criterion.

        Parameters
        ----------
        n_outputs : SIZE_t
            The number of targets to be predicted
        n_samples : SIZE_t
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
        const DOUBLE_t[:, ::1] y,
        const DOUBLE_t[:] sample_weight,
        double weighted_n_samples,
        const SIZE_t[:] sample_indices
    ) except -1 nogil:
        """Initialize the criterion."""
        # Initialize fields
        self.y = y
        self.sample_weight = sample_weight
        self.sample_indices = sample_indices
        self.weighted_n_samples = weighted_n_samples

        return 0

    cdef void set_sample_pointers(
        self,
        SIZE_t start,
        SIZE_t end
    ) noexcept nogil:
        """Set sample pointers in the criterion."""
        self.start = start
        self.end = end

        self.n_node_samples = end - start

        self.sq_sum_total = 0.0
        self.weighted_n_node_samples = 0.

        cdef SIZE_t i
        cdef SIZE_t p
        cdef SIZE_t k
        cdef DOUBLE_t y_ik
        cdef DOUBLE_t w_y_ik
        cdef DOUBLE_t w = 1.0

        memset(&self.sum_total[0], 0, self.n_outputs * sizeof(double))

        for p in range(start, end):
            i = self.sample_indices[p]

            if self.sample_weight is not None:
                w = self.sample_weight[i]

            for k in range(self.n_outputs):
                y_ik = self.y[i, k]
                w_y_ik = w * y_ik
                self.sum_total[k] += w_y_ik
                self.sq_sum_total += w_y_ik * y_ik

            self.weighted_n_node_samples += w

        # Reset to pos=start
        self.reset()

    cdef void init_sum_missing(self):
        """Init sum_missing to hold sums for missing values."""
        self.sum_missing = np.zeros(self.n_outputs, dtype=np.float64)

    cdef void init_missing(self, SIZE_t n_missing) noexcept nogil:
        """Initialize sum_missing if there are missing values.

        This method assumes that caller placed the missing samples in
        self.sample_indices[-n_missing:]
        """
        cdef SIZE_t i, p, k
        cdef DOUBLE_t y_ik
        cdef DOUBLE_t w_y_ik
        cdef DOUBLE_t w = 1.0

        self.n_missing = n_missing
        if n_missing == 0:
            return

        memset(&self.sum_missing[0], 0, self.n_outputs * sizeof(double))

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

    cdef int update(self, SIZE_t new_pos) except -1 nogil:
        """Updated statistics by moving sample_indices[pos:new_pos] to the left."""
        cdef const DOUBLE_t[:] sample_weight = self.sample_weight
        cdef const SIZE_t[:] sample_indices = self.sample_indices

        cdef SIZE_t pos = self.pos

        # The missing samples are assumed to be in
        # self.sample_indices[-self.n_missing:] that is
        # self.sample_indices[end_non_missing:self.end].
        cdef SIZE_t end_non_missing = self.end - self.n_missing
        cdef SIZE_t i
        cdef SIZE_t p
        cdef SIZE_t k
        cdef DOUBLE_t w = 1.0

        # Update statistics up to new_pos
        #
        # Given that
        #           sum_left[x] +  sum_right[x] = sum_total[x]
        # and that sum_total is known, we are going to update
        # sum_left from the direction that require the least amount
        # of computations, i.e. from pos to new_pos or from end to new_pos.
        if (new_pos - pos) <= (end_non_missing - new_pos):
            for p in range(pos, new_pos):
                i = sample_indices[p]

                if sample_weight is not None:
                    w = sample_weight[i]

                for k in range(self.n_outputs):
                    self.sum_left[k] += w * self.y[i, k]

                self.weighted_n_left += w
        else:
            self.reverse_reset()

            for p in range(end_non_missing - 1, new_pos - 1, -1):
                i = sample_indices[p]

                if sample_weight is not None:
                    w = sample_weight[i]

                for k in range(self.n_outputs):
                    self.sum_left[k] -= w * self.y[i, k]

                self.weighted_n_left -= w

        self.weighted_n_right = (self.weighted_n_node_samples -
                                 self.weighted_n_left)
        for k in range(self.n_outputs):
            self.sum_right[k] = self.sum_total[k] - self.sum_left[k]

        self.pos = new_pos
        return 0

    cdef double node_impurity(self) noexcept nogil:
        pass

    cdef void children_impurity(self, double* impurity_left,
                                double* impurity_right) noexcept nogil:
        pass

    cdef void node_value(self, double* dest) noexcept nogil:
        """Compute the node value of sample_indices[start:end] into dest."""
        cdef SIZE_t k

        for k in range(self.n_outputs):
            dest[k] = self.sum_total[k] / self.weighted_n_node_samples

    cdef inline void clip_node_value(self, double* dest, double lower_bound, double upper_bound) noexcept nogil:
        """Clip the value in dest between lower_bound and upper_bound for monotonic constraints."""
        if dest[0] < lower_bound:
            dest[0] = lower_bound
        elif dest[0] > upper_bound:
            dest[0] = upper_bound

    cdef double middle_value(self) noexcept nogil:
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
        double lower_bound,
        double upper_bound,
    ) noexcept nogil:
        """Check monotonicity constraint is satisfied at the current regression split"""
        cdef:
            double value_left = self.sum_left[0] / self.weighted_n_left
            double value_right = self.sum_right[0] / self.weighted_n_right

        return self._check_monotonicity(monotonic_cst, lower_bound, upper_bound, value_left, value_right)

cdef class MSE(RegressionCriterion):
    """Mean squared error impurity criterion.
        MSE = var_left + var_right
    """

    cdef double node_impurity(self) noexcept nogil:
        """Evaluate the impurity of the current node.

        Evaluate the MSE criterion as impurity of the current node,
        i.e. the impurity of sample_indices[start:end]. The smaller the impurity the
        better.
        """
        cdef double impurity
        cdef SIZE_t k

        impurity = self.sq_sum_total / self.weighted_n_node_samples
        for k in range(self.n_outputs):
            impurity -= (self.sum_total[k] / self.weighted_n_node_samples)**2.0

        return impurity / self.n_outputs

    cdef double proxy_impurity_improvement(self) noexcept nogil:
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
        cdef SIZE_t k
        cdef double proxy_impurity_left = 0.0
        cdef double proxy_impurity_right = 0.0

        for k in range(self.n_outputs):
            proxy_impurity_left += self.sum_left[k] * self.sum_left[k]
            proxy_impurity_right += self.sum_right[k] * self.sum_right[k]

        return (proxy_impurity_left / self.weighted_n_left +
                proxy_impurity_right / self.weighted_n_right)

    cdef void children_impurity(self, double* impurity_left,
                                double* impurity_right) noexcept nogil:
        """Evaluate the impurity in children nodes.

        i.e. the impurity of the left child (sample_indices[start:pos]) and the
        impurity the right child (sample_indices[pos:end]).
        """
        cdef const DOUBLE_t[:] sample_weight = self.sample_weight
        cdef const SIZE_t[:] sample_indices = self.sample_indices
        cdef SIZE_t pos = self.pos
        cdef SIZE_t start = self.start

        cdef DOUBLE_t y_ik

        cdef double sq_sum_left = 0.0
        cdef double sq_sum_right

        cdef SIZE_t i
        cdef SIZE_t p
        cdef SIZE_t k
        cdef DOUBLE_t w = 1.0

        for p in range(start, pos):
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


cdef class MAE(RegressionCriterion):
    r"""Mean absolute error impurity criterion.

       MAE = (1 / n)*(\sum_i |y_i - f_i|), where y_i is the true
       value and f_i is the predicted value."""

    cdef cnp.ndarray left_child
    cdef cnp.ndarray right_child
    cdef void** left_child_ptr
    cdef void** right_child_ptr
    cdef DOUBLE_t[::1] node_medians

    def __cinit__(self, SIZE_t n_outputs, SIZE_t n_samples):
        """Initialize parameters for this criterion.

        Parameters
        ----------
        n_outputs : SIZE_t
            The number of targets to be predicted

        n_samples : SIZE_t
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

        self.node_medians = np.zeros(n_outputs, dtype=np.float64)

        self.left_child = np.empty(n_outputs, dtype='object')
        self.right_child = np.empty(n_outputs, dtype='object')
        # initialize WeightedMedianCalculators
        for k in range(n_outputs):
            self.left_child[k] = WeightedMedianCalculator(n_samples)
            self.right_child[k] = WeightedMedianCalculator(n_samples)

        self.left_child_ptr = <void**> cnp.PyArray_DATA(self.left_child)
        self.right_child_ptr = <void**> cnp.PyArray_DATA(self.right_child)

    cdef int init(
        self,
        const DOUBLE_t[:, ::1] y,
        const DOUBLE_t[:] sample_weight,
        double weighted_n_samples,
        const SIZE_t[:] sample_indices
    ) except -1 nogil:
        """Initialize the criterion."""
        # Initialize fields
        self.y = y
        self.sample_weight = sample_weight
        self.sample_indices = sample_indices
        self.weighted_n_samples = weighted_n_samples

        return 0

    cdef void set_sample_pointers(
        self,
        SIZE_t start,
        SIZE_t end
    ) noexcept nogil:
        """Set sample pointers in the criterion."""
        cdef SIZE_t i, p, k
        cdef DOUBLE_t w = 1.0

        self.start = start
        self.end = end

        self.n_node_samples = end - start
        self.weighted_n_node_samples = 0.

        cdef void** left_child = self.left_child_ptr
        cdef void** right_child = self.right_child_ptr

        for k in range(self.n_outputs):
            (<WeightedMedianCalculator> left_child[k]).reset()
            (<WeightedMedianCalculator> right_child[k]).reset()

        for p in range(start, end):
            i = self.sample_indices[p]

            if self.sample_weight is not None:
                w = self.sample_weight[i]

            for k in range(self.n_outputs):
                # push method ends up calling safe_realloc, hence `except -1`
                # push all values to the right side,
                # since pos = start initially anyway
                (<WeightedMedianCalculator> right_child[k]).push(self.y[i, k], w)

            self.weighted_n_node_samples += w
        # calculate the node medians
        for k in range(self.n_outputs):
            self.node_medians[k] = (<WeightedMedianCalculator> right_child[k]).get_median()

        # Reset to pos=start
        self.reset()

    cdef void init_missing(self, SIZE_t n_missing) noexcept nogil:
        """Raise error if n_missing != 0."""
        if n_missing == 0:
            return
        with gil:
            raise ValueError("missing values is not supported for MAE.")

    cdef int reset(self) except -1 nogil:
        """Reset the criterion at pos=start.

        Returns -1 in case of failure to allocate memory (and raise MemoryError)
        or 0 otherwise.
        """
        cdef SIZE_t i, k
        cdef DOUBLE_t value
        cdef DOUBLE_t weight

        cdef void** left_child = self.left_child_ptr
        cdef void** right_child = self.right_child_ptr

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
                # push method ends up calling safe_realloc, hence `except -1`
                (<WeightedMedianCalculator> right_child[k]).push(value,
                                                                 weight)
        return 0

    cdef int reverse_reset(self) except -1 nogil:
        """Reset the criterion at pos=end.

        Returns -1 in case of failure to allocate memory (and raise MemoryError)
        or 0 otherwise.
        """
        self.weighted_n_right = 0.0
        self.weighted_n_left = self.weighted_n_node_samples
        self.pos = self.end

        cdef DOUBLE_t value
        cdef DOUBLE_t weight
        cdef void** left_child = self.left_child_ptr
        cdef void** right_child = self.right_child_ptr

        # reverse reset the WeightedMedianCalculators, right should have no
        # elements and left should have all elements.
        for k in range(self.n_outputs):
            # if right has no elements, it's already reset
            for i in range((<WeightedMedianCalculator> right_child[k]).size()):
                # remove everything from right and put it into left
                (<WeightedMedianCalculator> right_child[k]).pop(&value,
                                                                &weight)
                # push method ends up calling safe_realloc, hence `except -1`
                (<WeightedMedianCalculator> left_child[k]).push(value,
                                                                weight)
        return 0

    cdef int update(self, SIZE_t new_pos) except -1 nogil:
        """Updated statistics by moving sample_indices[pos:new_pos] to the left.

        Returns -1 in case of failure to allocate memory (and raise MemoryError)
        or 0 otherwise.
        """
        cdef const DOUBLE_t[:] sample_weight = self.sample_weight
        cdef const SIZE_t[:] sample_indices = self.sample_indices

        cdef void** left_child = self.left_child_ptr
        cdef void** right_child = self.right_child_ptr

        cdef SIZE_t pos = self.pos
        cdef SIZE_t end = self.end
        cdef SIZE_t i, p, k
        cdef DOUBLE_t w = 1.0

        # Update statistics up to new_pos
        #
        # We are going to update right_child and left_child
        # from the direction that require the least amount of
        # computations, i.e. from pos to new_pos or from end to new_pos.
        if (new_pos - pos) <= (end - new_pos):
            for p in range(pos, new_pos):
                i = sample_indices[p]

                if sample_weight is not None:
                    w = sample_weight[i]

                for k in range(self.n_outputs):
                    # remove y_ik and its weight w from right and add to left
                    (<WeightedMedianCalculator> right_child[k]).remove(self.y[i, k], w)
                    # push method ends up calling safe_realloc, hence except -1
                    (<WeightedMedianCalculator> left_child[k]).push(self.y[i, k], w)

                self.weighted_n_left += w
        else:
            self.reverse_reset()

            for p in range(end - 1, new_pos - 1, -1):
                i = sample_indices[p]

                if sample_weight is not None:
                    w = sample_weight[i]

                for k in range(self.n_outputs):
                    # remove y_ik and its weight w from left and add to right
                    (<WeightedMedianCalculator> left_child[k]).remove(self.y[i, k], w)
                    (<WeightedMedianCalculator> right_child[k]).push(self.y[i, k], w)

                self.weighted_n_left -= w

        self.weighted_n_right = (self.weighted_n_node_samples -
                                 self.weighted_n_left)
        self.pos = new_pos
        return 0

    cdef void node_value(self, double* dest) noexcept nogil:
        """Computes the node value of sample_indices[start:end] into dest."""
        cdef SIZE_t k
        for k in range(self.n_outputs):
            dest[k] = <double> self.node_medians[k]

    cdef inline double middle_value(self) noexcept nogil:
        """Compute the middle value of a split for monotonicity constraints as the simple average
        of the left and right children values.

        Monotonicity constraints are only supported for single-output trees we can safely assume
        n_outputs == 1.
        """
        return (
                (<WeightedMedianCalculator> self.left_child_ptr[0]).get_median() +
                (<WeightedMedianCalculator> self.right_child_ptr[0]).get_median()
        ) / 2

    cdef inline bint check_monotonicity(
        self,
        cnp.int8_t monotonic_cst,
        double lower_bound,
        double upper_bound,
    ) noexcept nogil:
        """Check monotonicity constraint is satisfied at the current regression split"""
        cdef:
            double value_left = (<WeightedMedianCalculator> self.left_child_ptr[0]).get_median()
            double value_right = (<WeightedMedianCalculator> self.right_child_ptr[0]).get_median()

        return self._check_monotonicity(monotonic_cst, lower_bound, upper_bound, value_left, value_right)

    cdef double node_impurity(self) noexcept nogil:
        """Evaluate the impurity of the current node.

        Evaluate the MAE criterion as impurity of the current node,
        i.e. the impurity of sample_indices[start:end]. The smaller the impurity the
        better.
        """
        cdef const DOUBLE_t[:] sample_weight = self.sample_weight
        cdef const SIZE_t[:] sample_indices = self.sample_indices
        cdef SIZE_t i, p, k
        cdef DOUBLE_t w = 1.0
        cdef DOUBLE_t impurity = 0.0

        for k in range(self.n_outputs):
            for p in range(self.start, self.end):
                i = sample_indices[p]

                if sample_weight is not None:
                    w = sample_weight[i]

                impurity += fabs(self.y[i, k] - self.node_medians[k]) * w

        return impurity / (self.weighted_n_node_samples * self.n_outputs)

    cdef void children_impurity(self, double* p_impurity_left,
                                double* p_impurity_right) noexcept nogil:
        """Evaluate the impurity in children nodes.

        i.e. the impurity of the left child (sample_indices[start:pos]) and the
        impurity the right child (sample_indices[pos:end]).
        """
        cdef const DOUBLE_t[:] sample_weight = self.sample_weight
        cdef const SIZE_t[:] sample_indices = self.sample_indices

        cdef SIZE_t start = self.start
        cdef SIZE_t pos = self.pos
        cdef SIZE_t end = self.end

        cdef SIZE_t i, p, k
        cdef DOUBLE_t median
        cdef DOUBLE_t w = 1.0
        cdef DOUBLE_t impurity_left = 0.0
        cdef DOUBLE_t impurity_right = 0.0

        cdef void** left_child = self.left_child_ptr
        cdef void** right_child = self.right_child_ptr

        for k in range(self.n_outputs):
            median = (<WeightedMedianCalculator> left_child[k]).get_median()
            for p in range(start, pos):
                i = sample_indices[p]

                if sample_weight is not None:
                    w = sample_weight[i]

                impurity_left += fabs(self.y[i, k] - median) * w
        p_impurity_left[0] = impurity_left / (self.weighted_n_left *
                                              self.n_outputs)

        for k in range(self.n_outputs):
            median = (<WeightedMedianCalculator> right_child[k]).get_median()
            for p in range(pos, end):
                i = sample_indices[p]

                if sample_weight is not None:
                    w = sample_weight[i]

                impurity_right += fabs(self.y[i, k] - median) * w
        p_impurity_right[0] = impurity_right / (self.weighted_n_right *
                                                self.n_outputs)


cdef class FriedmanMSE(MSE):
    """Mean squared error impurity criterion with improvement score by Friedman.

    Uses the formula (35) in Friedman's original Gradient Boosting paper:

        diff = mean_left - mean_right
        improvement = n_left * n_right * diff^2 / (n_left + n_right)
    """

    cdef double proxy_impurity_improvement(self) noexcept nogil:
        """Compute a proxy of the impurity reduction.

        This method is used to speed up the search for the best split.
        It is a proxy quantity such that the split that maximizes this value
        also maximizes the impurity improvement. It neglects all constant terms
        of the impurity decrease for a given split.

        The absolute impurity improvement is only computed by the
        impurity_improvement method once the best split has been found.
        """
        cdef double total_sum_left = 0.0
        cdef double total_sum_right = 0.0

        cdef SIZE_t k
        cdef double diff = 0.0

        for k in range(self.n_outputs):
            total_sum_left += self.sum_left[k]
            total_sum_right += self.sum_right[k]

        diff = (self.weighted_n_right * total_sum_left -
                self.weighted_n_left * total_sum_right)

        return diff * diff / (self.weighted_n_left * self.weighted_n_right)

    cdef double impurity_improvement(self, double impurity_parent, double
                                     impurity_left, double impurity_right) noexcept nogil:
        # Note: none of the arguments are used here
        cdef double total_sum_left = 0.0
        cdef double total_sum_right = 0.0

        cdef SIZE_t k
        cdef double diff = 0.0

        for k in range(self.n_outputs):
            total_sum_left += self.sum_left[k]
            total_sum_right += self.sum_right[k]

        diff = (self.weighted_n_right * total_sum_left -
                self.weighted_n_left * total_sum_right) / self.n_outputs

        return (diff * diff / (self.weighted_n_left * self.weighted_n_right *
                               self.weighted_n_node_samples))


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

    cdef double node_impurity(self) noexcept nogil:
        """Evaluate the impurity of the current node.

        Evaluate the Poisson criterion as impurity of the current node,
        i.e. the impurity of sample_indices[start:end]. The smaller the impurity the
        better.
        """
        return self.poisson_loss(self.start, self.end, self.sum_total,
                                 self.weighted_n_node_samples)

    cdef double proxy_impurity_improvement(self) noexcept nogil:
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
        cdef SIZE_t k
        cdef double proxy_impurity_left = 0.0
        cdef double proxy_impurity_right = 0.0
        cdef double y_mean_left = 0.
        cdef double y_mean_right = 0.

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

    cdef void children_impurity(self, double* impurity_left,
                                double* impurity_right) noexcept nogil:
        """Evaluate the impurity in children nodes.

        i.e. the impurity of the left child (sample_indices[start:pos]) and the
        impurity of the right child (sample_indices[pos:end]) for Poisson.
        """
        cdef SIZE_t start = self.start
        cdef SIZE_t pos = self.pos
        cdef SIZE_t end = self.end

        impurity_left[0] = self.poisson_loss(start, pos, self.sum_left,
                                             self.weighted_n_left)

        impurity_right[0] = self.poisson_loss(pos, end, self.sum_right,
                                              self.weighted_n_right)

    cdef inline DOUBLE_t poisson_loss(self,
                                      SIZE_t start,
                                      SIZE_t end,
                                      const double[::1] y_sum,
                                      DOUBLE_t weight_sum) noexcept nogil:
        """Helper function to compute Poisson loss (~deviance) of a given node.
        """
        cdef const DOUBLE_t[:, ::1] y = self.y
        cdef const DOUBLE_t[:] sample_weight = self.sample_weight
        cdef const SIZE_t[:] sample_indices = self.sample_indices

        cdef DOUBLE_t y_mean = 0.
        cdef DOUBLE_t poisson_loss = 0.
        cdef DOUBLE_t w = 1.0
        cdef SIZE_t i, k, p
        cdef SIZE_t n_outputs = self.n_outputs

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
