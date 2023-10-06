# Authors: Gilles Louppe <g.louppe@gmail.com>
#          Peter Prettenhofer <peter.prettenhofer@gmail.com>
#          Brian Holt <bdholt1@gmail.com>
#          Joel Nothman <joel.nothman@gmail.com>
#          Arnaud Joly <arnaud.v.joly@gmail.com>
#          Jacob Schreiber <jmschreiber91@gmail.com>
#
# License: BSD 3 clause

# See _criterion.pyx for implementation details.
cimport numpy as cnp

from ..utils._typedefs cimport float64_t, intp_t


cdef class Criterion:
    # The criterion computes the impurity of a node and the reduction of
    # impurity of a split on that node. It also computes the output statistics
    # such as the mean in regression and class probabilities in classification.

    # Internal structures
    cdef const float64_t[:, ::1] y         # Values of y
    cdef const float64_t[:] sample_weight  # Sample weights

    cdef const intp_t[:] sample_indices    # Sample indices in X, y
    cdef intp_t start                      # samples[start:pos] are the samples in the left node
    cdef intp_t pos                        # samples[pos:end] are the samples in the right node
    cdef intp_t end
    cdef intp_t n_missing                  # Number of missing values for the feature being evaluated
    cdef bint missing_go_to_left           # Whether missing values go to the left node

    cdef intp_t n_outputs                  # Number of outputs
    cdef intp_t n_samples                  # Number of samples
    cdef intp_t n_node_samples             # Number of samples in the node (end-start)
    cdef double weighted_n_samples         # Weighted number of samples (in total)
    cdef double weighted_n_node_samples    # Weighted number of samples in the node
    cdef double weighted_n_left            # Weighted number of samples in the left node
    cdef double weighted_n_right           # Weighted number of samples in the right node
    cdef double weighted_n_missing         # Weighted number of samples that are missing

    # The criterion object is maintained such that left and right collected
    # statistics correspond to samples[start:pos] and samples[pos:end].

    # Methods
    cdef int init(
        self,
        const float64_t[:, ::1] y,
        const float64_t[:] sample_weight,
        double weighted_n_samples,
        const intp_t[:] sample_indices,
        intp_t start,
        intp_t end
    ) except -1 nogil
    cdef void init_sum_missing(self)
    cdef void init_missing(self, intp_t n_missing) noexcept nogil
    cdef int reset(self) except -1 nogil
    cdef int reverse_reset(self) except -1 nogil
    cdef int update(self, intp_t new_pos) except -1 nogil
    cdef double node_impurity(self) noexcept nogil
    cdef void children_impurity(
        self,
        double* impurity_left,
        double* impurity_right
    ) noexcept nogil
    cdef void node_value(
        self,
        double* dest
    ) noexcept nogil
    cdef void clip_node_value(
        self,
        double* dest,
        double lower_bound,
        double upper_bound
    ) noexcept nogil
    cdef double middle_value(self) noexcept nogil
    cdef double impurity_improvement(
        self,
        double impurity_parent,
        double impurity_left,
        double impurity_right
    ) noexcept nogil
    cdef double proxy_impurity_improvement(self) noexcept nogil
    cdef bint check_monotonicity(
            self,
            cnp.int8_t monotonic_cst,
            double lower_bound,
            double upper_bound,
    ) noexcept nogil
    cdef inline bint _check_monotonicity(
            self,
            cnp.int8_t monotonic_cst,
            double lower_bound,
            double upper_bound,
            double sum_left,
            double sum_right,
    ) noexcept nogil

cdef class ClassificationCriterion(Criterion):
    """Abstract criterion for classification."""

    cdef intp_t[::1] n_classes
    cdef intp_t max_n_classes

    cdef double[:, ::1] sum_total    # The sum of the weighted count of each label.
    cdef double[:, ::1] sum_left     # Same as above, but for the left side of the split
    cdef double[:, ::1] sum_right    # Same as above, but for the right side of the split
    cdef double[:, ::1] sum_missing  # Same as above, but for missing values in X

cdef class RegressionCriterion(Criterion):
    """Abstract regression criterion."""

    cdef double sq_sum_total

    cdef double[::1] sum_total    # The sum of w*y.
    cdef double[::1] sum_left     # Same as above, but for the left side of the split
    cdef double[::1] sum_right    # Same as above, but for the right side of the split
    cdef double[::1] sum_missing  # Same as above, but for missing values in X
