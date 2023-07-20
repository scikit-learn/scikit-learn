# Authors: Gilles Louppe <g.louppe@gmail.com>
#          Peter Prettenhofer <peter.prettenhofer@gmail.com>
#          Brian Holt <bdholt1@gmail.com>
#          Joel Nothman <joel.nothman@gmail.com>
#          Arnaud Joly <arnaud.v.joly@gmail.com>
#          Jacob Schreiber <jmschreiber91@gmail.com>
#          Adam Li <adam2392@gmail.com>
#          Jong Shin <jshinm@gmail.com>
#
# License: BSD 3 clause

# See _splitter.pyx for details.
cimport numpy as cnp

from libcpp.vector cimport vector

from ._criterion cimport BaseCriterion, Criterion
from ._utils cimport DOUBLE_t  # Type of y, sample_weight
from ._utils cimport DTYPE_t  # Type of X
from ._utils cimport INT32_t  # Signed 32 bit integer
from ._utils cimport SIZE_t  # Type for indices and counters
from ._utils cimport UINT32_t  # Unsigned 32 bit integer
from ._utils cimport UINT64_t  # Unsigned 64 bit integer

from ._utils cimport SplitValue, SplitRecord, Node


cdef class BaseSplitter:
    """Abstract interface for splitter."""

    # The splitter searches in the input space for a feature and a threshold
    # to split the samples samples[start:end].
    #
    # The impurity computations are delegated to a criterion object.

    # Internal structures
    cdef public SIZE_t max_features      # Number of features to test
    cdef public SIZE_t min_samples_leaf  # Min samples in a leaf
    cdef public double min_weight_leaf   # Minimum weight in a leaf

    cdef object random_state             # Random state
    cdef UINT32_t rand_r_state           # sklearn_rand_r random number state

    cdef SIZE_t[::1] samples             # Sample indices in X, y
    cdef SIZE_t n_samples                # X.shape[0]
    cdef double weighted_n_samples       # Weighted number of samples
    cdef SIZE_t[::1] features            # Feature indices in X
    cdef SIZE_t[::1] constant_features   # Constant features indices
    cdef SIZE_t n_features               # X.shape[1]
    cdef DTYPE_t[::1] feature_values     # temp. array holding feature values

    cdef SIZE_t start                    # Start position for the current node
    cdef SIZE_t end                      # End position for the current node

    cdef const DOUBLE_t[:] sample_weight

    # The samples vector `samples` is maintained by the Splitter object such
    # that the samples contained in a node are contiguous. With this setting,
    # `node_split` reorganizes the node samples `samples[start:end]` in two
    # subsets `samples[start:pos]` and `samples[pos:end]`.

    # The 1-d  `features` array of size n_features contains the features
    # indices and allows fast sampling without replacement of features.

    # The 1-d `constant_features` array of size n_features holds in
    # `constant_features[:n_constant_features]` the feature ids with
    # constant values for all the samples that reached a specific node.
    # The value `n_constant_features` is given by the parent node to its
    # child nodes.  The content of the range `[n_constant_features:]` is left
    # undefined, but preallocated for performance reasons
    # This allows optimization with depth-based tree building.

    # Methods
    cdef int node_reset(
        self,
        SIZE_t start,
        SIZE_t end,
        double* weighted_n_node_samples
    ) except -1 nogil
    cdef int node_split(
        self,
        double impurity,   # Impurity of the node
        SplitRecord* split,
        SIZE_t* n_constant_features,
        double lower_bound,
        double upper_bound,
    ) except -1 nogil
    cdef void node_value(self, double* dest) noexcept nogil
    cdef double node_impurity(self) noexcept nogil
    cdef int pointer_size(self) noexcept nogil

cdef class Splitter(BaseSplitter):
    cdef public Criterion criterion     # Impurity criterion
    cdef const DOUBLE_t[:, ::1] y

    cdef INT32_t[:] n_categories        # (n_features,) array giving number of
    #                                   # categories (<0 for non-categorical)
    cdef UINT64_t[:] cat_cache          # Cache buffer for fast categorical split evaluation
    cdef bint breiman_shortcut          # Whether decision trees are allowed to use the
    #                                   # Breiman shortcut for categorical features
    #                                   # during binary classification.

    # Monotonicity constraints for each feature.
    # The encoding is as follows:
    #   -1: monotonic decrease
    #    0: no constraint
    #   +1: monotonic increase
    cdef const cnp.int8_t[:] monotonic_cst
    cdef bint with_monotonic_cst

    cdef int init(
        self,
        object X,
        const DOUBLE_t[:, ::1] y,
        const DOUBLE_t[:] sample_weight,
        const unsigned char[::1] missing_values_in_feature_mask,
        const INT32_t[:] n_categories,
    ) except -1

    cdef void node_samples(self, vector[vector[DOUBLE_t]]& dest) noexcept nogil

    # Methods that allow modifications to stopping conditions
    cdef bint check_presplit_conditions(
        self,
        SplitRecord current_split,
        SIZE_t n_missing,
        bint missing_go_to_left,
    ) noexcept nogil

    cdef bint check_postsplit_conditions(
        self
    ) noexcept nogil

    cdef void clip_node_value(
        self,
        double* dest,
        double lower_bound,
        double upper_bound
    ) noexcept nogil

    cdef void _breiman_sort_categories(
        self,
        SIZE_t start,
        SIZE_t end,
        INT32_t ncat,
        SIZE_t ncat_present,
        const INT32_t *cat_offset,
        SIZE_t *sorted_cat
    ) noexcept nogil
