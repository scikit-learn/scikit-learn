# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# See _splitter.pyx for details.

from ..utils._bitset cimport BITSET_INNER_DTYPE_C
from ..utils._typedefs cimport (
    float32_t, float64_t, int8_t, int32_t, intp_t, uint8_t, uint32_t, uint64_t
)
from ._utils cimport ParentInfo, SplitRecord, SplitValue
from ._criterion cimport Criterion


cdef class Splitter:
    # The splitter searches in the input space for a feature and a threshold
    # to split the samples samples[start:end].
    #
    # The impurity computations are delegated to a criterion object.

    # Internal structures
    cdef public Criterion criterion      # Impurity criterion
    cdef public intp_t max_features      # Number of features to test
    cdef public intp_t min_samples_leaf  # Min samples in a leaf
    cdef public float64_t min_weight_leaf   # Minimum weight in a leaf

    cdef object random_state             # Random state
    cdef uint32_t rand_r_state           # sklearn_rand_r random number state

    cdef intp_t[::1] samples             # Sample indices in X, y
    cdef intp_t n_samples                # X.shape[0]
    cdef float64_t weighted_n_samples       # Weighted number of samples
    cdef intp_t[::1] features            # Feature indices in X
    cdef intp_t[::1] constant_features   # Constant features indices
    cdef intp_t n_features               # X.shape[1]
    cdef float32_t[::1] feature_values   # temp. array holding feature values

    cdef intp_t start                    # Start position for the current node
    cdef intp_t end                      # End position for the current node

    cdef const float64_t[:, ::1] y
    # Monotonicity constraints for each feature.
    # The encoding is as follows:
    #   -1: monotonic decrease
    #    0: no constraint
    #   +1: monotonic increase
    cdef const int8_t[:] monotonic_cst
    cdef bint with_monotonic_cst
    cdef const float64_t[:] sample_weight

    # Whether or not to sort categories by probabilities to split categorical
    # features using the Breiman shortcut;
    # Used to accelerate finding a greedy categorical split for binary classification
    # with Gini Impurity
    # or univariate regression with MSE.
    # XXX: This could be a compile-time check with C++ 17's SFINAE
    cdef bint breiman_shortcut

    # We know the number of categories within our dataset across each feature.
    # If a feature index has -1, then it is not categorical
    cdef const int32_t[:] n_categories
    cdef int32_t max_n_categories
    cdef BITSET_INNER_DTYPE_C[:] cat_split

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
    cdef int init(
        self,
        object X,
        const float64_t[:, ::1] y,
        const float64_t[:] sample_weight,
        const uint8_t[::1] missing_values_in_feature_mask,
        const int32_t[::1] n_categories,
    ) except -1

    cdef int node_reset(
        self,
        intp_t start,
        intp_t end,
        float64_t* weighted_n_node_samples
    ) except -1 nogil

    cdef int node_split(
        self,
        ParentInfo* parent,
        SplitRecord* split,
    ) except -1 nogil

    cdef void node_value(self, float64_t* dest) noexcept nogil

    cdef void clip_node_value(self, float64_t* dest, float64_t lower_bound, float64_t upper_bound) noexcept nogil

    cdef float64_t node_impurity(self) noexcept nogil
