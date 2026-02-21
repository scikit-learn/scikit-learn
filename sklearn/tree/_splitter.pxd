# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# See _splitter.pyx for details.

from sklearn.utils._typedefs cimport (
    float32_t, float64_t, int8_t, int32_t, intp_t, uint8_t, uint32_t, uint64_t
)
from sklearn.utils._bitset cimport BITSET_INNER_DTYPE_C
from sklearn.tree._criterion cimport Criterion
from sklearn.tree._tree cimport ParentInfo


ctypedef union SplitValue:
    # Union type to generalize the concept of a threshold to categorical
    # features. The floating point view, i.e. ``SplitValue.split_value.threshold`` is used
    # for numerical features, where feature values less than or equal to the
    # threshold go left, and values greater than the threshold go right.
    #
    # For categorical features, the BITSET_INNER_DTYPE_C view (`SplitValue.cat_split``) is
    # used. It works in one of two ways, indicated by the value of its least
    # significant bit (LSB). If the LSB is 0, then cat_split acts as a bitfield
    # for up to 64 categories, sending samples left if the bit corresponding to
    # their category is 1 or right if it is 0. If the LSB is 1, then the most
    # significant 32 bits of cat_split make a random seed. To evaluate a
    # sample, use the random seed to flip a coin (category_value + 1) times and
    # send it left if the last flip gives 1; otherwise right. This second
    # method allows up to 2**31 category values, but can only be used for
    # RandomSplitter.
    float64_t threshold
    # BITSET_DTYPE_C cat_split
    uint64_t[8] categorical_bitset

cdef struct SplitRecord:
    # Data to track sample split
    intp_t feature         # Which feature to split on.
    intp_t pos             # Split samples array at the given position,
    #                      # i.e. count of samples below threshold for feature.
    #                      # pos is >= end if the node is a leaf.

    SplitValue split_value    # Generalized threshold for categorical and
    #                         # non-categorical features to split samples.

    float64_t improvement     # Impurity improvement given parent node.
    float64_t impurity_left   # Impurity of the left split.
    float64_t impurity_right  # Impurity of the right split.
    float64_t lower_bound     # Lower bound on value of both children for monotonicity
    float64_t upper_bound     # Upper bound on value of both children for monotonicity
    uint8_t missing_go_to_left  # Controls if missing values go to the left node.
    intp_t n_missing            # Number of missing values for the feature being split on

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
        const int32_t[::1] n_categories
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
