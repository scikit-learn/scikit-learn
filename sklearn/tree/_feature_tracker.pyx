# Feature Tracker used with the splitter to sample features
cimport cython

from libc.string cimport memcpy
from ._utils cimport rand_int

import numpy as np

@cython.final
cdef class FeatureTracker:
    """Feature Sampler using Fisher-Yates-based algorithm.

    Sample up to max_features without replacement using a Fisher-Yates-based algorithm
    (using the local variables `f_i` and `f_j` to compute a permutation of the
    `features` array).

    Skip the CPU intensive evaluation of the impurity criterion for features that were
    already detected as constant (hence not suitable for good splitting) by ancestor
    nodes and save the information on newly discovered constant features to spare
    computation on descendant nodes.
    """

    def __init__(self, SIZE_t n_features, SIZE_t max_features):
        self.features = np.arange(n_features, dtype=np.intp)
        self.constant_features = np.empty(n_features, dtype=np.intp)
        self.max_features = max_features

    cdef inline void reset(self, SIZE_t n_constant_features) nogil:
        """Reset for feature sampler."""
        self.f_i = self.features.shape[0]
        self.n_visited_features = 0
        self.n_found_constants = 0
        self.n_drawn_constants = 0
        self.n_known_constants = n_constant_features
        self.n_total_constants = n_constant_features

    cdef inline FeatureSample sample(self, UINT32_t* random_state) nogil:
        """Sample for new feature to check.

        This function can return a struct with a feature index, it's value, and status.
        """
        cdef:
            SIZE_t f_j
            bint should_continue
            SIZE_t[::1] features = self.features
            FeatureSample output

        should_continue = (
            self.f_i > self.n_total_constants and  # Stop early if remaining features are constant
            (self.n_visited_features < self.max_features or
            # At least one drawn features must be non constant
            self.n_visited_features <= self.n_found_constants + self.n_drawn_constants))

        if not should_continue:
            output.status = FeatureStatus.STOP
            return output

        self.n_visited_features += 1

        # Loop invariant: elements of features in
        # - [:n_drawn_constant[ holds drawn and known constant features;
        # - [n_drawn_constant:n_known_constant[ holds known constant
        #   features that haven't been drawn yet;
        # - [n_known_constant:n_total_constant[ holds newly found constant
        #   features;
        # - [n_total_constant:f_i[ holds features that haven't been drawn
        #   yet and aren't constant apriori.
        # - [f_i:n_features[ holds features that have been drawn
        #   and aren't constant.

        # Draw a feature at random
        f_j = rand_int(self.n_drawn_constants, self.f_i - self.n_found_constants,
                       random_state)

        if f_j < self.n_known_constants:
            # f_j in the interval [n_drawn_constants, n_known_constants[
            features[f_j], features[self.n_drawn_constants] = (
                features[self.n_drawn_constants], features[f_j]
            )

            self.n_drawn_constants += 1
            output.status = FeatureStatus.CONTINUE
            return output

        # f_j in the interval [n_known_constants, f_i - n_found_constants[
        f_j += self.n_found_constants
        # f_j in the interval [n_total_constants, f_i[
        output.f_j = f_j
        output.feature = features[f_j]
        output.status = FeatureStatus.EVALUTE
        return output

    cdef inline void update_found_constant(self, SIZE_t f_j) nogil:
        """Mark f_j as constant."""
        cdef SIZE_t[::1] features = self.features
        features[f_j], features[self.n_total_constants] = (
            features[self.n_total_constants], features[f_j]
        )
        self.n_found_constants += 1
        self.n_total_constants += 1

    cdef inline void update_drawn_feature(self, SIZE_t f_j) nogil:
        """Move f_j into features that have been drawn and are not constant."""
        cdef SIZE_t[::1] features = self.features
        self.f_i -= 1
        features[self.f_i], features[f_j] = features[f_j], features[self.f_i]

    cdef inline SIZE_t update_constant_features(self) nogil:
        """Move constant features.

        Respect invariant for constant features: the original order of element in
        features[:n_known_constants] must be preserved for sibling and child nodes.
        """
        cdef SIZE_t[::1] features = self.features
        cdef SIZE_t[::1] constant_features = self.constant_features
        memcpy(&features[0], &constant_features[0], sizeof(SIZE_t) * self.n_known_constants)
        # Copy newly found constant features
        memcpy(&constant_features[self.n_known_constants], &features[self.n_known_constants],
               sizeof(SIZE_t) * self.n_found_constants)
