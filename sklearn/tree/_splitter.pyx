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
#
# License: BSD 3 clause

from ._criterion cimport Criterion

from libc.stdlib cimport qsort
from libc.string cimport memcpy
from cython cimport final

import numpy as np

from scipy.sparse import csc_matrix

from ._utils cimport log
from ._utils cimport rand_int
from ._utils cimport rand_uniform
from ._utils cimport RAND_R_MAX

cdef double INFINITY = np.inf

# Mitigate precision differences between 32 bit and 64 bit
cdef DTYPE_t FEATURE_THRESHOLD = 1e-7

# Constant to switch between algorithm non zero value extract algorithm
# in SparsePartitioner
cdef DTYPE_t EXTRACT_NNZ_SWITCH = 0.1

cdef inline void _init_split(SplitRecord* self, SIZE_t start_pos) noexcept nogil:
    self.impurity_left = INFINITY
    self.impurity_right = INFINITY
    self.pos = start_pos
    self.feature = 0
    self.threshold = 0.
    self.improvement = -INFINITY

cdef class Splitter:
    """Abstract splitter class.

    Splitters are called by tree builders to find the best splits on both
    sparse and dense data, one split at a time.
    """

    def __cinit__(self, Criterion criterion, SIZE_t max_features,
                  SIZE_t min_samples_leaf, double min_weight_leaf,
                  object random_state):
        """
        Parameters
        ----------
        criterion : Criterion
            The criterion to measure the quality of a split.

        max_features : SIZE_t
            The maximal number of randomly selected features which can be
            considered for a split.

        min_samples_leaf : SIZE_t
            The minimal number of samples each leaf can have, where splits
            which would result in having less samples in a leaf are not
            considered.

        min_weight_leaf : double
            The minimal weight each leaf can have, where the weight is the sum
            of the weights of each sample in it.

        random_state : object
            The user inputted random state to be used for pseudo-randomness
        """

        self.criterion = criterion

        self.n_samples = 0
        self.n_features = 0

        self.max_features = max_features
        self.min_samples_leaf = min_samples_leaf
        self.min_weight_leaf = min_weight_leaf
        self.random_state = random_state

    def __getstate__(self):
        return {}

    def __setstate__(self, d):
        pass

    def __reduce__(self):
        return (type(self), (self.criterion,
                             self.max_features,
                             self.min_samples_leaf,
                             self.min_weight_leaf,
                             self.random_state), self.__getstate__())

    cdef int init(
        self,
        object X,
        const DOUBLE_t[:, ::1] y,
        const DOUBLE_t[:] sample_weight
    ) except -1:
        """Initialize the splitter.

        Take in the input data X, the target Y, and optional sample weights.

        Returns -1 in case of failure to allocate memory (and raise MemoryError)
        or 0 otherwise.

        Parameters
        ----------
        X : object
            This contains the inputs. Usually it is a 2d numpy array.

        y : ndarray, dtype=DOUBLE_t
            This is the vector of targets, or true labels, for the samples represented
            as a Cython memoryview.

        sample_weight : ndarray, dtype=DOUBLE_t
            The weights of the samples, where higher weighted samples are fit
            closer than lower weight samples. If not provided, all samples
            are assumed to have uniform weight. This is represented
            as a Cython memoryview.
        """

        self.rand_r_state = self.random_state.randint(0, RAND_R_MAX)
        cdef SIZE_t n_samples = X.shape[0]

        # Create a new array which will be used to store nonzero
        # samples from the feature of interest
        self.samples = np.empty(n_samples, dtype=np.intp)
        cdef SIZE_t[::1] samples = self.samples

        cdef SIZE_t i, j
        cdef double weighted_n_samples = 0.0
        j = 0

        for i in range(n_samples):
            # Only work with positively weighted samples
            if sample_weight is None or sample_weight[i] != 0.0:
                samples[j] = i
                j += 1

            if sample_weight is not None:
                weighted_n_samples += sample_weight[i]
            else:
                weighted_n_samples += 1.0

        # Number of samples is number of positively weighted samples
        self.n_samples = j
        self.weighted_n_samples = weighted_n_samples

        cdef SIZE_t n_features = X.shape[1]
        self.features = np.arange(n_features, dtype=np.intp)
        self.n_features = n_features

        self.feature_values = np.empty(n_samples, dtype=np.float32)
        self.constant_features = np.empty(n_features, dtype=np.intp)

        self.y = y

        self.sample_weight = sample_weight
        return 0

    cdef int node_reset(self, SIZE_t start, SIZE_t end,
                        double* weighted_n_node_samples) except -1 nogil:
        """Reset splitter on node samples[start:end].

        Returns -1 in case of failure to allocate memory (and raise MemoryError)
        or 0 otherwise.

        Parameters
        ----------
        start : SIZE_t
            The index of the first sample to consider
        end : SIZE_t
            The index of the last sample to consider
        weighted_n_node_samples : ndarray, dtype=double pointer
            The total weight of those samples
        """

        self.start = start
        self.end = end

        self.criterion.init(
            self.y,
            self.sample_weight,
            self.weighted_n_samples,
            self.samples,
            start,
            end
        )

        weighted_n_node_samples[0] = self.criterion.weighted_n_node_samples
        return 0

    cdef int node_split(self, double impurity, SplitRecord* split,
                        SIZE_t* n_constant_features) except -1 nogil:
        """Find the best split on node samples[start:end].

        This is a placeholder method. The majority of computation will be done
        here.

        It should return -1 upon errors.
        """

        pass

    cdef void node_value(self, double* dest) noexcept nogil:
        """Copy the value of node samples[start:end] into dest."""

        self.criterion.node_value(dest)

    cdef double node_impurity(self) noexcept nogil:
        """Return the impurity of the current node."""

        return self.criterion.node_impurity()

# Introduce a fused-class to make it possible to share the split implementation
# between the dense and sparse cases in the node_split_best and node_split_random
# functions. The alternative would have been to use inheritance-based polymorphism
# but it would have resulted in a ~10% overall tree fitting performance
# degradation caused by the overhead frequent virtual method lookups.
ctypedef fused Partitioner:
    DensePartitioner
    SparsePartitioner

cdef inline int node_split_best(
    Splitter splitter,
    Partitioner partitioner,
    Criterion criterion,
    double impurity,
    SplitRecord* split,
    SIZE_t* n_constant_features,
) except -1 nogil:
    """Find the best split on node samples[start:end]

    Returns -1 in case of failure to allocate memory (and raise MemoryError)
    or 0 otherwise.
    """
    # Find the best split
    cdef SIZE_t start = splitter.start
    cdef SIZE_t end = splitter.end

    cdef SIZE_t[::1] features = splitter.features
    cdef SIZE_t[::1] constant_features = splitter.constant_features
    cdef SIZE_t n_features = splitter.n_features

    cdef DTYPE_t[::1] feature_values = splitter.feature_values
    cdef SIZE_t max_features = splitter.max_features
    cdef SIZE_t min_samples_leaf = splitter.min_samples_leaf
    cdef double min_weight_leaf = splitter.min_weight_leaf
    cdef UINT32_t* random_state = &splitter.rand_r_state

    cdef SplitRecord best_split, current_split
    cdef double current_proxy_improvement = -INFINITY
    cdef double best_proxy_improvement = -INFINITY

    cdef SIZE_t f_i = n_features
    cdef SIZE_t f_j
    cdef SIZE_t p
    cdef SIZE_t p_prev

    cdef SIZE_t n_visited_features = 0
    # Number of features discovered to be constant during the split search
    cdef SIZE_t n_found_constants = 0
    # Number of features known to be constant and drawn without replacement
    cdef SIZE_t n_drawn_constants = 0
    cdef SIZE_t n_known_constants = n_constant_features[0]
    # n_total_constants = n_known_constants + n_found_constants
    cdef SIZE_t n_total_constants = n_known_constants

    _init_split(&best_split, end)

    partitioner.init_node_split(start, end)

    # Sample up to max_features without replacement using a
    # Fisher-Yates-based algorithm (using the local variables `f_i` and
    # `f_j` to compute a permutation of the `features` array).
    #
    # Skip the CPU intensive evaluation of the impurity criterion for
    # features that were already detected as constant (hence not suitable
    # for good splitting) by ancestor nodes and save the information on
    # newly discovered constant features to spare computation on descendant
    # nodes.
    while (f_i > n_total_constants and  # Stop early if remaining features
                                        # are constant
            (n_visited_features < max_features or
             # At least one drawn features must be non constant
             n_visited_features <= n_found_constants + n_drawn_constants)):

        n_visited_features += 1

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
        f_j = rand_int(n_drawn_constants, f_i - n_found_constants,
                       random_state)

        if f_j < n_known_constants:
            # f_j in the interval [n_drawn_constants, n_known_constants[
            features[n_drawn_constants], features[f_j] = features[f_j], features[n_drawn_constants]

            n_drawn_constants += 1
            continue

        # f_j in the interval [n_known_constants, f_i - n_found_constants[
        f_j += n_found_constants
        # f_j in the interval [n_total_constants, f_i[
        current_split.feature = features[f_j]
        partitioner.sort_samples_and_feature_values(current_split.feature)

        if feature_values[end - 1] <= feature_values[start] + FEATURE_THRESHOLD:
            features[f_j], features[n_total_constants] = features[n_total_constants], features[f_j]

            n_found_constants += 1
            n_total_constants += 1
            continue

        f_i -= 1
        features[f_i], features[f_j] = features[f_j], features[f_i]

        # Evaluate all splits
        # At this point, the criterion has a view into the samples that was sorted
        # by the partitioner. The criterion will use that ordering when evaluating the splits.
        criterion.reset()
        p = start

        while p < end:
            partitioner.next_p(&p_prev, &p)

            if p >= end:
                continue

            current_split.pos = p

            # Reject if min_samples_leaf is not guaranteed
            if (((current_split.pos - start) < min_samples_leaf) or
                    ((end - current_split.pos) < min_samples_leaf)):
                continue

            criterion.update(current_split.pos)

            # Reject if min_weight_leaf is not satisfied
            if ((criterion.weighted_n_left < min_weight_leaf) or
                    (criterion.weighted_n_right < min_weight_leaf)):
                continue

            current_proxy_improvement = criterion.proxy_impurity_improvement()

            if current_proxy_improvement > best_proxy_improvement:
                best_proxy_improvement = current_proxy_improvement
                # sum of halves is used to avoid infinite value
                current_split.threshold = (
                    feature_values[p_prev] / 2.0 + feature_values[p] / 2.0
                )

                if (
                    current_split.threshold == feature_values[p] or
                    current_split.threshold == INFINITY or
                    current_split.threshold == -INFINITY
                ):
                    current_split.threshold = feature_values[p_prev]

                # This creates a SplitRecord copy
                best_split = current_split

    # Reorganize into samples[start:best_split.pos] + samples[best_split.pos:end]
    if best_split.pos < end:
        partitioner.partition_samples_final(
            best_split.pos,
            best_split.threshold,
            best_split.feature
        )
        criterion.reset()
        criterion.update(best_split.pos)
        criterion.children_impurity(
            &best_split.impurity_left, &best_split.impurity_right
        )
        best_split.improvement = criterion.impurity_improvement(
            impurity,
            best_split.impurity_left,
            best_split.impurity_right
        )

    # Respect invariant for constant features: the original order of
    # element in features[:n_known_constants] must be preserved for sibling
    # and child nodes
    memcpy(&features[0], &constant_features[0], sizeof(SIZE_t) * n_known_constants)

    # Copy newly found constant features
    memcpy(&constant_features[n_known_constants],
           &features[n_known_constants],
           sizeof(SIZE_t) * n_found_constants)

    # Return values
    split[0] = best_split
    n_constant_features[0] = n_total_constants
    return 0


# Sort n-element arrays pointed to by feature_values and samples, simultaneously,
# by the values in feature_values. Algorithm: Introsort (Musser, SP&E, 1997).
cdef inline void sort(DTYPE_t* feature_values, SIZE_t* samples, SIZE_t n) noexcept nogil:
    if n == 0:
      return
    cdef int maxd = 2 * <int>log(n)
    introsort(feature_values, samples, n, maxd)


cdef inline void swap(DTYPE_t* feature_values, SIZE_t* samples,
        SIZE_t i, SIZE_t j) noexcept nogil:
    # Helper for sort
    feature_values[i], feature_values[j] = feature_values[j], feature_values[i]
    samples[i], samples[j] = samples[j], samples[i]


cdef inline DTYPE_t median3(DTYPE_t* feature_values, SIZE_t n) noexcept nogil:
    # Median of three pivot selection, after Bentley and McIlroy (1993).
    # Engineering a sort function. SP&E. Requires 8/3 comparisons on average.
    cdef DTYPE_t a = feature_values[0], b = feature_values[n / 2], c = feature_values[n - 1]
    if a < b:
        if b < c:
            return b
        elif a < c:
            return c
        else:
            return a
    elif b < c:
        if a < c:
            return a
        else:
            return c
    else:
        return b


# Introsort with median of 3 pivot selection and 3-way partition function
# (robust to repeated elements, e.g. lots of zero features).
cdef void introsort(DTYPE_t* feature_values, SIZE_t *samples,
                    SIZE_t n, int maxd) noexcept nogil:
    cdef DTYPE_t pivot
    cdef SIZE_t i, l, r

    while n > 1:
        if maxd <= 0:   # max depth limit exceeded ("gone quadratic")
            heapsort(feature_values, samples, n)
            return
        maxd -= 1

        pivot = median3(feature_values, n)

        # Three-way partition.
        i = l = 0
        r = n
        while i < r:
            if feature_values[i] < pivot:
                swap(feature_values, samples, i, l)
                i += 1
                l += 1
            elif feature_values[i] > pivot:
                r -= 1
                swap(feature_values, samples, i, r)
            else:
                i += 1

        introsort(feature_values, samples, l, maxd)
        feature_values += r
        samples += r
        n -= r


cdef inline void sift_down(DTYPE_t* feature_values, SIZE_t* samples,
                           SIZE_t start, SIZE_t end) noexcept nogil:
    # Restore heap order in feature_values[start:end] by moving the max element to start.
    cdef SIZE_t child, maxind, root

    root = start
    while True:
        child = root * 2 + 1

        # find max of root, left child, right child
        maxind = root
        if child < end and feature_values[maxind] < feature_values[child]:
            maxind = child
        if child + 1 < end and feature_values[maxind] < feature_values[child + 1]:
            maxind = child + 1

        if maxind == root:
            break
        else:
            swap(feature_values, samples, root, maxind)
            root = maxind


cdef void heapsort(DTYPE_t* feature_values, SIZE_t* samples, SIZE_t n) noexcept nogil:
    cdef SIZE_t start, end

    # heapify
    start = (n - 2) / 2
    end = n
    while True:
        sift_down(feature_values, samples, start, end)
        if start == 0:
            break
        start -= 1

    # sort by shrinking the heap, putting the max element immediately after it
    end = n - 1
    while end > 0:
        swap(feature_values, samples, 0, end)
        sift_down(feature_values, samples, 0, end)
        end = end - 1

cdef inline int node_split_random(
    Splitter splitter,
    Partitioner partitioner,
    Criterion criterion,
    double impurity,
    SplitRecord* split,
    SIZE_t* n_constant_features
) except -1 nogil:
    """Find the best random split on node samples[start:end]

    Returns -1 in case of failure to allocate memory (and raise MemoryError)
    or 0 otherwise.
    """
    # Draw random splits and pick the best
    cdef SIZE_t start = splitter.start
    cdef SIZE_t end = splitter.end

    cdef SIZE_t[::1] features = splitter.features
    cdef SIZE_t[::1] constant_features = splitter.constant_features
    cdef SIZE_t n_features = splitter.n_features

    cdef SIZE_t max_features = splitter.max_features
    cdef SIZE_t min_samples_leaf = splitter.min_samples_leaf
    cdef double min_weight_leaf = splitter.min_weight_leaf
    cdef UINT32_t* random_state = &splitter.rand_r_state

    cdef SplitRecord best_split, current_split
    cdef double current_proxy_improvement = - INFINITY
    cdef double best_proxy_improvement = - INFINITY

    cdef SIZE_t f_i = n_features
    cdef SIZE_t f_j
    # Number of features discovered to be constant during the split search
    cdef SIZE_t n_found_constants = 0
    # Number of features known to be constant and drawn without replacement
    cdef SIZE_t n_drawn_constants = 0
    cdef SIZE_t n_known_constants = n_constant_features[0]
    # n_total_constants = n_known_constants + n_found_constants
    cdef SIZE_t n_total_constants = n_known_constants
    cdef SIZE_t n_visited_features = 0
    cdef DTYPE_t min_feature_value
    cdef DTYPE_t max_feature_value

    _init_split(&best_split, end)

    partitioner.init_node_split(start, end)

    # Sample up to max_features without replacement using a
    # Fisher-Yates-based algorithm (using the local variables `f_i` and
    # `f_j` to compute a permutation of the `features` array).
    #
    # Skip the CPU intensive evaluation of the impurity criterion for
    # features that were already detected as constant (hence not suitable
    # for good splitting) by ancestor nodes and save the information on
    # newly discovered constant features to spare computation on descendant
    # nodes.
    while (f_i > n_total_constants and  # Stop early if remaining features
                                        # are constant
            (n_visited_features < max_features or
             # At least one drawn features must be non constant
             n_visited_features <= n_found_constants + n_drawn_constants)):
        n_visited_features += 1

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
        f_j = rand_int(n_drawn_constants, f_i - n_found_constants,
                       random_state)

        if f_j < n_known_constants:
            # f_j in the interval [n_drawn_constants, n_known_constants[
            features[n_drawn_constants], features[f_j] = features[f_j], features[n_drawn_constants]
            n_drawn_constants += 1
            continue

        # f_j in the interval [n_known_constants, f_i - n_found_constants[
        f_j += n_found_constants
        # f_j in the interval [n_total_constants, f_i[

        current_split.feature = features[f_j]

        # Find min, max
        partitioner.find_min_max(
            current_split.feature, &min_feature_value, &max_feature_value
        )

        if max_feature_value <= min_feature_value + FEATURE_THRESHOLD:
            features[f_j], features[n_total_constants] = features[n_total_constants], current_split.feature

            n_found_constants += 1
            n_total_constants += 1
            continue

        f_i -= 1
        features[f_i], features[f_j] = features[f_j], features[f_i]

        # Draw a random threshold
        current_split.threshold = rand_uniform(
            min_feature_value,
            max_feature_value,
            random_state,
        )

        if current_split.threshold == max_feature_value:
            current_split.threshold = min_feature_value

        # Partition
        current_split.pos = partitioner.partition_samples(current_split.threshold)

        # Reject if min_samples_leaf is not guaranteed
        if (((current_split.pos - start) < min_samples_leaf) or
                ((end - current_split.pos) < min_samples_leaf)):
            continue

        # Evaluate split
        # At this point, the criterion has a view into the samples that was partitioned
        # by the partitioner. The criterion will use the parition to evaluating the split.
        criterion.reset()
        criterion.update(current_split.pos)

        # Reject if min_weight_leaf is not satisfied
        if ((criterion.weighted_n_left < min_weight_leaf) or
                (criterion.weighted_n_right < min_weight_leaf)):
            continue

        current_proxy_improvement = criterion.proxy_impurity_improvement()

        if current_proxy_improvement > best_proxy_improvement:
            best_proxy_improvement = current_proxy_improvement
            best_split = current_split  # copy

    # Reorganize into samples[start:best_split.pos] + samples[best_split.pos:end]
    if best_split.pos < end:
        if current_split.feature != best_split.feature:
            partitioner.partition_samples_final(
                best_split.pos, best_split.threshold, best_split.feature
            )

        criterion.reset()
        criterion.update(best_split.pos)
        criterion.children_impurity(
            &best_split.impurity_left, &best_split.impurity_right
        )
        best_split.improvement = criterion.impurity_improvement(
            impurity, best_split.impurity_left, best_split.impurity_right
        )

    # Respect invariant for constant features: the original order of
    # element in features[:n_known_constants] must be preserved for sibling
    # and child nodes
    memcpy(&features[0], &constant_features[0], sizeof(SIZE_t) * n_known_constants)

    # Copy newly found constant features
    memcpy(&constant_features[n_known_constants],
           &features[n_known_constants],
           sizeof(SIZE_t) * n_found_constants)

    # Return values
    split[0] = best_split
    n_constant_features[0] = n_total_constants
    return 0


@final
cdef class DensePartitioner:
    """Partitioner specialized for dense data.

    Note that this partitioner is agnostic to the splitting strategy (best vs. random).
    """
    cdef:
        const DTYPE_t[:, :] X
        cdef SIZE_t[::1] samples
        cdef DTYPE_t[::1] feature_values
        cdef SIZE_t start
        cdef SIZE_t end

    def __init__(
        self,
        const DTYPE_t[:, :] X,
        SIZE_t[::1] samples,
        DTYPE_t[::1] feature_values,
    ):
        self.X = X
        self.samples = samples
        self.feature_values = feature_values

    cdef inline void init_node_split(self, SIZE_t start, SIZE_t end) noexcept nogil:
        """Initialize splitter at the beginning of node_split."""
        self.start = start
        self.end = end

    cdef inline void sort_samples_and_feature_values(
        self, SIZE_t current_feature
    ) noexcept nogil:
        """Simultaneously sort based on the feature_values."""
        cdef:
            SIZE_t i
            DTYPE_t[::1] feature_values = self.feature_values
            const DTYPE_t[:, :] X = self.X
            SIZE_t[::1] samples = self.samples

        # Sort samples along that feature; by
        # copying the values into an array and
        # sorting the array in a manner which utilizes the cache more
        # effectively.
        for i in range(self.start, self.end):
            feature_values[i] = X[samples[i], current_feature]
        sort(&feature_values[self.start], &samples[self.start], self.end - self.start)

    cdef inline void find_min_max(
        self,
        SIZE_t current_feature,
        DTYPE_t* min_feature_value_out,
        DTYPE_t* max_feature_value_out,
    ) noexcept nogil:
        """Find the minimum and maximum value for current_feature."""
        cdef:
            SIZE_t p
            DTYPE_t current_feature_value
            const DTYPE_t[:, :] X = self.X
            SIZE_t[::1] samples = self.samples
            DTYPE_t min_feature_value = X[samples[self.start], current_feature]
            DTYPE_t max_feature_value = min_feature_value
            DTYPE_t[::1] feature_values = self.feature_values

        feature_values[self.start] = min_feature_value

        for p in range(self.start + 1, self.end):
            current_feature_value = X[samples[p], current_feature]
            feature_values[p] = current_feature_value

            if current_feature_value < min_feature_value:
                min_feature_value = current_feature_value
            elif current_feature_value > max_feature_value:
                max_feature_value = current_feature_value

        min_feature_value_out[0] = min_feature_value
        max_feature_value_out[0] = max_feature_value

    cdef inline void next_p(self, SIZE_t* p_prev, SIZE_t* p) noexcept nogil:
        """Compute the next p_prev and p for iteratiing over feature values."""
        cdef DTYPE_t[::1] feature_values = self.feature_values

        while (
            p[0] + 1 < self.end and
            feature_values[p[0] + 1] <= feature_values[p[0]] + FEATURE_THRESHOLD
        ):
            p[0] += 1

        p_prev[0] = p[0]

        # By adding 1, we have
        # (feature_values[p] >= end) or (feature_values[p] > feature_values[p - 1])
        p[0] += 1

    cdef inline SIZE_t partition_samples(self, double current_threshold) noexcept nogil:
        """Partition samples for feature_values at the current_threshold."""
        cdef:
            SIZE_t p = self.start
            SIZE_t partition_end = self.end
            SIZE_t[::1] samples = self.samples
            DTYPE_t[::1] feature_values = self.feature_values

        while p < partition_end:
            if feature_values[p] <= current_threshold:
                p += 1
            else:
                partition_end -= 1

                feature_values[p], feature_values[partition_end] = (
                    feature_values[partition_end], feature_values[p]
                )
                samples[p], samples[partition_end] = samples[partition_end], samples[p]

        return partition_end

    cdef inline void partition_samples_final(
        self,
        SIZE_t best_pos,
        double best_threshold,
        SIZE_t best_feature,
    ) noexcept nogil:
        """Partition samples for X at the best_threshold and best_feature."""
        cdef:
            SIZE_t p = self.start
            SIZE_t partition_end = self.end
            SIZE_t[::1] samples = self.samples
            const DTYPE_t[:, :] X = self.X

        while p < partition_end:
            if X[samples[p], best_feature] <= best_threshold:
                p += 1
            else:
                partition_end -= 1
                samples[p], samples[partition_end] = samples[partition_end], samples[p]

@final
cdef class SparsePartitioner:
    """Partitioner specialized for sparse CSC data.

    Note that this partitioner is agnostic to the splitting strategy (best vs. random).
    """
    cdef SIZE_t[::1] samples
    cdef DTYPE_t[::1] feature_values
    cdef SIZE_t start
    cdef SIZE_t end

    cdef const DTYPE_t[::1] X_data
    cdef const INT32_t[::1] X_indices
    cdef const INT32_t[::1] X_indptr

    cdef SIZE_t n_total_samples

    cdef SIZE_t[::1] index_to_samples
    cdef SIZE_t[::1] sorted_samples

    cdef SIZE_t start_positive
    cdef SIZE_t end_negative
    cdef bint is_samples_sorted

    def __init__(
        self,
        object X,
        SIZE_t[::1] samples,
        SIZE_t n_samples,
        DTYPE_t[::1] feature_values,
    ):
        if not isinstance(X, csc_matrix):
            raise ValueError("X should be in csc format")

        self.samples = samples
        self.feature_values = feature_values

        # Initialize X
        cdef SIZE_t n_total_samples = X.shape[0]

        self.X_data = X.data
        self.X_indices = X.indices
        self.X_indptr = X.indptr
        self.n_total_samples = n_total_samples

        # Initialize auxiliary array used to perform split
        self.index_to_samples = np.full(n_total_samples, fill_value=-1, dtype=np.intp)
        self.sorted_samples = np.empty(n_samples, dtype=np.intp)

        cdef SIZE_t p
        for p in range(n_samples):
            self.index_to_samples[samples[p]] = p

    cdef inline void init_node_split(self, SIZE_t start, SIZE_t end) noexcept nogil:
        """Initialize splitter at the beginning of node_split."""
        self.start = start
        self.end = end
        self.is_samples_sorted = 0

    cdef inline void sort_samples_and_feature_values(
        self, SIZE_t current_feature
    ) noexcept nogil:
        """Simultaneously sort based on the feature_values."""
        cdef:
            DTYPE_t[::1] feature_values = self.feature_values
            SIZE_t[::1] index_to_samples = self.index_to_samples
            SIZE_t[::1] samples = self.samples

        self.extract_nnz(current_feature)
        # Sort the positive and negative parts of `feature_values`
        sort(&feature_values[self.start], &samples[self.start], self.end_negative - self.start)
        if self.start_positive < self.end:
            sort(&feature_values[self.start_positive], &samples[self.start_positive],
                    self.end - self.start_positive)

        # Update index_to_samples to take into account the sort
        for p in range(self.start, self.end_negative):
            index_to_samples[samples[p]] = p
        for p in range(self.start_positive, self.end):
            index_to_samples[samples[p]] = p

        # Add one or two zeros in feature_values, if there is any
        if self.end_negative < self.start_positive:
            self.start_positive -= 1
            feature_values[self.start_positive] = 0.

            if self.end_negative != self.start_positive:
                feature_values[self.end_negative] = 0.
                self.end_negative += 1

    cdef inline void find_min_max(
        self,
        SIZE_t current_feature,
        DTYPE_t* min_feature_value_out,
        DTYPE_t* max_feature_value_out,
    ) noexcept nogil:
        """Find the minimum and maximum value for current_feature."""
        cdef:
            SIZE_t p
            DTYPE_t current_feature_value, min_feature_value, max_feature_value
            DTYPE_t[::1] feature_values = self.feature_values

        self.extract_nnz(current_feature)

        if self.end_negative != self.start_positive:
            # There is a zero
            min_feature_value = 0
            max_feature_value = 0
        else:
            min_feature_value = feature_values[self.start]
            max_feature_value = min_feature_value

        # Find min, max in feature_values[start:end_negative]
        for p in range(self.start, self.end_negative):
            current_feature_value = feature_values[p]

            if current_feature_value < min_feature_value:
                min_feature_value = current_feature_value
            elif current_feature_value > max_feature_value:
                max_feature_value = current_feature_value

        # Update min, max given feature_values[start_positive:end]
        for p in range(self.start_positive, self.end):
            current_feature_value = feature_values[p]

            if current_feature_value < min_feature_value:
                min_feature_value = current_feature_value
            elif current_feature_value > max_feature_value:
                max_feature_value = current_feature_value

        min_feature_value_out[0] = min_feature_value
        max_feature_value_out[0] = max_feature_value

    cdef inline void next_p(self, SIZE_t* p_prev, SIZE_t* p) noexcept nogil:
        """Compute the next p_prev and p for iteratiing over feature values."""
        cdef:
            SIZE_t p_next
            DTYPE_t[::1] feature_values = self.feature_values

        if p[0] + 1 != self.end_negative:
            p_next = p[0] + 1
        else:
            p_next = self.start_positive

        while (p_next < self.end and
                feature_values[p_next] <= feature_values[p[0]] + FEATURE_THRESHOLD):
            p[0] = p_next
            if p[0] + 1 != self.end_negative:
                p_next = p[0] + 1
            else:
                p_next = self.start_positive

        p_prev[0] = p[0]
        p[0] = p_next

    cdef inline SIZE_t partition_samples(self, double current_threshold) noexcept nogil:
        """Partition samples for feature_values at the current_threshold."""
        return self._partition(current_threshold, self.start_positive)

    cdef inline void partition_samples_final(
        self,
        SIZE_t best_pos,
        double best_threshold,
        SIZE_t best_feature,
    ) noexcept nogil:
        """Partition samples for X at the best_threshold and best_feature."""
        self.extract_nnz(best_feature)
        self._partition(best_threshold, best_pos)

    cdef inline SIZE_t _partition(self, double threshold, SIZE_t zero_pos) noexcept nogil:
        """Partition samples[start:end] based on threshold."""
        cdef:
            SIZE_t p, partition_end
            SIZE_t[::1] index_to_samples = self.index_to_samples
            DTYPE_t[::1] feature_values = self.feature_values
            SIZE_t[::1] samples = self.samples

        if threshold < 0.:
            p = self.start
            partition_end = self.end_negative
        elif threshold > 0.:
            p = self.start_positive
            partition_end = self.end
        else:
            # Data are already split
            return zero_pos

        while p < partition_end:
            if feature_values[p] <= threshold:
                p += 1

            else:
                partition_end -= 1

                feature_values[p], feature_values[partition_end] = (
                    feature_values[partition_end], feature_values[p]
                )
                sparse_swap(index_to_samples, samples, p, partition_end)

        return partition_end

    cdef inline void extract_nnz(self, SIZE_t feature) noexcept nogil:
        """Extract and partition values for a given feature.

        The extracted values are partitioned between negative values
        feature_values[start:end_negative[0]] and positive values
        feature_values[start_positive[0]:end].
        The samples and index_to_samples are modified according to this
        partition.

        The extraction corresponds to the intersection between the arrays
        X_indices[indptr_start:indptr_end] and samples[start:end].
        This is done efficiently using either an index_to_samples based approach
        or binary search based approach.

        Parameters
        ----------
        feature : SIZE_t,
            Index of the feature we want to extract non zero value.
        """
        cdef SIZE_t[::1] samples = self.samples
        cdef DTYPE_t[::1] feature_values = self.feature_values
        cdef SIZE_t indptr_start = self.X_indptr[feature],
        cdef SIZE_t indptr_end = self.X_indptr[feature + 1]
        cdef SIZE_t n_indices = <SIZE_t>(indptr_end - indptr_start)
        cdef SIZE_t n_samples = self.end - self.start
        cdef SIZE_t[::1] index_to_samples = self.index_to_samples
        cdef SIZE_t[::1] sorted_samples = self.sorted_samples
        cdef const INT32_t[::1] X_indices = self.X_indices
        cdef const DTYPE_t[::1] X_data = self.X_data

        # Use binary search if n_samples * log(n_indices) <
        # n_indices and index_to_samples approach otherwise.
        # O(n_samples * log(n_indices)) is the running time of binary
        # search and O(n_indices) is the running time of index_to_samples
        # approach.
        if ((1 - self.is_samples_sorted) * n_samples * log(n_samples) +
                n_samples * log(n_indices) < EXTRACT_NNZ_SWITCH * n_indices):
            extract_nnz_binary_search(X_indices, X_data,
                                      indptr_start, indptr_end,
                                      samples, self.start, self.end,
                                      index_to_samples,
                                      feature_values,
                                      &self.end_negative, &self.start_positive,
                                      sorted_samples, &self.is_samples_sorted)

        # Using an index to samples  technique to extract non zero values
        # index_to_samples is a mapping from X_indices to samples
        else:
            extract_nnz_index_to_samples(X_indices, X_data,
                                         indptr_start, indptr_end,
                                         samples, self.start, self.end,
                                         index_to_samples,
                                         feature_values,
                                         &self.end_negative, &self.start_positive)


cdef int compare_SIZE_t(const void* a, const void* b) noexcept nogil:
    """Comparison function for sort."""
    return <int>((<SIZE_t*>a)[0] - (<SIZE_t*>b)[0])


cdef inline void binary_search(const INT32_t[::1] sorted_array,
                               INT32_t start, INT32_t end,
                               SIZE_t value, SIZE_t* index,
                               INT32_t* new_start) noexcept nogil:
    """Return the index of value in the sorted array.

    If not found, return -1. new_start is the last pivot + 1
    """
    cdef INT32_t pivot
    index[0] = -1
    while start < end:
        pivot = start + (end - start) / 2

        if sorted_array[pivot] == value:
            index[0] = pivot
            start = pivot + 1
            break

        if sorted_array[pivot] < value:
            start = pivot + 1
        else:
            end = pivot
    new_start[0] = start


cdef inline void extract_nnz_index_to_samples(const INT32_t[::1] X_indices,
                                              const DTYPE_t[::1] X_data,
                                              INT32_t indptr_start,
                                              INT32_t indptr_end,
                                              SIZE_t[::1] samples,
                                              SIZE_t start,
                                              SIZE_t end,
                                              SIZE_t[::1] index_to_samples,
                                              DTYPE_t[::1] feature_values,
                                              SIZE_t* end_negative,
                                              SIZE_t* start_positive) noexcept nogil:
    """Extract and partition values for a feature using index_to_samples.

    Complexity is O(indptr_end - indptr_start).
    """
    cdef INT32_t k
    cdef SIZE_t index
    cdef SIZE_t end_negative_ = start
    cdef SIZE_t start_positive_ = end

    for k in range(indptr_start, indptr_end):
        if start <= index_to_samples[X_indices[k]] < end:
            if X_data[k] > 0:
                start_positive_ -= 1
                feature_values[start_positive_] = X_data[k]
                index = index_to_samples[X_indices[k]]
                sparse_swap(index_to_samples, samples, index, start_positive_)


            elif X_data[k] < 0:
                feature_values[end_negative_] = X_data[k]
                index = index_to_samples[X_indices[k]]
                sparse_swap(index_to_samples, samples, index, end_negative_)
                end_negative_ += 1

    # Returned values
    end_negative[0] = end_negative_
    start_positive[0] = start_positive_


cdef inline void extract_nnz_binary_search(const INT32_t[::1] X_indices,
                                           const DTYPE_t[::1] X_data,
                                           INT32_t indptr_start,
                                           INT32_t indptr_end,
                                           SIZE_t[::1] samples,
                                           SIZE_t start,
                                           SIZE_t end,
                                           SIZE_t[::1] index_to_samples,
                                           DTYPE_t[::1] feature_values,
                                           SIZE_t* end_negative,
                                           SIZE_t* start_positive,
                                           SIZE_t[::1] sorted_samples,
                                           bint* is_samples_sorted) noexcept nogil:
    """Extract and partition values for a given feature using binary search.

    If n_samples = end - start and n_indices = indptr_end - indptr_start,
    the complexity is

        O((1 - is_samples_sorted[0]) * n_samples * log(n_samples) +
          n_samples * log(n_indices)).
    """
    cdef SIZE_t n_samples

    if not is_samples_sorted[0]:
        n_samples = end - start
        memcpy(&sorted_samples[start], &samples[start],
               n_samples * sizeof(SIZE_t))
        qsort(&sorted_samples[start], n_samples, sizeof(SIZE_t),
              compare_SIZE_t)
        is_samples_sorted[0] = 1

    while (indptr_start < indptr_end and
           sorted_samples[start] > X_indices[indptr_start]):
        indptr_start += 1

    while (indptr_start < indptr_end and
           sorted_samples[end - 1] < X_indices[indptr_end - 1]):
        indptr_end -= 1

    cdef SIZE_t p = start
    cdef SIZE_t index
    cdef SIZE_t k
    cdef SIZE_t end_negative_ = start
    cdef SIZE_t start_positive_ = end

    while (p < end and indptr_start < indptr_end):
        # Find index of sorted_samples[p] in X_indices
        binary_search(X_indices, indptr_start, indptr_end,
                      sorted_samples[p], &k, &indptr_start)

        if k != -1:
             # If k != -1, we have found a non zero value

            if X_data[k] > 0:
                start_positive_ -= 1
                feature_values[start_positive_] = X_data[k]
                index = index_to_samples[X_indices[k]]
                sparse_swap(index_to_samples, samples, index, start_positive_)


            elif X_data[k] < 0:
                feature_values[end_negative_] = X_data[k]
                index = index_to_samples[X_indices[k]]
                sparse_swap(index_to_samples, samples, index, end_negative_)
                end_negative_ += 1
        p += 1

    # Returned values
    end_negative[0] = end_negative_
    start_positive[0] = start_positive_


cdef inline void sparse_swap(SIZE_t[::1] index_to_samples, SIZE_t[::1] samples,
                             SIZE_t pos_1, SIZE_t pos_2) noexcept nogil:
    """Swap sample pos_1 and pos_2 preserving sparse invariant."""
    samples[pos_1], samples[pos_2] =  samples[pos_2], samples[pos_1]
    index_to_samples[samples[pos_1]] = pos_1
    index_to_samples[samples[pos_2]] = pos_2


cdef class BestSplitter(Splitter):
    """Splitter for finding the best split on dense data."""
    cdef DensePartitioner partitioner
    cdef int init(
        self,
        object X,
        const DOUBLE_t[:, ::1] y,
        const DOUBLE_t[:] sample_weight
    ) except -1:
        Splitter.init(self, X, y, sample_weight)
        self.partitioner = DensePartitioner(X, self.samples, self.feature_values)

    cdef int node_split(self, double impurity, SplitRecord* split,
                        SIZE_t* n_constant_features) except -1 nogil:
        return node_split_best(
            self,
            self.partitioner,
            self.criterion,
            impurity,
            split,
            n_constant_features,
        )

cdef class BestSparseSplitter(Splitter):
    """Splitter for finding the best split, using the sparse data."""
    cdef SparsePartitioner partitioner
    cdef int init(
        self,
        object X,
        const DOUBLE_t[:, ::1] y,
        const DOUBLE_t[:] sample_weight
    ) except -1:
        Splitter.init(self, X, y, sample_weight)
        self.partitioner = SparsePartitioner(
            X, self.samples, self.n_samples, self.feature_values
        )

    cdef int node_split(self, double impurity, SplitRecord* split,
                        SIZE_t* n_constant_features) except -1 nogil:
        return node_split_best(
            self,
            self.partitioner,
            self.criterion,
            impurity,
            split,
            n_constant_features,
        )

cdef class RandomSplitter(Splitter):
    """Splitter for finding the best random split on dense data."""
    cdef DensePartitioner partitioner
    cdef int init(
        self,
        object X,
        const DOUBLE_t[:, ::1] y,
        const DOUBLE_t[:] sample_weight
    ) except -1:
        Splitter.init(self, X, y, sample_weight)
        self.partitioner = DensePartitioner(X, self.samples, self.feature_values)

    cdef int node_split(self, double impurity, SplitRecord* split,
                        SIZE_t* n_constant_features) except -1 nogil:
        return node_split_random(
            self,
            self.partitioner,
            self.criterion,
            impurity,
            split,
            n_constant_features,
        )

cdef class RandomSparseSplitter(Splitter):
    """Splitter for finding the best random split, using the sparse data."""
    cdef SparsePartitioner partitioner
    cdef int init(
        self,
        object X,
        const DOUBLE_t[:, ::1] y,
        const DOUBLE_t[:] sample_weight
    ) except -1:
        Splitter.init(self, X, y, sample_weight)
        self.partitioner = SparsePartitioner(
            X, self.samples, self.n_samples, self.feature_values
        )

    cdef int node_split(self, double impurity, SplitRecord* split,
                        SIZE_t* n_constant_features) except -1 nogil:
        return node_split_random(
            self,
            self.partitioner,
            self.criterion,
            impurity,
            split,
            n_constant_features,
        )
