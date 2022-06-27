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

cimport cython
from ._criterion cimport Criterion

from libc.stdlib cimport free
from libc.stdlib cimport qsort
from libc.string cimport memcpy
from libc.string cimport memset

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
# in SparseSplitter
cdef DTYPE_t EXTRACT_NNZ_SWITCH = 0.1

@cython.final
cdef class FeatureTracker:
    """Feature Sampler using Fisher-Yates-based algorithm.

    Sample up to max_features without replacement using a
    Fisher-Yates-based algorithm (using the local variables `f_i` and
    `f_j` to compute a permutation of the `features` array).

    Skip the CPU intensive evaluation of the impurity criterion for
    features that were already detected as constant (hence not suitable
    for good splitting) by ancestor nodes and save the information on
    newly discovered constant features to spare computation on descendant
    nodes.
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

        This function can return a struct with a feature index, it's value,
        and status.
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

    cdef inline void is_constant(self, SIZE_t f_j) nogil:
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

        Respect invariant for constant features: the original order of
        element in features[:n_known_constants] must be preserved for sibling
        and child nodes.
        """
        cdef:
            SIZE_t[::1] features = self.features
            SIZE_t[::1] constant_features = self.constant_features
        memcpy(&features[0], &constant_features[0],
               sizeof(SIZE_t) * self.n_known_constants)

        # Copy newly found constant features
        memcpy(&constant_features[self.n_known_constants],
               &features[self.n_known_constants],
               sizeof(SIZE_t) * self.n_found_constants)


cdef inline void _init_split(SplitRecord* self, SIZE_t start_pos) nogil:
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

        self.sample_weight = NULL

        self.max_features = max_features
        self.min_samples_leaf = min_samples_leaf
        self.min_weight_leaf = min_weight_leaf
        self.random_state = random_state

    def __getstate__(self):
        return {}

    def __setstate__(self, d):
        pass

    cdef int init(self,
                   object X,
                   const DOUBLE_t[:, ::1] y,
                   DOUBLE_t* sample_weight) except -1:
        """Initialize the splitter.

        Take in the input data X, the target Y, and optional sample weights.

        Returns -1 in case of failure to allocate memory (and raise MemoryError)
        or 0 otherwise.

        Parameters
        ----------
        X : object
            This contains the inputs. Usually it is a 2d numpy array.

        y : ndarray, dtype=DOUBLE_t
            This is the vector of targets, or true labels, for the samples

        sample_weight : DOUBLE_t*
            The weights of the samples, where higher weighted samples are fit
            closer than lower weight samples. If not provided, all samples
            are assumed to have uniform weight.
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
            if sample_weight == NULL or sample_weight[i] != 0.0:
                samples[j] = i
                j += 1

            if sample_weight != NULL:
                weighted_n_samples += sample_weight[i]
            else:
                weighted_n_samples += 1.0

        # Number of samples is number of positively weighted samples
        self.n_samples = j
        self.weighted_n_samples = weighted_n_samples

        cdef SIZE_t n_features = X.shape[1]
        self.feature_tracker = FeatureTracker(X.shape[1], self.max_features)
        self.feature_values = np.empty(n_samples, dtype=np.float32)

        self.y = y

        self.sample_weight = sample_weight
        return 0

    cdef int node_reset(self, SIZE_t start, SIZE_t end,
                        double* weighted_n_node_samples) nogil except -1:
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

        self.criterion.init(self.y,
                            self.sample_weight,
                            self.weighted_n_samples,
                            &self.samples[0],
                            start,
                            end)

        weighted_n_node_samples[0] = self.criterion.weighted_n_node_samples
        return 0

    cdef int node_split(self, double impurity, SplitRecord* split,
                        SIZE_t* n_constant_features) nogil except -1:
        """Find the best split on node samples[start:end].

        This is a placeholder method. The majority of computation will be done
        here.

        It should return -1 upon errors.
        """

        pass

    cdef void node_value(self, double* dest) nogil:
        """Copy the value of node samples[start:end] into dest."""

        self.criterion.node_value(dest)

    cdef double node_impurity(self) nogil:
        """Return the impurity of the current node."""

        return self.criterion.node_impurity()


cdef class BaseDenseSplitter(Splitter):
    cdef const DTYPE_t[:, :] X

    cdef SIZE_t n_total_samples

    cdef int init(self,
                  object X,
                  const DOUBLE_t[:, ::1] y,
                  DOUBLE_t* sample_weight) except -1:
        """Initialize the splitter

        Returns -1 in case of failure to allocate memory (and raise MemoryError)
        or 0 otherwise.
        """

        # Call parent init
        Splitter.init(self, X, y, sample_weight)

        self.X = X
        return 0


cdef class BestSplitter(BaseDenseSplitter):
    """Splitter for finding the best split."""
    def __reduce__(self):
        return (BestSplitter, (self.criterion,
                               self.max_features,
                               self.min_samples_leaf,
                               self.min_weight_leaf,
                               self.random_state), self.__getstate__())

    cdef int node_split(self, double impurity, SplitRecord* split,
                        SIZE_t* n_constant_features) nogil except -1:
        """Find the best split on node samples[start:end]

        Returns -1 in case of failure to allocate memory (and raise MemoryError)
        or 0 otherwise.
        """
        # Find the best split
        cdef SIZE_t[::1] samples = self.samples
        cdef SIZE_t start = self.start
        cdef SIZE_t end = self.end

        cdef DTYPE_t[::1] Xf = self.feature_values
        cdef SIZE_t max_features = self.max_features
        cdef SIZE_t min_samples_leaf = self.min_samples_leaf
        cdef double min_weight_leaf = self.min_weight_leaf
        cdef UINT32_t* random_state = &self.rand_r_state

        cdef SplitRecord best, current
        cdef double current_proxy_improvement = -INFINITY
        cdef double best_proxy_improvement = -INFINITY

        cdef FeatureSample feature_sample
        cdef SIZE_t f_j
        cdef SIZE_t p
        cdef SIZE_t feature_idx_offset
        cdef SIZE_t feature_offset
        cdef SIZE_t i
        cdef SIZE_t j

        cdef SIZE_t partition_end

        _init_split(&best, end)
        self.feature_tracker.reset(n_constant_features[0])

        while True:
            feature_sample = self.feature_tracker.sample(random_state)
            if feature_sample.status == FeatureStatus.STOP:
                break
            if feature_sample.status == FeatureStatus.CONTINUE:
                continue

            f_j = feature_sample.f_j
            current.feature = feature_sample.feature

            # Sort samples along that feature; by
            # copying the values into an array and
            # sorting the array in a manner which utilizes the cache more
            # effectively.
            for i in range(start, end):
                Xf[i] = self.X[samples[i], current.feature]

            sort(&Xf[start], &samples[start], end - start)

            if Xf[end - 1] <= Xf[start] + FEATURE_THRESHOLD:
                self.feature_tracker.is_constant(f_j)
                continue

            self.feature_tracker.update_drawn_feature(f_j)

            # Evaluate all splits
            self.criterion.reset()
            p = start

            while p < end:
                while p + 1 < end and Xf[p + 1] <= Xf[p] + FEATURE_THRESHOLD:
                    p += 1

                # (p + 1 >= end) or (X[samples[p + 1], current.feature] >
                #                    X[samples[p], current.feature])
                p += 1
                # (p >= end) or (X[samples[p], current.feature] >
                #                X[samples[p - 1], current.feature])

                if p >= end:
                    continue

                current.pos = p

                # Reject if min_samples_leaf is not guaranteed
                if (((current.pos - start) < min_samples_leaf) or
                        ((end - current.pos) < min_samples_leaf)):
                    continue

                self.criterion.update(current.pos)

                # Reject if min_weight_leaf is not satisfied
                if ((self.criterion.weighted_n_left < min_weight_leaf) or
                        (self.criterion.weighted_n_right < min_weight_leaf)):
                    continue

                current_proxy_improvement = self.criterion.proxy_impurity_improvement()

                if current_proxy_improvement > best_proxy_improvement:
                    best_proxy_improvement = current_proxy_improvement
                    # sum of halves is used to avoid infinite value
                    current.threshold = Xf[p - 1] / 2.0 + Xf[p] / 2.0

                    if (
                        current.threshold == Xf[p] or
                        current.threshold == INFINITY or
                        current.threshold == -INFINITY
                    ):
                        current.threshold = Xf[p - 1]

                    best = current  # copy

        # Reorganize into samples[start:best.pos] + samples[best.pos:end]
        if best.pos < end:
            partition_end = end
            p = start

            while p < partition_end:
                if self.X[samples[p], best.feature] <= best.threshold:
                    p += 1

                else:
                    partition_end -= 1

                    samples[p], samples[partition_end] = samples[partition_end], samples[p]

            self.criterion.reset()
            self.criterion.update(best.pos)
            self.criterion.children_impurity(&best.impurity_left,
                                             &best.impurity_right)
            best.improvement = self.criterion.impurity_improvement(
                impurity, best.impurity_left, best.impurity_right)

        self.feature_tracker.update_constant_features()

        # Return values
        split[0] = best
        n_constant_features[0] = self.feature_tracker.n_total_constants
        return 0


# Sort n-element arrays pointed to by Xf and samples, simultaneously,
# by the values in Xf. Algorithm: Introsort (Musser, SP&E, 1997).
cdef inline void sort(DTYPE_t* Xf, SIZE_t* samples, SIZE_t n) nogil:
    if n == 0:
      return
    cdef int maxd = 2 * <int>log(n)
    introsort(Xf, samples, n, maxd)


cdef inline void swap(DTYPE_t* Xf, SIZE_t* samples,
        SIZE_t i, SIZE_t j) nogil:
    # Helper for sort
    Xf[i], Xf[j] = Xf[j], Xf[i]
    samples[i], samples[j] = samples[j], samples[i]


cdef inline DTYPE_t median3(DTYPE_t* Xf, SIZE_t n) nogil:
    # Median of three pivot selection, after Bentley and McIlroy (1993).
    # Engineering a sort function. SP&E. Requires 8/3 comparisons on average.
    cdef DTYPE_t a = Xf[0], b = Xf[n / 2], c = Xf[n - 1]
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
cdef void introsort(DTYPE_t* Xf, SIZE_t *samples,
                    SIZE_t n, int maxd) nogil:
    cdef DTYPE_t pivot
    cdef SIZE_t i, l, r

    while n > 1:
        if maxd <= 0:   # max depth limit exceeded ("gone quadratic")
            heapsort(Xf, samples, n)
            return
        maxd -= 1

        pivot = median3(Xf, n)

        # Three-way partition.
        i = l = 0
        r = n
        while i < r:
            if Xf[i] < pivot:
                swap(Xf, samples, i, l)
                i += 1
                l += 1
            elif Xf[i] > pivot:
                r -= 1
                swap(Xf, samples, i, r)
            else:
                i += 1

        introsort(Xf, samples, l, maxd)
        Xf += r
        samples += r
        n -= r


cdef inline void sift_down(DTYPE_t* Xf, SIZE_t* samples,
                           SIZE_t start, SIZE_t end) nogil:
    # Restore heap order in Xf[start:end] by moving the max element to start.
    cdef SIZE_t child, maxind, root

    root = start
    while True:
        child = root * 2 + 1

        # find max of root, left child, right child
        maxind = root
        if child < end and Xf[maxind] < Xf[child]:
            maxind = child
        if child + 1 < end and Xf[maxind] < Xf[child + 1]:
            maxind = child + 1

        if maxind == root:
            break
        else:
            swap(Xf, samples, root, maxind)
            root = maxind


cdef void heapsort(DTYPE_t* Xf, SIZE_t* samples, SIZE_t n) nogil:
    cdef SIZE_t start, end

    # heapify
    start = (n - 2) / 2
    end = n
    while True:
        sift_down(Xf, samples, start, end)
        if start == 0:
            break
        start -= 1

    # sort by shrinking the heap, putting the max element immediately after it
    end = n - 1
    while end > 0:
        swap(Xf, samples, 0, end)
        sift_down(Xf, samples, 0, end)
        end = end - 1


cdef class RandomSplitter(BaseDenseSplitter):
    """Splitter for finding the best random split."""
    def __reduce__(self):
        return (RandomSplitter, (self.criterion,
                                 self.max_features,
                                 self.min_samples_leaf,
                                 self.min_weight_leaf,
                                 self.random_state), self.__getstate__())

    cdef int node_split(self, double impurity, SplitRecord* split,
                        SIZE_t* n_constant_features) nogil except -1:
        """Find the best random split on node samples[start:end]

        Returns -1 in case of failure to allocate memory (and raise MemoryError)
        or 0 otherwise.
        """
        # Draw random splits and pick the best
        cdef SIZE_t[::1] samples = self.samples
        cdef SIZE_t start = self.start
        cdef SIZE_t end = self.end

        cdef DTYPE_t[::1] Xf = self.feature_values
        cdef SIZE_t max_features = self.max_features
        cdef SIZE_t min_samples_leaf = self.min_samples_leaf
        cdef double min_weight_leaf = self.min_weight_leaf
        cdef UINT32_t* random_state = &self.rand_r_state

        cdef SplitRecord best, current
        cdef double current_proxy_improvement = - INFINITY
        cdef double best_proxy_improvement = - INFINITY

        cdef FeatureSample feature_sample
        cdef SIZE_t f_j
        cdef SIZE_t p
        cdef SIZE_t partition_end
        cdef SIZE_t feature_stride
        cdef DTYPE_t min_feature_value
        cdef DTYPE_t max_feature_value
        cdef DTYPE_t current_feature_value

        _init_split(&best, end)
        self.feature_tracker.reset(n_constant_features[0])

        while True:
            feature_sample = self.feature_tracker.sample(random_state)
            if feature_sample.status == FeatureStatus.STOP:
                break
            if feature_sample.status == FeatureStatus.CONTINUE:
                continue

            f_j = feature_sample.f_j
            current.feature = feature_sample.feature

            # Find min, max
            min_feature_value = self.X[samples[start], current.feature]
            max_feature_value = min_feature_value
            Xf[start] = min_feature_value

            for p in range(start + 1, end):
                current_feature_value = self.X[samples[p], current.feature]
                Xf[p] = current_feature_value

                if current_feature_value < min_feature_value:
                    min_feature_value = current_feature_value
                elif current_feature_value > max_feature_value:
                    max_feature_value = current_feature_value

            if max_feature_value <= min_feature_value + FEATURE_THRESHOLD:
                self.feature_tracker.is_constant(f_j)
                continue

            self.feature_tracker.update_drawn_feature(f_j)

            # Draw a random threshold
            current.threshold = rand_uniform(min_feature_value,
                                             max_feature_value,
                                             random_state)

            if current.threshold == max_feature_value:
                current.threshold = min_feature_value

            # Partition
            p, partition_end = start, end
            while p < partition_end:
                if Xf[p] <= current.threshold:
                    p += 1
                else:
                    partition_end -= 1

                    Xf[p], Xf[partition_end] = Xf[partition_end], Xf[p]
                    samples[p], samples[partition_end] = samples[partition_end], samples[p]

            current.pos = partition_end

            # Reject if min_samples_leaf is not guaranteed
            if (((current.pos - start) < min_samples_leaf) or
                    ((end - current.pos) < min_samples_leaf)):
                continue

            # Evaluate split
            self.criterion.reset()
            self.criterion.update(current.pos)

            # Reject if min_weight_leaf is not satisfied
            if ((self.criterion.weighted_n_left < min_weight_leaf) or
                    (self.criterion.weighted_n_right < min_weight_leaf)):
                continue

            current_proxy_improvement = self.criterion.proxy_impurity_improvement()

            if current_proxy_improvement > best_proxy_improvement:
                best_proxy_improvement = current_proxy_improvement
                best = current  # copy

        # Reorganize into samples[start:best.pos] + samples[best.pos:end]
        if best.pos < end:
            if current.feature != best.feature:
                p, partition_end = start, end

                while p < partition_end:
                    if self.X[samples[p], best.feature] <= best.threshold:
                        p += 1
                    else:
                        partition_end -= 1

                        samples[p], samples[partition_end] = samples[partition_end], samples[p]

            self.criterion.reset()
            self.criterion.update(best.pos)
            self.criterion.children_impurity(&best.impurity_left,
                                             &best.impurity_right)
            best.improvement = self.criterion.impurity_improvement(
                impurity, best.impurity_left, best.impurity_right)

        self.feature_tracker.update_constant_features()

        # Return values
        split[0] = best
        n_constant_features[0] = self.feature_tracker.n_total_constants
        return 0


cdef class BaseSparseSplitter(Splitter):
    # The sparse splitter works only with csc sparse matrix format
    cdef DTYPE_t[::1] X_data
    cdef INT32_t[::1] X_indices
    cdef INT32_t[::1] X_indptr

    cdef SIZE_t n_total_samples

    cdef SIZE_t[::1] index_to_samples
    cdef SIZE_t[::1] sorted_samples

    def __cinit__(self, Criterion criterion, SIZE_t max_features,
                  SIZE_t min_samples_leaf, double min_weight_leaf,
                  object random_state):
        # Parent __cinit__ is automatically called
        self.n_total_samples = 0

    cdef int init(self,
                  object X,
                  const DOUBLE_t[:, ::1] y,
                  DOUBLE_t* sample_weight) except -1:
        """Initialize the splitter

        Returns -1 in case of failure to allocate memory (and raise MemoryError)
        or 0 otherwise.
        """
        # Call parent init
        Splitter.init(self, X, y, sample_weight)

        if not isinstance(X, csc_matrix):
            raise ValueError("X should be in csc format")

        cdef SIZE_t[::1] samples = self.samples
        cdef SIZE_t n_samples = self.n_samples

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
        return 0

    cdef inline SIZE_t _partition(self, double threshold,
                                  SIZE_t end_negative, SIZE_t start_positive,
                                  SIZE_t zero_pos) nogil:
        """Partition samples[start:end] based on threshold."""

        cdef SIZE_t p
        cdef SIZE_t partition_end

        cdef DTYPE_t[::1] Xf = self.feature_values
        cdef SIZE_t[::1] samples = self.samples
        cdef SIZE_t[::1] index_to_samples = self.index_to_samples

        if threshold < 0.:
            p = self.start
            partition_end = end_negative
        elif threshold > 0.:
            p = start_positive
            partition_end = self.end
        else:
            # Data are already split
            return zero_pos

        while p < partition_end:
            if Xf[p] <= threshold:
                p += 1

            else:
                partition_end -= 1

                Xf[p], Xf[partition_end] = Xf[partition_end], Xf[p]
                sparse_swap(index_to_samples, samples, p, partition_end)

        return partition_end

    cdef inline void extract_nnz(self, SIZE_t feature,
                                 SIZE_t* end_negative, SIZE_t* start_positive,
                                 bint* is_samples_sorted) nogil:
        """Extract and partition values for a given feature.

        The extracted values are partitioned between negative values
        Xf[start:end_negative[0]] and positive values Xf[start_positive[0]:end].
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


        end_negative, start_positive : SIZE_t*, SIZE_t*,
            Return extracted non zero values in self.samples[start:end] where
            negative values are in self.feature_values[start:end_negative[0]]
            and positive values are in
            self.feature_values[start_positive[0]:end].

        is_samples_sorted : bint*,
            If is_samples_sorted, then self.sorted_samples[start:end] will be
            the sorted version of self.samples[start:end].

        """
        cdef SIZE_t indptr_start = self.X_indptr[feature],
        cdef SIZE_t indptr_end = self.X_indptr[feature + 1]
        cdef SIZE_t n_indices = <SIZE_t>(indptr_end - indptr_start)
        cdef SIZE_t n_samples = self.end - self.start
        cdef SIZE_t[::1] samples = self.samples
        cdef DTYPE_t[::1] feature_values = self.feature_values
        cdef SIZE_t[::1] index_to_samples = self.index_to_samples
        cdef SIZE_t[::1] sorted_samples = self.sorted_samples
        cdef INT32_t[::1] X_indices = self.X_indices
        cdef DTYPE_t[::1] X_data = self.X_data

        # Use binary search if n_samples * log(n_indices) <
        # n_indices and index_to_samples approach otherwise.
        # O(n_samples * log(n_indices)) is the running time of binary
        # search and O(n_indices) is the running time of index_to_samples
        # approach.
        if ((1 - is_samples_sorted[0]) * n_samples * log(n_samples) +
                n_samples * log(n_indices) < EXTRACT_NNZ_SWITCH * n_indices):
            extract_nnz_binary_search(X_indices, X_data,
                                      indptr_start, indptr_end,
                                      samples, self.start, self.end,
                                      index_to_samples,
                                      feature_values,
                                      end_negative, start_positive,
                                      sorted_samples, is_samples_sorted)

        # Using an index to samples  technique to extract non zero values
        # index_to_samples is a mapping from X_indices to samples
        else:
            extract_nnz_index_to_samples(X_indices, X_data,
                                         indptr_start, indptr_end,
                                         samples, self.start, self.end,
                                         index_to_samples,
                                         feature_values,
                                         end_negative, start_positive)


cdef int compare_SIZE_t(const void* a, const void* b) nogil:
    """Comparison function for sort."""
    return <int>((<SIZE_t*>a)[0] - (<SIZE_t*>b)[0])


cdef inline void binary_search(INT32_t[::1] sorted_array,
                               INT32_t start, INT32_t end,
                               SIZE_t value, SIZE_t* index,
                               INT32_t* new_start) nogil:
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


cdef inline void extract_nnz_index_to_samples(INT32_t[::1] X_indices,
                                              DTYPE_t[::1] X_data,
                                              INT32_t indptr_start,
                                              INT32_t indptr_end,
                                              SIZE_t[::1] samples,
                                              SIZE_t start,
                                              SIZE_t end,
                                              SIZE_t[::1] index_to_samples,
                                              DTYPE_t[::1] Xf,
                                              SIZE_t* end_negative,
                                              SIZE_t* start_positive) nogil:
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
                Xf[start_positive_] = X_data[k]
                index = index_to_samples[X_indices[k]]
                sparse_swap(index_to_samples, samples, index, start_positive_)


            elif X_data[k] < 0:
                Xf[end_negative_] = X_data[k]
                index = index_to_samples[X_indices[k]]
                sparse_swap(index_to_samples, samples, index, end_negative_)
                end_negative_ += 1

    # Returned values
    end_negative[0] = end_negative_
    start_positive[0] = start_positive_


cdef inline void extract_nnz_binary_search(INT32_t[::1] X_indices,
                                           DTYPE_t[::1] X_data,
                                           INT32_t indptr_start,
                                           INT32_t indptr_end,
                                           SIZE_t[::1] samples,
                                           SIZE_t start,
                                           SIZE_t end,
                                           SIZE_t[::1] index_to_samples,
                                           DTYPE_t[::1] Xf,
                                           SIZE_t* end_negative,
                                           SIZE_t* start_positive,
                                           SIZE_t[::1] sorted_samples,
                                           bint* is_samples_sorted) nogil:
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
                Xf[start_positive_] = X_data[k]
                index = index_to_samples[X_indices[k]]
                sparse_swap(index_to_samples, samples, index, start_positive_)


            elif X_data[k] < 0:
                Xf[end_negative_] = X_data[k]
                index = index_to_samples[X_indices[k]]
                sparse_swap(index_to_samples, samples, index, end_negative_)
                end_negative_ += 1
        p += 1

    # Returned values
    end_negative[0] = end_negative_
    start_positive[0] = start_positive_


cdef inline void sparse_swap(SIZE_t[::1] index_to_samples, SIZE_t[::1] samples,
                             SIZE_t pos_1, SIZE_t pos_2) nogil:
    """Swap sample pos_1 and pos_2 preserving sparse invariant."""
    samples[pos_1], samples[pos_2] =  samples[pos_2], samples[pos_1]
    index_to_samples[samples[pos_1]] = pos_1
    index_to_samples[samples[pos_2]] = pos_2


cdef class BestSparseSplitter(BaseSparseSplitter):
    """Splitter for finding the best split, using the sparse data."""

    def __reduce__(self):
        return (BestSparseSplitter, (self.criterion,
                                     self.max_features,
                                     self.min_samples_leaf,
                                     self.min_weight_leaf,
                                     self.random_state), self.__getstate__())

    cdef int node_split(self, double impurity, SplitRecord* split,
                        SIZE_t* n_constant_features) nogil except -1:
        """Find the best split on node samples[start:end], using sparse features

        Returns -1 in case of failure to allocate memory (and raise MemoryError)
        or 0 otherwise.
        """
        # Find the best split
        cdef SIZE_t[::1] samples = self.samples
        cdef SIZE_t start = self.start
        cdef SIZE_t end = self.end

        cdef DTYPE_t[::1] Xf = self.feature_values
        cdef SIZE_t[::1] index_to_samples = self.index_to_samples
        cdef SIZE_t max_features = self.max_features
        cdef SIZE_t min_samples_leaf = self.min_samples_leaf
        cdef double min_weight_leaf = self.min_weight_leaf
        cdef UINT32_t* random_state = &self.rand_r_state

        cdef SplitRecord best, current
        _init_split(&best, end)
        cdef double current_proxy_improvement = - INFINITY
        cdef double best_proxy_improvement = - INFINITY

        cdef FeatureSample feature_sample
        cdef SIZE_t f_j, p
        cdef DTYPE_t current_feature_value

        cdef SIZE_t p_next
        cdef SIZE_t p_prev
        cdef bint is_samples_sorted = 0  # indicate is sorted_samples is
                                         # inititialized

        # We assume implicitly that end_positive = end and
        # start_negative = start
        cdef SIZE_t start_positive
        cdef SIZE_t end_negative

        self.feature_tracker.reset(n_constant_features[0])

        while True:
            feature_sample = self.feature_tracker.sample(random_state)
            if feature_sample.status == FeatureStatus.STOP:
                break
            if feature_sample.status == FeatureStatus.CONTINUE:
                continue

            f_j = feature_sample.f_j
            current.feature = feature_sample.feature

            self.extract_nnz(current.feature, &end_negative, &start_positive,
                             &is_samples_sorted)
            # Sort the positive and negative parts of `Xf`
            sort(&Xf[start], &samples[start], end_negative - start)
            if start_positive < end:
                sort(&Xf[start_positive], &samples[start_positive],
                     end - start_positive)

            # Update index_to_samples to take into account the sort
            for p in range(start, end_negative):
                index_to_samples[samples[p]] = p
            for p in range(start_positive, end):
                index_to_samples[samples[p]] = p

            # Add one or two zeros in Xf, if there is any
            if end_negative < start_positive:
                start_positive -= 1
                Xf[start_positive] = 0.

                if end_negative != start_positive:
                    Xf[end_negative] = 0.
                    end_negative += 1

            if Xf[end - 1] <= Xf[start] + FEATURE_THRESHOLD:
                self.feature_tracker.is_constant(f_j)
                continue

            self.feature_tracker.update_drawn_feature(f_j)

            # Evaluate all splits
            self.criterion.reset()
            p = start

            while p < end:
                if p + 1 != end_negative:
                    p_next = p + 1
                else:
                    p_next = start_positive

                while (p_next < end and
                        Xf[p_next] <= Xf[p] + FEATURE_THRESHOLD):
                    p = p_next
                    if p + 1 != end_negative:
                        p_next = p + 1
                    else:
                        p_next = start_positive


                # (p_next >= end) or (X[samples[p_next], current.feature] >
                #                     X[samples[p], current.feature])
                p_prev = p
                p = p_next
                # (p >= end) or (X[samples[p], current.feature] >
                #                X[samples[p_prev], current.feature])

                if p >= end:
                    continue

                current.pos = p

                # Reject if min_samples_leaf is not guaranteed
                if (((current.pos - start) < min_samples_leaf) or
                        ((end - current.pos) < min_samples_leaf)):
                    continue

                self.criterion.update(current.pos)

                # Reject if min_weight_leaf is not satisfied
                if ((self.criterion.weighted_n_left < min_weight_leaf) or
                        (self.criterion.weighted_n_right < min_weight_leaf)):
                    continue

                current_proxy_improvement = self.criterion.proxy_impurity_improvement()

                if current_proxy_improvement > best_proxy_improvement:
                    best_proxy_improvement = current_proxy_improvement
                    # sum of halves used to avoid infinite values
                    current.threshold = Xf[p_prev] / 2.0 + Xf[p] / 2.0

                    if (
                        current.threshold == Xf[p] or
                        current.threshold == INFINITY or
                        current.threshold == -INFINITY
                    ):
                        current.threshold = Xf[p_prev]

                    best = current

        # Reorganize into samples[start:best.pos] + samples[best.pos:end]
        if best.pos < end:
            self.extract_nnz(best.feature, &end_negative, &start_positive,
                             &is_samples_sorted)

            self._partition(best.threshold, end_negative, start_positive,
                            best.pos)

            self.criterion.reset()
            self.criterion.update(best.pos)
            self.criterion.children_impurity(&best.impurity_left,
                                             &best.impurity_right)
            best.improvement = self.criterion.impurity_improvement(
                impurity, best.impurity_left, best.impurity_right)

        self.feature_tracker.update_constant_features()

        # Return values
        split[0] = best
        n_constant_features[0] = self.feature_tracker.n_total_constants
        return 0


cdef class RandomSparseSplitter(BaseSparseSplitter):
    """Splitter for finding a random split, using the sparse data."""

    def __reduce__(self):
        return (RandomSparseSplitter, (self.criterion,
                                       self.max_features,
                                       self.min_samples_leaf,
                                       self.min_weight_leaf,
                                       self.random_state), self.__getstate__())

    cdef int node_split(self, double impurity, SplitRecord* split,
                        SIZE_t* n_constant_features) nogil except -1:
        """Find a random split on node samples[start:end], using sparse features

        Returns -1 in case of failure to allocate memory (and raise MemoryError)
        or 0 otherwise.
        """
        # Find the best split
        cdef SIZE_t[::1] samples = self.samples
        cdef SIZE_t start = self.start
        cdef SIZE_t end = self.end

        cdef DTYPE_t[::1] Xf = self.feature_values
        cdef SIZE_t[::1] index_to_samples = self.index_to_samples
        cdef SIZE_t max_features = self.max_features
        cdef SIZE_t min_samples_leaf = self.min_samples_leaf
        cdef double min_weight_leaf = self.min_weight_leaf
        cdef UINT32_t* random_state = &self.rand_r_state

        cdef SplitRecord best, current
        _init_split(&best, end)
        cdef double current_proxy_improvement = - INFINITY
        cdef double best_proxy_improvement = - INFINITY

        cdef DTYPE_t current_feature_value

        cdef FeatureSample feature_sample
        cdef SIZE_t f_j, p
        cdef SIZE_t partition_end

        cdef DTYPE_t min_feature_value
        cdef DTYPE_t max_feature_value

        cdef bint is_samples_sorted = 0  # indicate that sorted_samples is
                                         # inititialized

        # We assume implicitly that end_positive = end and
        # start_negative = start
        cdef SIZE_t start_positive
        cdef SIZE_t end_negative

        self.feature_tracker.reset(n_constant_features[0])
        while True:
            feature_sample = self.feature_tracker.sample(random_state)
            if feature_sample.status == FeatureStatus.STOP:
                break
            if feature_sample.status == FeatureStatus.CONTINUE:
                continue

            f_j = feature_sample.f_j
            current.feature = feature_sample.feature

            self.extract_nnz(current.feature,
                             &end_negative, &start_positive,
                             &is_samples_sorted)

            if end_negative != start_positive:
                # There is a zero
                min_feature_value = 0
                max_feature_value = 0
            else:
                min_feature_value = Xf[start]
                max_feature_value = min_feature_value

            # Find min, max in Xf[start:end_negative]
            for p in range(start, end_negative):
                current_feature_value = Xf[p]

                if current_feature_value < min_feature_value:
                    min_feature_value = current_feature_value
                elif current_feature_value > max_feature_value:
                    max_feature_value = current_feature_value

            # Update min, max given Xf[start_positive:end]
            for p in range(start_positive, end):
                current_feature_value = Xf[p]

                if current_feature_value < min_feature_value:
                    min_feature_value = current_feature_value
                elif current_feature_value > max_feature_value:
                    max_feature_value = current_feature_value

            if max_feature_value <= min_feature_value + FEATURE_THRESHOLD:
                self.feature_tracker.is_constant(f_j)
                continue

            self.feature_tracker.update_drawn_feature(f_j)

            # Draw a random threshold
            current.threshold = rand_uniform(min_feature_value,
                                             max_feature_value,
                                             random_state)

            if current.threshold == max_feature_value:
                current.threshold = min_feature_value

            # Partition
            current.pos = self._partition(current.threshold,
                                          end_negative,
                                          start_positive,
                                          start_positive)

            # Reject if min_samples_leaf is not guaranteed
            if (((current.pos - start) < min_samples_leaf) or
                    ((end - current.pos) < min_samples_leaf)):
                continue

            # Evaluate split
            self.criterion.reset()
            self.criterion.update(current.pos)

            # Reject if min_weight_leaf is not satisfied
            if ((self.criterion.weighted_n_left < min_weight_leaf) or
                    (self.criterion.weighted_n_right < min_weight_leaf)):
                continue

            current_proxy_improvement = self.criterion.proxy_impurity_improvement()

            if current_proxy_improvement > best_proxy_improvement:
                best_proxy_improvement = current_proxy_improvement
                self.criterion.children_impurity(&current.impurity_left,
                                                 &current.impurity_right)
                current.improvement = self.criterion.impurity_improvement(
                    impurity, current.impurity_left, current.impurity_right)
                best = current

        # Reorganize into samples[start:best.pos] + samples[best.pos:end]
        if best.pos < end:
            if current.feature != best.feature:
                self.extract_nnz(best.feature, &end_negative, &start_positive,
                                 &is_samples_sorted)

                self._partition(best.threshold, end_negative, start_positive,
                                best.pos)

            self.criterion.reset()
            self.criterion.update(best.pos)
            self.criterion.children_impurity(&best.impurity_left,
                                             &best.impurity_right)
            best.improvement = self.criterion.impurity_improvement(
                impurity, best.impurity_left, best.impurity_right)

        self.feature_tracker.update_constant_features()

        # Return values
        split[0] = best
        n_constant_features[0] = self.feature_tracker.n_total_constants
        return 0
