#distutils: language=c++
#cython: language_level=3
#cython: boundscheck=False
#cython: wraparound=False
#cython: profile=True

cimport cython
import numpy as np
cimport numpy as np
np.import_array()

from ._criterion cimport Criterion

from libc.stdlib cimport malloc
from libc.stdlib cimport free
from libc.stdlib cimport qsort
from libc.string cimport memcpy
from libc.string cimport memset
from libc.stdio cimport printf
from libcpp.vector cimport vector

# TODO: allow sparse operations
from scipy.sparse import csc_matrix

from cython.operator cimport dereference as deref
from cython.parallel import prange

from ._utils cimport log
from ._utils cimport rand_int
from ._utils cimport rand_uniform
from ._utils cimport RAND_R_MAX
from ._utils cimport safe_realloc

cdef double INFINITY = np.inf

# Mitigate precision differences between 32 bit and 64 bit
cdef DTYPE_t FEATURE_THRESHOLD = 1e-7

# Constant to switch between algorithm non zero value extract algorithm
# in SparseSplitter
cdef DTYPE_t EXTRACT_NNZ_SWITCH = 0.1


cdef inline void _init_split(ObliqueSplitRecord* self, SIZE_t start_pos) nogil:
    self.impurity_left = INFINITY
    self.impurity_right = INFINITY
    self.pos = start_pos
    self.feature = 0
    self.threshold = 0.
    self.improvement = -INFINITY

cdef class ObliqueSplitter(Splitter):
    """Abstract oblique splitter class.

    Splitters are called by tree builders to find the best splits on 
    both sparse and dense data, one split at a time.
    """

    def __cinit__(self, Criterion criterion, SIZE_t max_features,
                  SIZE_t min_samples_leaf, double min_weight_leaf,
                  double feature_combinations, object random_state):
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
        
        feature_combinations : double
            The average number of features to combine in an oblique split. 
            Each feature is independently included with probability
            ``feature_combination`` / ``n_features``.

        random_state : object
            The user inputted random state to be used for pseudo-randomness
        """
        self.criterion = criterion

        self.samples = NULL
        self.n_samples = 0
        self.features = NULL
        self.n_features = 0
        self.feature_values = NULL

        self.sample_weight = NULL

        # Max features = output dimensionality of projection vectors
        self.max_features = max_features
        self.min_samples_leaf = min_samples_leaf
        self.min_weight_leaf = min_weight_leaf
        self.random_state = random_state

        # Oblique tree parameteres
        self.feature_combinations = feature_combinations

        # Sparse max_features x n_features projection matrix
        self.proj_mat_weights = vector[vector[DTYPE_t]](self.max_features)
        self.proj_mat_indices = vector[vector[SIZE_t]](self.max_features)

        # or max w/ 1...
        self.n_non_zeros = max(int(self.max_features * self.feature_combinations), 1)

    def __getstate__(self):
        return {}

    def __setstate__(self, d):
        pass

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
                            self.samples,
                            start,
                            end)

        weighted_n_node_samples[0] = self.criterion.weighted_n_node_samples

        # Clear all projection vectors
        for i in range(self.max_features):
            self.proj_mat_weights[i].clear()
            self.proj_mat_indices[i].clear()

    cdef void sample_proj_mat(self, 
                              vector[vector[DTYPE_t]]& proj_mat_weights, 
                              vector[vector[SIZE_t]]& proj_mat_indices) nogil:
        """ Sample the projection vector. 
        
        This is a placeholder method. 

        """

        pass

    cdef int pointer_size(self) nogil:
        """Get size of a pointer to record for ObliqueSplitter."""

        return sizeof(ObliqueSplitRecord)

cdef class BaseDenseObliqueSplitter(ObliqueSplitter):
    
    cdef const DTYPE_t[:, :] X

    def __cinit__(self, Criterion criterion, SIZE_t max_features,
                  SIZE_t min_samples_leaf, double min_weight_leaf,
                  double feature_combinations,
                  object random_state):
        self.max_features = max_features # number of proj_vecs
        self.feature_combinations = feature_combinations

    cdef int init(self,
                  object X,
                  const DOUBLE_t[:, ::1] y,
                  DOUBLE_t* sample_weight) except -1:
        """Initialize the splitter

        Returns -1 in case of failure to allocate memory (and raise MemoryError)
        or 0 otherwise.
        """

        # Call parent init
        ObliqueSplitter.init(self, X, y, sample_weight)

        self.X = X

cdef class BestObliqueSplitter(BaseDenseObliqueSplitter):
    def __reduce__(self):
        """Enable pickling the splitter."""
        return (BestObliqueSplitter, (self.criterion,
                               self.max_features,
                               self.min_samples_leaf,
                               self.min_weight_leaf,
                               self.feature_combinations,
                               self.random_state), self.__getstate__())

    # NOTE: vectors are passed by value, so & is needed to pass by reference
    cdef void sample_proj_mat(self, 
                              vector[vector[DTYPE_t]]& proj_mat_weights,
                              vector[vector[SIZE_t]]& proj_mat_indices) nogil:
        """
        Sparse Oblique Projection matrix.
        Randomly sample features to put in randomly sampled projection vectors
        weight = 1 or -1 with probability 0.5 
        """
 
        cdef SIZE_t n_features = self.n_features
        cdef SIZE_t max_features = self.max_features
        cdef SIZE_t n_non_zeros = self.n_non_zeros
        cdef UINT32_t* random_state = &self.rand_r_state

        cdef int i, feat_i, proj_i
        cdef DTYPE_t weight

        for i in range(n_non_zeros):
            proj_i = rand_int(0, max_features, random_state)
            feat_i = rand_int(0, n_features, random_state)
            weight = 1 if (rand_int(0, 2, random_state) == 1) else -1

            proj_mat_indices[proj_i].push_back(feat_i)  # Store index of nonzero
            proj_mat_weights[proj_i].push_back(weight)  # Store weight of nonzero

    cdef int node_split(self, double impurity, SplitRecord* split,
                        SIZE_t* n_constant_features) nogil except -1:
        """Find the best split on node samples[start:end]

        Returns -1 in case of failure to allocate memory (and raise MemoryError)
        or 0 otherwise.
        """
        # typecast the pointer to an ObliqueSplitRecord
        cdef ObliqueSplitRecord* oblique_split = <ObliqueSplitRecord*>(split)

        # Draw random splits and pick the best
        cdef SIZE_t* samples = self.samples
        cdef SIZE_t start = self.start
        cdef SIZE_t end = self.end

        cdef SIZE_t* features = self.features
        cdef SIZE_t* constant_features = self.constant_features
        cdef SIZE_t n_features = self.n_features

        # pointer array to store feature values to split on
        cdef DTYPE_t* Xf = self.feature_values
        cdef SIZE_t max_features = self.max_features
        cdef SIZE_t min_samples_leaf = self.min_samples_leaf
        cdef double min_weight_leaf = self.min_weight_leaf
        cdef UINT32_t* random_state = &self.rand_r_state

        # keep track of split record for current node and the best split
        # found among the sampled projection vectors
        cdef ObliqueSplitRecord best, current
        cdef double current_proxy_improvement = -INFINITY
        cdef double best_proxy_improvement = -INFINITY

        cdef SIZE_t feat_i, p       # index over computed features and start/end
        cdef SIZE_t idx, jdx        # index over max_feature, and 
        cdef SIZE_t partition_end
        cdef DTYPE_t temp_d         # to compute a projection feature value

        # instantiate the split records
        _init_split(&best, end)

        # Sample the projection matrix
        self.sample_proj_mat(self.proj_mat_weights, self.proj_mat_indices)

        # define vector pointers to the projection weights and indices
        cdef vector[DTYPE_t]* proj_vec_weights
        cdef vector[SIZE_t]* proj_vec_indices

        # For every vector in the projection matrix
        for feat_i in range(max_features):
            # Projection vector has no nonzeros
            if self.proj_mat_weights[feat_i].empty():
                continue

            # XXX: 'feature' is not actually used in oblique split records
            # Just indicates which split was sampled
            current.feature = feat_i
            current.proj_vec_weights = &self.proj_mat_weights[feat_i]
            current.proj_vec_indices = &self.proj_mat_indices[feat_i]

            # fill memory block with 0's for Xf from start to end with DTYPE_t
            memset(Xf + start, 0, (end - start) * sizeof(DTYPE_t))

            # Compute linear combination of features and then
            # sort samples according to the feature values.
            for idx in range(start, end):
                for jdx in range(0, current.proj_vec_indices.size()):
                    Xf[idx] += self.X[samples[idx], deref(current.proj_vec_indices)[jdx]] * deref(current.proj_vec_weights)[jdx]

            # Sort the samples
            sort(Xf + start, samples + start, end - start)

            # Evaluate all splits
            self.criterion.reset()
            p = start
            while p < end:
                while (p + 1 < end and
                        Xf[p + 1] <= Xf[p] + FEATURE_THRESHOLD):
                    p += 1

                p += 1

                if p < end:
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

                        if ((current.threshold == Xf[p]) or
                            (current.threshold == INFINITY) or
                            (current.threshold == -INFINITY)):
                            current.threshold = Xf[p - 1]

                        best = current  # copy

        # Reorganize into samples[start:best.pos] + samples[best.pos:end]
        if best.pos < end:
            partition_end = end
            p = start

            while p < partition_end:
                
                # Account for projection vector
                temp_d = 0.0
                for j in range(best.proj_vec_indices.size()):
                    temp_d += self.X[samples[p], deref(best.proj_vec_indices)[j]] * deref(best.proj_vec_weights)[j]

                if temp_d <= best.threshold:
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

        # Return values
        deref(oblique_split).proj_vec_indices = best.proj_vec_indices
        deref(oblique_split).proj_vec_weights = best.proj_vec_weights
        deref(oblique_split).feature = best.feature
        deref(oblique_split).pos = best.pos
        deref(oblique_split).threshold = best.threshold
        deref(oblique_split).improvement = best.improvement
        deref(oblique_split).impurity_left = best.impurity_left
        deref(oblique_split).impurity_right = best.impurity_right

        # n_constant_features[0] = n_total_constants
        return 0

cdef class RandomObliqueSplitter(BaseDenseObliqueSplitter):
    """Splitter for finding the best random split."""
    def __reduce__(self):
        return (RandomObliqueSplitter, (self.criterion,
                                 self.max_features,
                                 self.min_samples_leaf,
                                 self.min_weight_leaf,
                                 self.feature_combinations,
                                 self.random_state), self.__getstate__())

    cdef int node_split(self, double impurity, SplitRecord* split,
                        SIZE_t* n_constant_features) nogil except -1:
        """Find the best random split on node samples[start:end]

        Returns -1 in case of failure to allocate memory (and raise MemoryError)
        or 0 otherwise.
        """
        # typecast the pointer to an ObliqueSplitRecord
        cdef ObliqueSplitRecord* oblique_split = <ObliqueSplitRecord*>(split)

        # Draw random splits and pick the best
        cdef SIZE_t* samples = self.samples
        cdef SIZE_t start = self.start
        cdef SIZE_t end = self.end

        cdef SIZE_t* features = self.features
        cdef SIZE_t* constant_features = self.constant_features
        cdef SIZE_t n_features = self.n_features

        # pointer array to store feature values to split on
        cdef DTYPE_t* Xf = self.feature_values
        cdef SIZE_t max_features = self.max_features
        cdef SIZE_t min_samples_leaf = self.min_samples_leaf
        cdef double min_weight_leaf = self.min_weight_leaf
        cdef UINT32_t* random_state = &self.rand_r_state

        # keep track of split record for current node and the best split
        # found among the sampled projection vectors
        cdef ObliqueSplitRecord best, current
        cdef double current_proxy_improvement = -INFINITY
        cdef double best_proxy_improvement = -INFINITY

        cdef SIZE_t feat_i, p       # index over computed features and start/end
        cdef SIZE_t idx, jdx        # index over max_feature, and 
        cdef SIZE_t partition_end
        cdef DTYPE_t temp_d         # to compute a projection feature value
        cdef DTYPE_t min_feature_value
        cdef DTYPE_t max_feature_value
        cdef DTYPE_t current_feature_value

        # instantiate the split records
        _init_split(&best, end)

        # Sample the projection matrix
        self.sample_proj_mat(self.proj_mat_weights, self.proj_mat_indices)

        # define vector pointers to the projection weights and indices
        cdef vector[DTYPE_t]* proj_vec_weights
        cdef vector[SIZE_t]* proj_vec_indices

        # For every vector in the projection matrix
        for feat_i in range(max_features):
            # Projection vector has no nonzeros
            if self.proj_mat_weights[feat_i].empty():
                continue

            # XXX: 'feature' is not actually used in oblique split records
            # Just indicates which split was sampled
            current.feature = feat_i
            current.proj_vec_weights = &self.proj_mat_weights[feat_i]
            current.proj_vec_indices = &self.proj_mat_indices[feat_i]

            # fill memory block with 0's for Xf from start to end with DTYPE_t
            memset(Xf + start, 0, (end - start) * sizeof(DTYPE_t))

            # Find min and max
            # First, compute the projection for the first sample
            for jdx in range(0, current.proj_vec_indices.size()):
                Xf[start] += self.X[samples[start], deref(current.proj_vec_indices)[jdx]] * deref(current.proj_vec_weights)[jdx]
            min_feature_value = Xf[start]
            max_feature_value = min_feature_value

            # Second, compute the projection for the rest of the samples
            # and then keep track of the minimum and maximum values
            for p in range(start + 1, end):
                for jdx in range(0, current.proj_vec_indices.size()):
                    Xf[p] += self.X[samples[p], deref(current.proj_vec_indices)[jdx]] * deref(current.proj_vec_weights)[jdx]

                current_feature_value = Xf[p]
                Xf[p] = current_feature_value

                if current_feature_value < min_feature_value:
                    min_feature_value = current_feature_value
                elif current_feature_value > max_feature_value:
                    max_feature_value = current_feature_value

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
                    # Account for projection vector
                    temp_d = 0.0
                    for j in range(best.proj_vec_indices.size()):
                        temp_d += self.X[samples[p], deref(best.proj_vec_indices)[j]] * deref(best.proj_vec_weights)[j]

                    if temp_d <= best.threshold:
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

        # Return values
        deref(oblique_split).proj_vec_indices = best.proj_vec_indices
        deref(oblique_split).proj_vec_weights = best.proj_vec_weights
        deref(oblique_split).feature = best.feature
        deref(oblique_split).pos = best.pos
        deref(oblique_split).threshold = best.threshold
        deref(oblique_split).improvement = best.improvement
        deref(oblique_split).impurity_left = best.impurity_left
        deref(oblique_split).impurity_right = best.impurity_right

        return 0