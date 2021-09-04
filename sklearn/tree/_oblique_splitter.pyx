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

# allow sparse operations
# from scipy.sparse import csc_matrixfrom ._criterion cimport Criterion
# from scipy.sparse import csc_matrix

from cython.operator cimport dereference as deref
from cython.parallel import prange

from ._splitter cimport (sort, INFINITY, FEATURE_THRESHOLD, EXTRACT_NNZ_SWITCH)

from ._utils cimport log
from ._utils cimport rand_int
from ._utils cimport rand_uniform
from ._utils cimport RAND_R_MAX
from ._utils cimport safe_realloc


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

        # SPORF parameters
        self.feature_combinations = feature_combinations

        # Sparse max_features x n_features projection matrix
        self.proj_mat_weights = vector[vector[DTYPE_t]](self.max_features)
        self.proj_mat_indices = vector[vector[SIZE_t]](self.max_features)

        self.n_non_zeros = max(int(self.max_features * self.feature_combinations), 1)

    def __dealloc__(self):
        """Destructor."""

        free(self.samples)
        # print("freed samples")

        free(self.features)
        # print("freed features")

        free(self.constant_features)
        # print("freed constant_features")

        free(self.feature_values)
        # print("freed feature_values")

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
        # print('inside split init...')
        self.rand_r_state = self.random_state.randint(0, RAND_R_MAX)
        cdef SIZE_t n_samples = X.shape[0]

        # Create a new array which will be used to store nonzero
        # samples from the feature of interest
        cdef SIZE_t* samples = safe_realloc(&self.samples, n_samples)

        cdef SIZE_t i, j
        cdef double weighted_n_samples = 0.0
        j = 0

        # print('Initializing sample weights...')
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
        cdef SIZE_t* features = safe_realloc(&self.features, n_features)

        for i in range(n_features):
            features[i] = i

        self.n_features = n_features

        # print('About to safe realloc...')
        safe_realloc(&self.feature_values, n_samples)
        # print('Safe reallocated feature values..')
        safe_realloc(&self.constant_features, n_features)
        # print('After second...')
        self.y = y

        self.sample_weight = sample_weight
        # print('Finished init...')

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
                            self.samples,
                            start,
                            end)

        weighted_n_node_samples[0] = self.criterion.weighted_n_node_samples

        # Clear all projection vectors
        for i in range(self.max_features):
            self.proj_mat_weights[i].clear()
            self.proj_mat_indices[i].clear()

    cdef int oblique_node_split(self, double impurity, ObliqueSplitRecord* split,
                                SIZE_t* n_constant_features) nogil except -1:
        """Find the best split on node samples[start:end].

        This is a placeholder method. The majority of computation will be done
        here.

        It should return -1 upon errors.
        """
        pass

    cdef void sample_proj_mat(self,
                              vector[vector[DTYPE_t]]& proj_mat_weights,
                              vector[vector[SIZE_t]]& proj_mat_indices) nogil:
        """ Sample the projection vector from the set of features.

        This is a placeholder method.
        """
        pass


cdef class DenseObliqueSplitter(ObliqueSplitter):
    cdef const DTYPE_t[:, :] X

    cdef np.ndarray X_idx_sorted
    cdef INT32_t* X_idx_sorted_ptr
    cdef SIZE_t X_idx_sorted_stride
    cdef SIZE_t n_total_samples
    cdef SIZE_t* sample_mask

    def __cinit__(self, Criterion criterion, SIZE_t max_features,
                  SIZE_t min_samples_leaf, double min_weight_leaf,
                  double feature_combinations,
                  object random_state):
        self.X_idx_sorted_ptr = NULL
        self.X_idx_sorted_stride = 0
        self.sample_mask = NULL
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


cdef class BestObliqueSplitter(DenseObliqueSplitter):
    def __reduce__(self):
        """Enable pickling the splitter."""
        return (ObliqueSplitter, (self.criterion,
                               self.max_features,
                               self.min_samples_leaf,
                               self.min_weight_leaf,
                               self.feature_combinations,
                               self.random_state), self.__getstate__())

    cdef void sample_proj_mat(self,
                              vector[vector[DTYPE_t]]& proj_mat_weights,
                              vector[vector[SIZE_t]]& proj_mat_indices) nogil:
        """
        Oblique sparse projection matrix.

        Randomly sample features to put in randomly sampled projection vectors
        weight = 1 or -1 with probability 0.5

        NOTE: vectors are passed by value, so & is needed to pass by reference.
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

    cdef int oblique_node_split(self, double impurity, ObliqueSplitRecord* split,
                                SIZE_t* n_constant_features) nogil except -1:
        """Find the best split on node samples[start:end]

        Returns -1 in case of failure to allocate memory (and raise MemoryError)
        or 0 otherwise.
        """

        cdef SIZE_t* samples = self.samples
        cdef SIZE_t start = self.start
        cdef SIZE_t end = self.end

        cdef SIZE_t* features = self.features
        cdef SIZE_t* constant_features = self.constant_features
        cdef SIZE_t n_features = self.n_features

        cdef DTYPE_t* Xf = self.feature_values
        cdef SIZE_t max_features = self.max_features
        cdef SIZE_t min_samples_leaf = self.min_samples_leaf
        cdef double min_weight_leaf = self.min_weight_leaf
        cdef UINT32_t* random_state = &self.rand_r_state

        cdef INT32_t* X_idx_sorted = self.X_idx_sorted_ptr
        cdef SIZE_t* sample_mask = self.sample_mask

        # keep track of split record for current node and the best split
        # found among the sampled projection vectors
        cdef ObliqueSplitRecord best, current

        cdef double current_proxy_improvement = -INFINITY
        cdef double best_proxy_improvement = -INFINITY

        cdef SIZE_t f, p, i, j
        cdef SIZE_t partition_end
        cdef DTYPE_t temp_d

        # instantiate the split records
        _init_split(&best, end)

        # Sample the projection matrix
        self.sample_proj_mat(self.proj_mat_weights, self.proj_mat_indices)

        # cdef DTYPE_t** X_proj = self.X_proj
        cdef vector[DTYPE_t]* proj_vec_weights
        cdef vector[SIZE_t]* proj_vec_indices

        # For every vector in the projection matrix
        for f in range(max_features):

            # Projection vector has no nonzeros
            if self.proj_mat_weights[f].empty():
                continue

            current.feature = f
            current.proj_vec_weights = &self.proj_mat_weights[f]
            current.proj_vec_indices = &self.proj_mat_indices[f]

            # Compute linear combination of features
            memset(Xf + start, 0, (end - start) * sizeof(DTYPE_t))
            for i in range(start, end):
                for j in range(0, current.proj_vec_indices.size()):
                    Xf[i] += self.X[samples[i], deref(current.proj_vec_indices)[j]] * deref(current.proj_vec_weights)[j]

            # Sort the samples
            sort(Xf + start, samples + start, end - start)

            # Evaluate all splits
            self.criterion.reset()
            for p in range(start, end-1):

                # invalid split
                if (Xf[p + 1] <= Xf[p] + FEATURE_THRESHOLD):
                    continue

                current.pos = p

                # reject if min_samples_leaf not guaranteed
                if ((current.pos - start) < min_samples_leaf or
                    (end - current.pos) < min_samples_leaf):
                    continue

                self.criterion.update(current.pos)

                # reject if min_weight_leaf not satisfied
                if (self.criterion.weighted_n_left < min_weight_leaf or
                    self.criterion.weighted_n_right < min_weight_leaf):
                    continue

                current_proxy_improvement = self.criterion.proxy_impurity_improvement()

                if current_proxy_improvement > best_proxy_improvement:
                    best_proxy_improvement = current_proxy_improvement
                    # sum of halves is used to avoid infinite value
                    current.threshold = Xf[p - 1] / 2.0 + Xf[p] / 2.0

                    if (current.threshold == Xf[p] or
                        current.threshold == INFINITY or
                        current.threshold == -INFINITY):

                        current.threshold = Xf[p-1]

                    best = current

        # Reorganize into samples[start:best.pos] + samples[best.pos:end]
        if best.pos < end:
            partition_end = end
            p = start

            while p < partition_end:

                # Account for projection vector
                temp_d = 0
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
        split[0] = best
        # n_constant_features[0] = n_total_constants
        return 0
