#cython: language_level=3
#cython: boundscheck=False
#cython: wraparound=False

cimport cython
import numpy as np
cimport numpy as np
np.import_array()

from ._criterion cimport Criterion
# from sklearn.tree._criterion cimport Criterion

from libc.stdlib cimport malloc
from libc.stdlib cimport free
from libc.stdlib cimport qsort
from libc.string cimport memcpy
from libc.string cimport memset
from libc.stdio cimport printf

# allow sparse operations
# from scipy.sparse import csc_matrixfrom ._criterion cimport Criterion
# from scipy.sparse import csc_matrix

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
# cdef DTYPE_t EXTRACT_NNZ_SWITCH = 0.1

cdef void _init_split(ObliqueSplitRecord* self, SIZE_t start_pos) nogil:
    self.impurity_left = INFINITY
    self.impurity_right = INFINITY
    self.pos = start_pos
    self.feature = 0
    self.threshold = 0.
    self.improvement = -INFINITY
    self.proj_vec = NULL


cdef class BaseObliqueSplitter:
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
        self.proj_mat = NULL # max_features x n_features matrix. There are max_features vectors.

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

        cdef SIZE_t i
        if self.proj_mat:
            for i in range(self.max_features):
                free(self.proj_mat[i])
            # print("freed proj_mat vectors")
        free(self.proj_mat)
        # print("freed proj_mat")

    def __getstate__(self):
        return {}

    def __setstate__(self, d):
        pass

    cdef int init(self,
                   object X,
                   const DOUBLE_t[:, ::1] y,
                   DOUBLE_t* sample_weight,
                   np.ndarray X_idx_sorted=None) except -1:
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

        X_idx_sorted : ndarray, default=None
            The indexes of the sorted training input samples
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
        
        # TODO: throw memory error if this fails!
        # Reset projection matrix to 0
        cdef SIZE_t i
        for i in range(self.max_features):
            # I dont think we need to reallocate the memory. Just reset to 0.
            #safe_realloc(&self.proj_mat[i], self.n_features)
            memset(self.proj_mat[i], 0, self.n_features * sizeof(DTYPE_t))

    cdef int node_split(self, double impurity, ObliqueSplitRecord* split,
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

    cdef void sample_proj_mat(self, DTYPE_t** proj_mat) nogil:
        """ Sample the projection vector. 
        
        This is a placeholder method. 
        """

        pass


cdef class DenseObliqueSplitter(BaseObliqueSplitter):
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
                  DOUBLE_t* sample_weight,
                  np.ndarray X_idx_sorted=None) except -1:
        """Initialize the splitter

        Returns -1 in case of failure to allocate memory (and raise MemoryError)
        or 0 otherwise.
        """

        # Call parent init
        BaseObliqueSplitter.init(self, X, y, sample_weight)

        self.X = X

        cdef SIZE_t i
        cdef SIZE_t n_features = self.n_features
        cdef SIZE_t max_features = self.max_features

        self.proj_mat = <DTYPE_t**> malloc(max_features * sizeof(DTYPE_t*))
        for i in range(max_features):
            self.proj_mat[i] = <DTYPE_t*> malloc(n_features * sizeof(DTYPE_t))
            memset(self.proj_mat[i], 0, n_features * sizeof(DTYPE_t))


cdef class ObliqueSplitter(DenseObliqueSplitter):
    def __reduce__(self):
        """Enable pickling the splitter."""
        return (ObliqueSplitter, (self.criterion,
                               self.max_features,
                               self.min_samples_leaf,
                               self.min_weight_leaf,
                               self.feature_combinations,
                               self.random_state), self.__getstate__())

    cdef void sample_proj_mat(self, DTYPE_t** proj_mat) nogil:
        """
        SPORF Projection matrix.
        Randomly sample features to put in randomly sampled projection vectors
        weight = 1 or -1 with probability 0.5 
        """

        cdef SIZE_t n_features = self.n_features
        cdef SIZE_t max_features = self.max_features
        cdef SIZE_t n_non_zeros = self.n_non_zeros
        cdef UINT32_t* random_state = &self.rand_r_state

        cdef int i, feat_i, proj_i
        cdef DTYPE_t weight

        for i in range(0, n_non_zeros):

            proj_i = rand_int(0, max_features, random_state)
            feat_i = rand_int(0, n_features, random_state)
            weight = 1 if (rand_int(0, 2, random_state) == 1) else -1

            proj_mat[proj_i][feat_i] = weight

    cdef int node_split(self, double impurity, ObliqueSplitRecord* split,
                        SIZE_t* n_constant_features) nogil except -1:
        """Find the best split on node samples[start:end]
        Returns -1 in case of failure to allocate memory (and raise MemoryError)
        or 0 otherwise.
        """
        
        cdef SIZE_t n_samples = self.n_samples
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

        cdef SIZE_t f, i, j, p, partition_end
        cdef DTYPE_t temp_d

        # instantiate the projection matrix and a 
        # point for projection vectors to pass the selected projection vector
        cdef DTYPE_t** proj_mat = self.proj_mat
        cdef DTYPE_t* proj_vec

        # instantiate the split records
        _init_split(&best, end)

        # Sample the projection matrix
        self.sample_proj_mat(proj_mat)
        
        # For every vector in the projection matrix
        for f in range(max_features):
            proj_vec = proj_mat[f]
            current.feature = f
            current.proj_vec = proj_mat[f]

            # Compute linear combination of features
            memset(Xf + start, 0, (end - start) * sizeof(DTYPE_t))
            for i in range(start, end):
                for j in range(n_features):
                    Xf[i] += self.X[samples[i], j] * proj_vec[j]

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
                for j in range(n_features):
                    temp_d += self.X[samples[p], j] * best.proj_vec[j]

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

        # Skipping over constant features part cause its irrelevant

        # Return values
        split[0] = best
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

# XXX n/2 compile error - adding extra variable
cdef inline DTYPE_t median3(DTYPE_t* Xf, SIZE_t n) nogil:
    # Median of three pivot selection, after Bentley and McIlroy (1993).
    # Engineering a sort function. SP&E. Requires 8/3 comparisons on average.
    cdef SIZE_t mid = int(n/2)
    cdef DTYPE_t a = Xf[0], b = Xf[mid], c = Xf[n - 1]
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

# XXX Same problem here, casting to int
cdef void heapsort(DTYPE_t* Xf, SIZE_t* samples, SIZE_t n) nogil:
    cdef SIZE_t start, end

    # heapify
    start = int((n - 2) / 2)
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
