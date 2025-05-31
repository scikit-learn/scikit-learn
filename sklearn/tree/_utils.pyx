# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from libc.math cimport isnan
from libc.math cimport log as ln
from libc.stdlib cimport free, realloc

import numpy as np

cimport numpy as cnp

cnp.import_array()

from ..utils._random cimport our_rand_r

# =============================================================================
# Helper functions
# =============================================================================

cdef int safe_realloc(realloc_ptr* p, size_t nelems) except -1 nogil:
    # sizeof(realloc_ptr[0]) would be more like idiomatic C, but causes Cython
    # 0.20.1 to crash.
    cdef size_t nbytes = nelems * sizeof(p[0][0])
    if nbytes / sizeof(p[0][0]) != nelems:
        # Overflow in the multiplication
        raise MemoryError(f"could not allocate ({nelems} * {sizeof(p[0][0])}) bytes")

    cdef realloc_ptr tmp = <realloc_ptr>realloc(p[0], nbytes)
    if tmp == NULL:
        raise MemoryError(f"could not allocate {nbytes} bytes")

    p[0] = tmp
    return 0


def _realloc_test():
    # Helper for tests. Tries to allocate <size_t>(-1) / 2 * sizeof(size_t)
    # bytes, which will always overflow.
    cdef intp_t* p = NULL
    safe_realloc(&p, <size_t>(-1) / 2)
    if p != NULL:
        free(p)
        assert False


cdef inline cnp.ndarray sizet_ptr_to_ndarray(intp_t* data, intp_t size):
    """Return copied data as 1D numpy array of intp's."""
    cdef cnp.npy_intp shape[1]
    shape[0] = <cnp.npy_intp> size
    return cnp.PyArray_SimpleNewFromData(1, shape, cnp.NPY_INTP, data).copy()


cdef inline cnp.ndarray int32t_ptr_to_ndarray(int32_t* data, intp_t size):
    """Encapsulate data into a 1D numpy array of int32's."""
    cdef cnp.npy_intp shape[1]
    shape[0] = <cnp.npy_intp> size
    return cnp.PyArray_SimpleNewFromData(1, shape, cnp.NPY_INT32, data)

cdef inline intp_t rand_int(intp_t low, intp_t high,
                            uint32_t* random_state) noexcept nogil:
    """Generate a random integer in [low; end)."""
    return low + our_rand_r(random_state) % (high - low)


cdef inline float64_t rand_uniform(float64_t low, float64_t high,
                                   uint32_t* random_state) noexcept nogil:
    """Generate a random float64_t in [low; high)."""
    return ((high - low) * <float64_t> our_rand_r(random_state) /
            <float64_t> RAND_R_MAX) + low


cdef inline float64_t log(float64_t x) noexcept nogil:
    return ln(x) / ln(2.0)

# =============================================================================
# WeightedPQueue data structure
# =============================================================================

cdef class WeightedPQueue:
    """A priority queue class, always sorted in increasing order.

    Attributes
    ----------
    capacity : intp_t
        The capacity of the priority queue.

    array_ptr : intp_t
        The water mark of the priority queue; the priority queue grows from
        left to right in the array ``array_``. ``array_ptr`` is always
        less than ``capacity``.

    array_ : WeightedPQueueRecord*
        The array of priority queue records. The minimum element is on the
        left at index 0, and the maximum element is on the right at index
        ``array_ptr-1``.
    """

    def __cinit__(self, intp_t capacity):
        self.capacity = capacity
        self.array_ptr = 0
        safe_realloc(&self.array_, capacity)

    def __dealloc__(self):
        free(self.array_)

    cdef int reset(self) except -1 nogil:
        """Reset the WeightedPQueue to its state at construction

        Return -1 in case of failure to allocate memory (and raise MemoryError)
        or 0 otherwise.
        """
        self.array_ptr = 0
        # Since safe_realloc can raise MemoryError, use `except -1`
        safe_realloc(&self.array_, self.capacity)
        return 0

    cdef bint is_empty(self) noexcept nogil:
        return self.array_ptr <= 0

    cdef intp_t size(self) noexcept nogil:
        return self.array_ptr

    cdef int push(self, float64_t data, float64_t weight) except -1 nogil:
        """Push record on the array.

        Return -1 in case of failure to allocate memory (and raise MemoryError)
        or 0 otherwise.
        """
        cdef intp_t array_ptr = self.array_ptr
        cdef WeightedPQueueRecord* array = NULL
        cdef intp_t i

        # Resize if capacity not sufficient
        if array_ptr >= self.capacity:
            self.capacity *= 2
            # Since safe_realloc can raise MemoryError, use `except -1`
            safe_realloc(&self.array_, self.capacity)

        # Put element as last element of array
        array = self.array_
        array[array_ptr].data = data
        array[array_ptr].weight = weight

        # bubble last element up according until it is sorted
        # in ascending order
        i = array_ptr
        while(i != 0 and array[i].data < array[i-1].data):
            array[i], array[i-1] = array[i-1], array[i]
            i -= 1

        # Increase element count
        self.array_ptr = array_ptr + 1
        return 0

    cdef int remove(self, float64_t data, float64_t weight) noexcept nogil:
        """Remove a specific value/weight record from the array.
        Returns 0 if successful, -1 if record not found."""
        cdef intp_t array_ptr = self.array_ptr
        cdef WeightedPQueueRecord* array = self.array_
        cdef intp_t idx_to_remove = -1
        cdef intp_t i

        if array_ptr <= 0:
            return -1

        # find element to remove
        for i in range(array_ptr):
            if array[i].data == data and array[i].weight == weight:
                idx_to_remove = i
                break

        if idx_to_remove == -1:
            return -1

        # shift the elements after the removed element
        # to the left.
        for i in range(idx_to_remove, array_ptr-1):
            array[i] = array[i+1]

        self.array_ptr = array_ptr - 1
        return 0

    cdef int pop(self, float64_t* data, float64_t* weight) noexcept nogil:
        """Remove the top (minimum) element from array.
        Returns 0 if successful, -1 if nothing to remove."""
        cdef intp_t array_ptr = self.array_ptr
        cdef WeightedPQueueRecord* array = self.array_
        cdef intp_t i

        if array_ptr <= 0:
            return -1

        data[0] = array[0].data
        weight[0] = array[0].weight

        # shift the elements after the removed element
        # to the left.
        for i in range(0, array_ptr-1):
            array[i] = array[i+1]

        self.array_ptr = array_ptr - 1
        return 0

    cdef int peek(self, float64_t* data, float64_t* weight) noexcept nogil:
        """Write the top element from array to a pointer.
        Returns 0 if successful, -1 if nothing to write."""
        cdef WeightedPQueueRecord* array = self.array_
        if self.array_ptr <= 0:
            return -1
        # Take first value
        data[0] = array[0].data
        weight[0] = array[0].weight
        return 0

    cdef float64_t get_weight_from_index(self, intp_t index) noexcept nogil:
        """Given an index between [0,self.current_capacity], access
        the appropriate heap and return the requested weight"""
        cdef WeightedPQueueRecord* array = self.array_

        # get weight at index
        return array[index].weight

    cdef float64_t get_value_from_index(self, intp_t index) noexcept nogil:
        """Given an index between [0,self.current_capacity], access
        the appropriate heap and return the requested value"""
        cdef WeightedPQueueRecord* array = self.array_

        # get value at index
        return array[index].data

# =============================================================================
# WeightedMedianCalculator data structure
# =============================================================================

cdef class WeightedMedianCalculator:
    """A class to handle calculation of the weighted median from streams of
    data. To do so, it maintains a parameter ``k`` such that the sum of the
    weights in the range [0,k) is greater than or equal to half of the total
    weight. By minimizing the value of ``k`` that fulfills this constraint,
    calculating the median is done by either taking the value of the sample
    at index ``k-1`` of ``samples`` (samples[k-1].data) or the average of
    the samples at index ``k-1`` and ``k`` of ``samples``
    ((samples[k-1] + samples[k]) / 2).

    Attributes
    ----------
    initial_capacity : intp_t
        The initial capacity of the WeightedMedianCalculator.

    samples : WeightedPQueue
        Holds the samples (consisting of values and their weights) used in the
        weighted median calculation.

    total_weight : float64_t
        The sum of the weights of items in ``samples``. Represents the total
        weight of all samples used in the median calculation.

    k : intp_t
        Index used to calculate the median.

    sum_w_0_k : float64_t
        The sum of the weights from samples[0:k]. Used in the weighted
        median calculation; minimizing the value of ``k`` such that
        ``sum_w_0_k`` >= ``total_weight / 2`` provides a mechanism for
        calculating the median in constant time.

    """

    def __cinit__(self, intp_t initial_capacity):
        self.initial_capacity = initial_capacity
        self.samples = WeightedPQueue(initial_capacity)
        self.total_weight = 0
        self.k = 0
        self.sum_w_0_k = 0

    cdef intp_t size(self) noexcept nogil:
        """Return the number of samples in the
        WeightedMedianCalculator"""
        return self.samples.size()

    cdef int reset(self) except -1 nogil:
        """Reset the WeightedMedianCalculator to its state at construction

        Return -1 in case of failure to allocate memory (and raise MemoryError)
        or 0 otherwise.
        """
        # samples.reset (WeightedPQueue.reset) uses safe_realloc, hence
        # except -1
        self.samples.reset()
        self.total_weight = 0
        self.k = 0
        self.sum_w_0_k = 0
        return 0

    cdef int push(self, float64_t data, float64_t weight) except -1 nogil:
        """Push a value and its associated weight to the WeightedMedianCalculator

        Return -1 in case of failure to allocate memory (and raise MemoryError)
        or 0 otherwise.
        """
        cdef int return_value
        cdef float64_t original_median = 0.0

        if self.size() != 0:
            original_median = self.get_median()
        # samples.push (WeightedPQueue.push) uses safe_realloc, hence except -1
        return_value = self.samples.push(data, weight)
        self.update_median_parameters_post_push(data, weight,
                                                original_median)
        return return_value

    cdef int update_median_parameters_post_push(
            self, float64_t data, float64_t weight,
            float64_t original_median) noexcept nogil:
        """Update the parameters used in the median calculation,
        namely `k` and `sum_w_0_k` after an insertion"""

        # trivial case of one element.
        if self.size() == 1:
            self.k = 1
            self.total_weight = weight
            self.sum_w_0_k = self.total_weight
            return 0

        # get the original weighted median
        self.total_weight += weight

        if data < original_median:
            # inserting below the median, so increment k and
            # then update self.sum_w_0_k accordingly by adding
            # the weight that was added.
            self.k += 1
            # update sum_w_0_k by adding the weight added
            self.sum_w_0_k += weight

            # minimize k such that sum(W[0:k]) >= total_weight / 2
            # minimum value of k is 1
            while(self.k > 1 and ((self.sum_w_0_k -
                                   self.samples.get_weight_from_index(self.k-1))
                                  >= self.total_weight / 2.0)):
                self.k -= 1
                self.sum_w_0_k -= self.samples.get_weight_from_index(self.k)
            return 0

        if data >= original_median:
            # inserting above or at the median
            # minimize k such that sum(W[0:k]) >= total_weight / 2
            while(self.k < self.samples.size() and
                  (self.sum_w_0_k < self.total_weight / 2.0)):
                self.k += 1
                self.sum_w_0_k += self.samples.get_weight_from_index(self.k-1)
            return 0

    cdef int remove(self, float64_t data, float64_t weight) noexcept nogil:
        """Remove a value from the MedianHeap, removing it
        from consideration in the median calculation
        """
        cdef int return_value
        cdef float64_t original_median = 0.0

        if self.size() != 0:
            original_median = self.get_median()

        return_value = self.samples.remove(data, weight)
        self.update_median_parameters_post_remove(data, weight,
                                                  original_median)
        return return_value

    cdef int pop(self, float64_t* data, float64_t* weight) noexcept nogil:
        """Pop a value from the MedianHeap, starting from the
        left and moving to the right.
        """
        cdef int return_value
        cdef float64_t original_median = 0.0

        if self.size() != 0:
            original_median = self.get_median()

        # no elements to pop
        if self.samples.size() == 0:
            return -1

        return_value = self.samples.pop(data, weight)
        self.update_median_parameters_post_remove(data[0],
                                                  weight[0],
                                                  original_median)
        return return_value

    cdef int update_median_parameters_post_remove(
            self, float64_t data, float64_t weight,
            float64_t original_median) noexcept nogil:
        """Update the parameters used in the median calculation,
        namely `k` and `sum_w_0_k` after a removal"""
        # reset parameters because it there are no elements
        if self.samples.size() == 0:
            self.k = 0
            self.total_weight = 0
            self.sum_w_0_k = 0
            return 0

        # trivial case of one element.
        if self.samples.size() == 1:
            self.k = 1
            self.total_weight -= weight
            self.sum_w_0_k = self.total_weight
            return 0

        # get the current weighted median
        self.total_weight -= weight

        if data < original_median:
            # removing below the median, so decrement k and
            # then update self.sum_w_0_k accordingly by subtracting
            # the removed weight

            self.k -= 1
            # update sum_w_0_k by removing the weight at index k
            self.sum_w_0_k -= weight

            # minimize k such that sum(W[0:k]) >= total_weight / 2
            # by incrementing k and updating sum_w_0_k accordingly
            # until the condition is met.
            while(self.k < self.samples.size() and
                  (self.sum_w_0_k < self.total_weight / 2.0)):
                self.k += 1
                self.sum_w_0_k += self.samples.get_weight_from_index(self.k-1)
            return 0

        if data >= original_median:
            # removing above the median
            # minimize k such that sum(W[0:k]) >= total_weight / 2
            while(self.k > 1 and ((self.sum_w_0_k -
                                   self.samples.get_weight_from_index(self.k-1))
                                  >= self.total_weight / 2.0)):
                self.k -= 1
                self.sum_w_0_k -= self.samples.get_weight_from_index(self.k)
            return 0

    cdef float64_t get_median(self) noexcept nogil:
        """Write the median to a pointer, taking into account
        sample weights."""
        if self.sum_w_0_k == (self.total_weight / 2.0):
            # split median
            return (self.samples.get_value_from_index(self.k) +
                    self.samples.get_value_from_index(self.k-1)) / 2.0
        if self.sum_w_0_k > (self.total_weight / 2.0):
            # whole median
            return self.samples.get_value_from_index(self.k-1)


def _any_isnan_axis0(const float32_t[:, :] X):
    """Same as np.any(np.isnan(X), axis=0)"""
    cdef:
        intp_t i, j
        intp_t n_samples = X.shape[0]
        intp_t n_features = X.shape[1]
        uint8_t[::1] isnan_out = np.zeros(X.shape[1], dtype=np.bool_)

    with nogil:
        for i in range(n_samples):
            for j in range(n_features):
                if isnan_out[j]:
                    continue
                if isnan(X[i, j]):
                    isnan_out[j] = True
                    break
    return np.asarray(isnan_out)


cdef inline void copy_memview_to_array(
    BITSET_INNER_DTYPE_C[:] memview,
    BITSET_DTYPE_C arr
) noexcept nogil:
    for i in range(8):
        arr[i] = memview[i]

cdef inline void copy_array_to_memview(
    BITSET_INNER_DTYPE_C[:] memview,
    BITSET_DTYPE_C arr
) noexcept nogil:
    for i in range(8):
        memview[i] = arr[i]


cdef inline bint goes_left(
    float32_t feature_value,
    SplitValue split,
    int32_t n_categories,
    BITSET_INNER_DTYPE_C[:] cat_bitsets
) noexcept nogil:
    """Determine whether a sample goes to the left or right child node.

    For numerical features, ``(-inf, split.threshold]`` is the left child, and
    ``(split.threshold, inf)`` the right child.
    For categorical features, if the corresponding bit for the category is set
    in cachebits, the left child isused, and if not set, the right child. If
    the given input category is larger than the ``n_categories``, the right
    child is assumed.

    Attributes
    ----------
    feature_value : float32_t
        The value of the feature for which the decision needs to be made.
    split : SplitValue
        The union (of float64_t and BITSET_INNER_DTYPE_C) indicating the split. However, it
        is used (as a float64_t) only for numerical features.
    n_categories : int32_t
        The number of categories present in the feature in question. The
        feature is considered a numerical one and not a categorical one if
        n_categories is negative.
    cat_bitsets : BITSET_INNER_DTYPE_C[:]
        The array containing the expansion of split.cat_split. The function
        setup_cat_cache is the one filling it.

    Returns
    -------
    result : bint
        Indicating whether the left branch should be used.
    """
    cdef intp_t idx

    if n_categories < 0:
        # Non-categorical feature
        return feature_value <= split.threshold
    else:
        # Categorical feature, using bit cache
        if (<uint8_t> feature_value) < n_categories:
            # idx = (<intp_t> feature_value) // 64
            # offset = (<intp_t> feature_value) % 64
            # return bs_get(cat_cache[idx], offset)
            # TODO: allow uint16_t
            return in_bitset_memoryview(cat_bitsets, <uint8_t> feature_value)
        else:
            return 0


# cdef inline void setup_cat_cache_memoryview(
#     BITSET_INNER_DTYPE_C[:] cachebits,
#     BITSET_DTYPE_C cat_split,
#     int32_t n_categories,
# ) noexcept nogil:
#     cdef int32_t j
#     cdef uint32_t rng_seed, rand_set_bit
#     cdef intp_t cache_size = (n_categories + 31) // 32
#     if n_categories > 0:
#         if cat_split & 1:
#             # random categorical split
#             for j in range(cache_size):
#                 cachebits[j] = 0

#             # recover the random seed by shifting right by 32 bits
#             rng_seed = cat_split >> 32
#             for j in range(n_categories):
#                 rand_set_bit = rand_int(0, 2, &rng_seed)
#                 if not rand_set_bit:
#                     continue
#                 set_bitset_memoryview(cachebits, j)
#         else:
#             # best split, so just copy back over the split
#             copy_array_to_memview(cachebits, cat_split)

# cdef inline void setup_cat_cache(
#     BITSET_t[:] cachebits,
#     BITSET_t cat_split,
#     int32_t n_categories
# ) noexcept nogil:
#     """Populate the bits of the category cache from a split.
#     Attributes
#     ----------
#     cachebits : BITSET_t[::1]
#         This is a pointer to the output array. The size of the array should be
#         ``ceil(n_categories / 64)``. This function assumes the required
#         memory is allocated for the array by the caller.
#     cat_split : BITSET_t
#         If ``least significant bit == 0``:
#             It stores the split of the maximum 64 categories in its bits.
#             This is used in `BestSplitter`, and without loss of generality it
#             is assumed to be even, i.e. for any odd value there is an
#             equivalent even ``cat_split``.
#         If ``least significant bit == 1``:
#             It is a random split, and the 32 most significant bits of
#             ``cat_split`` contain the random seed of the split. The
#             ``n_categories`` lowest bits of ``cachebits`` are then filled with
#             random zeros and ones given the random seed.
#     n_categories : int32_t
#         The number of categories.
#     """
#     cdef int32_t j
#     cdef uint32_t rng_seed, val
#     cdef intp_t cache_size = (n_categories + 63) // 64
#     if n_categories > 0:
#         if cat_split & 1:
#             # RandomSplitter
#             for j in range(cache_size):
#                 cachebits[j] = 0
#             rng_seed = cat_split >> 32
#             for j in range(n_categories):
#                 val = rand_int(0, 2, &rng_seed)
#                 if not val:
#                     continue
#                 cachebits[j // 64] = bs_set(cachebits[j // 64], j % 64)
#         else:
#             # BestSplitter
#             # In practice, cache_size here should ALWAYS be 1
#             # XXX TODO: check cache_size == 1?
#             cachebits[0] = cat_split
