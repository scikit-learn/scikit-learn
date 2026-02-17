from cython cimport floating
from libc.math cimport log2


# This file contains 2 Cython implementation of:
#   def simultaneous_sort(dist, idx):
#       i = np.argsort(dist)
#       return dist[i], idx[i]

# Algorithm: Introsort (Musser, SP&E, 1997) with two variants for the quicksort part:
# - "2-way" partitioning: at each step, partition the current array in... 3 parts:
#   [x <= pivot] [pivot] [x >= pivot] (the middle part is only one element)
# - 3-way partitioning: at each step, partition the current array in 3 parts:
#   [x < pivot] [x == pivot] [x > pivot]. This variant is faster when working
#   with many duplicate values

# Note:
# Arrays are manipulated via a pointer to there first element and their size
# as to ease the processing of dynamically allocated buffers.

# TODO: In order to support discrete distance metrics, we need to have a
# simultaneous sort which breaks ties on indices when distances are identical.
# The best might be using a std::stable_sort and a Comparator which might need
# an Array of Structures (AoS) instead of the Structure of Arrays (SoA)
# currently used.
# An alternative would be to implement a stable sort ourselves, like the
# radix sort for instance


cdef void simultaneous_sort(floating* dist, intp_t* idx, intp_t size) noexcept nogil:
    if size == 0:
        return
    cdef intp_t maxd = 2 * <intp_t>log2(size)
    introsort_2way(dist, idx, size, maxd)


cdef void sort(floating* values, intp_t* indices, intp_t n) noexcept nogil:
    if n == 0:
        return
    cdef intp_t maxd = 2 * <intp_t>log2(n)
    introsort_3way(values, indices, n, maxd)


def _py_sort(
    floating[::1] values,
    intp_t[::1] indices,
    intp_t n,
    bint use_three_way_partition=False,
):
    """Python wrapper used for testing."""
    if use_three_way_partition:
        sort(&values[0], &indices[0], n)
    else:
        simultaneous_sort(&values[0], &indices[0], n)


cdef void introsort_2way(
    floating* values,
    intp_t* indices,
    intp_t n,
    intp_t maxd,
) noexcept nogil:
    cdef floating pivot
    cdef intp_t pivot_idx, i, j

    # in the small-array case, do things efficiently
    while n > 1:
        if maxd <= 0:   # max depth limit exceeded ("gone quadratic")
            heapsort(values, indices, n)
            return
        maxd -= 1
        # Determine the pivot using the median-of-three rule:
        pivot = inplace_median3(values, indices, n)

        i = 1  # the median3 step ensures values[0] <= pivot
        j = n - 2  # the median3 step ensures values[-1] >= pivot
        while True:
            # Find element >= pivot from left
            while i <= j and values[i] < pivot:
                i += 1
            # Find element <= pivot from right
            while i <= j and values[j] > pivot:
                j -= 1
            if i >= j:
                break
            swap(values, indices, i, j)
            i += 1
            j -= 1

        # Put pivot at pivot_idx
        pivot_idx = i
        swap(values, indices, pivot_idx, n - 1)

        # Recursively sort left side of the pivot
        introsort_2way(values, indices, pivot_idx, maxd)

        # Continue with right side:
        values += pivot_idx + 1
        indices += pivot_idx + 1
        n -= pivot_idx + 1


cdef void introsort_3way(
    floating* values, intp_t *indices,
    intp_t n, intp_t maxd
) noexcept nogil:
    """
    Introsort with median of 3 pivot selection and 3-way partition function
    (fast for repeated elements, e.g. lots of zeros).
    """
    cdef floating pivot
    cdef intp_t i, l, r

    while n > 1:
        if maxd <= 0:   # max depth limit exceeded ("gone quadratic")
            heapsort(values, indices, n)
            return
        maxd -= 1

        pivot = median3(values, n)

        # Three-way partition.
        i = l = 0
        r = n
        while i < r:
            if values[i] < pivot:
                swap(values, indices, i, l)
                i += 1
                l += 1
            elif values[i] > pivot:
                r -= 1
                swap(values, indices, i, r)
            else:
                i += 1

        introsort_3way(values, indices, l, maxd)
        values += r
        indices += r
        n -= r


cdef inline void sift_down(floating* feature_values, intp_t* samples,
                           intp_t start, intp_t end) noexcept nogil:
    # Restore heap order in feature_values[start:end] by moving the max element to start.
    cdef intp_t child, maxind, root

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


cdef void heapsort(floating* feature_values, intp_t* samples, intp_t n) noexcept nogil:
    cdef intp_t start, end

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


cdef inline void swap(floating* values, intp_t* indices,
                      intp_t i, intp_t j) noexcept nogil:
    # Helper for sort
    values[i], values[j] = values[j], values[i]
    indices[i], indices[j] = indices[j], indices[i]


cdef inline void insertion_sort(
    floating* values, intp_t *indices, intp_t n
) noexcept nogil:
    # Note: this sort is stable:
    cdef intp_t i = 1
    cdef intp_t j, temp_idx
    cdef floating temp_val

    for i in range(n):
        # Store current element to be inserted in sorted position
        temp_val = values[i]
        temp_idx = indices[i]

        # Find correct position by shifting larger elements right
        j = i  # Start comparing from current position
        # Shift elements that are larger than temp_val to the right
        # Continue while we haven't reached the start and
        # the previous element is larger than our value to insert
        while j > 0 and values[j - 1] > temp_val:
            values[j] = values[j - 1]        # Shift value right
            indices[j] = indices[j - 1]  # Shift index right
            j -= 1  # Move comparison position left

        # Insert the stored element at its correct sorted position
        values[j] = temp_val
        indices[j] = temp_idx


cdef inline floating inplace_median3(floating* values, intp_t* indices, intp_t n) noexcept nogil:
    # # Median of three pivot selection
    # The smallest of the three is moved to the beginning of the array,
    # the middle (the pivot value) is moved to the end, and the largest
    # is moved to the pivot index.
    pivot_idx = n // 2
    if values[0] > values[n - 1]:
        swap(values, indices, 0, n - 1)
    if values[n - 1] > values[pivot_idx]:
        swap(values, indices, n - 1, pivot_idx)
        if values[0] > values[n - 1]:
            swap(values, indices, 0, n - 1)
    return values[n - 1]


cdef inline floating median3(floating* feature_values, intp_t n) noexcept nogil:
    # Median of three pivot selection, after Bentley and McIlroy (1993).
    # Engineering a sort function. SP&E. Requires 8/3 comparisons on average.
    cdef floating a = feature_values[0], b = feature_values[n / 2], c = feature_values[n - 1]
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
