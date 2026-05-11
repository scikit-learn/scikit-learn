from libc.math cimport log2

from cython cimport floating

from sklearn.utils._typedefs cimport float32_t, intp_t


cdef inline void dual_swap(
    floating* darr,
    intp_t *iarr,
    intp_t a,
    intp_t b,
) noexcept nogil:
    """Swap the values at index a and b of both darr and iarr"""
    darr[a], darr[b] = darr[b], darr[a]
    iarr[a], iarr[b] = iarr[b], iarr[a]


cdef int simultaneous_sort(
    floating* values,
    intp_t* indices,
    intp_t size,
) noexcept nogil:
    """
    Perform a recursive quicksort on the values array as to sort them ascendingly.
    This simultaneously performs the swaps on both the values and the indices arrays.

    The numpy equivalent is:

        def simultaneous_sort(dist, idx):
             i = np.argsort(dist)
             return dist[i], idx[i]

    Notes
    -----
    Arrays are manipulated via a pointer to there first element and their size
    as to ease the processing of dynamically allocated buffers.
    """
    # TODO: In order to support discrete distance metrics, we need to have a
    # simultaneous sort which breaks ties on indices when distances are identical.
    # The best might be using a std::stable_sort and a Comparator which might need
    # an Array of Structures (AoS) instead of the Structure of Arrays (SoA)
    # currently used.
    cdef:
        intp_t pivot_idx, i, store_idx
        floating pivot_val

    # in the small-array case, do things efficiently
    if size <= 1:
        pass
    elif size == 2:
        if values[0] > values[1]:
            dual_swap(values, indices, 0, 1)
    elif size == 3:
        if values[0] > values[1]:
            dual_swap(values, indices, 0, 1)
        if values[1] > values[2]:
            dual_swap(values, indices, 1, 2)
            if values[0] > values[1]:
                dual_swap(values, indices, 0, 1)
    else:
        # Determine the pivot using the median-of-three rule.
        # The smallest of the three is moved to the beginning of the array,
        # the middle (the pivot value) is moved to the end, and the largest
        # is moved to the pivot index.
        pivot_idx = size // 2
        if values[0] > values[size - 1]:
            dual_swap(values, indices, 0, size - 1)
        if values[size - 1] > values[pivot_idx]:
            dual_swap(values, indices, size - 1, pivot_idx)
            if values[0] > values[size - 1]:
                dual_swap(values, indices, 0, size - 1)
        pivot_val = values[size - 1]

        # Partition indices about pivot.  At the end of this operation,
        # pivot_idx will contain the pivot value, everything to the left
        # will be smaller, and everything to the right will be larger.
        store_idx = 0
        for i in range(size - 1):
            if values[i] < pivot_val:
                dual_swap(values, indices, i, store_idx)
                store_idx += 1
        dual_swap(values, indices, store_idx, size - 1)
        pivot_idx = store_idx

        # Recursively sort each side of the pivot
        if pivot_idx > 1:
            simultaneous_sort(values, indices, pivot_idx)
        if pivot_idx + 2 < size:
            simultaneous_sort(values + pivot_idx + 1,
                              indices + pivot_idx + 1,
                              size - pivot_idx - 1)
    return 0


def _py_sort(float32_t[::1] feature_values, intp_t[::1] samples, intp_t n):
    """Used for testing sort."""
    sort(&feature_values[0], &samples[0], n)


# Sort n-element arrays pointed to by feature_values and samples, simultaneously,
# by the values in feature_values. Algorithm: Introsort (Musser, SP&E, 1997).
cdef void sort(floating* feature_values, intp_t* samples, intp_t n) noexcept nogil:
    if n == 0:
        return
    cdef intp_t maxd = 2 * <intp_t>log2(n)
    introsort(feature_values, samples, n, maxd)


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


# Introsort with median of 3 pivot selection and 3-way partition function
# (robust to repeated elements, e.g. lots of zero features).
cdef void introsort(floating* feature_values, intp_t *samples,
                    intp_t n, intp_t maxd) noexcept nogil:
    cdef floating pivot
    cdef intp_t i, l, r

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
                dual_swap(feature_values, samples, i, l)
                i += 1
                l += 1
            elif feature_values[i] > pivot:
                r -= 1
                dual_swap(feature_values, samples, i, r)
            else:
                i += 1

        introsort(feature_values, samples, l, maxd)
        feature_values += r
        samples += r
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
            dual_swap(feature_values, samples, root, maxind)
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
        dual_swap(feature_values, samples, 0, end)
        sift_down(feature_values, samples, 0, end)
        end = end - 1
