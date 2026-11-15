from libc.math cimport log2

from cython cimport floating

from sklearn.utils._typedefs cimport intp_t


cdef void simultaneous_sort(
    floating* values,
    intp_t* indices,
    intp_t n,
    bint use_three_way_partition=False,
) noexcept nogil:
    """Sort values and indices simultaneously by values.

    The numpy equivalent is:
        def simultaneous_sort(values, indices):
             i = np.argsort(values)
             return values[i], indices[i]

    Algorithm: Introsort (Musser, SP&E, 1997) with two variants for the
    quicksort part:

    - If use_three_way_partition is True, use 3-way partitioning:
      [x < pivot] [x == pivot] [x > pivot]. This variant is fast when
      working with many duplicate values, otherwise it's slower.
    - If use_three_way_partition is False, use 2-way partitioning:
      [x <= pivot] [pivot] [x >= pivot]. There are three parts too, but the middle
      part is only the selected pivot element, not all values equal to the pivot.

    Notes
    -----
    Arrays are manipulated via a pointer to their first element and their size
    to ease the processing of dynamically allocated buffers.

    TODO: In order to support discrete distance metrics, we need to have a
    simultaneous sort which breaks ties on indices when distances are
    identical. The best might be using a std::stable_sort and a Comparator
    which might need an Array of Structures (AoS) instead of the Structure of
    Arrays (SoA) currently used. An alternative would be to implement a stable
    sort ourselves, like the radix sort for instance.
    """
    if n == 0:
        return
    cdef intp_t maxd = 2 * <intp_t>log2(n)
    if use_three_way_partition:
        introsort_3way(values, indices, n, maxd)
    else:
        introsort_2way(values, indices, n, maxd)


def _py_simultaneous_sort(
    floating[::1] values,
    intp_t[::1] indices,
    intp_t n,
    *,
    bint use_three_way_partition,
):
    """Python wrapper used for testing."""
    simultaneous_sort(&values[0], &indices[0], n, use_three_way_partition)


cdef void introsort_2way(
    floating* values,
    intp_t* indices,
    intp_t n,
    intp_t maxd,
) noexcept nogil:
    cdef floating pivot
    cdef intp_t pivot_idx, i, j

    while n > 15:
        if maxd <= 0:   # max depth limit exceeded ("gone quadratic")
            heapsort(values, indices, n)
            return
        maxd -= 1

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

    # in the small-array case, insertion sort is faster
    insertion_sort(values, indices, n)


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

    while n > 15:
        if maxd <= 0:   # max depth limit exceeded ("gone quadratic")
            heapsort(values, indices, n)
            return
        maxd -= 1

        pivot = median3(values, n)

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

        # Three-way partition:
        # - values[:l] contains elements < pivot
        # - values[l:r] contains elements == pivot
        # - values[r:] contains elements > pivot

        # Recursively sort left side:
        introsort_3way(values, indices, l, maxd)

        # Continue with right side:
        values += r
        indices += r
        n -= r

    # in the small-array case, insertion sort is faster
    insertion_sort(values, indices, n)

# ------------ HEAP SORT -------------

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


# ------------ HELPERS -------------

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


cdef inline void swap(floating* values, intp_t* indices,
                      intp_t i, intp_t j) noexcept nogil:
    # Helper for sort
    values[i], values[j] = values[j], values[i]
    indices[i], indices[j] = indices[j], indices[i]


cdef inline void insertion_sort(
    floating* values, intp_t *indices, intp_t n
) noexcept nogil:
    cdef intp_t i, j, temp_idx
    cdef floating temp_val

    for i in range(1, n):
        temp_val = values[i]
        temp_idx = indices[i]

        j = i
        while j > 0 and values[j - 1] > temp_val:
            values[j] = values[j - 1]
            indices[j] = indices[j - 1]
            j -= 1

        values[j] = temp_val
        indices[j] = temp_idx
