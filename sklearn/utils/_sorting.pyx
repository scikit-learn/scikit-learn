cimport numpy as cnp

from cython cimport floating
from libc.math cimport log as ln

# TODO: In order to support discrete distance metrics, we need to have a
# simultaneous sort which breaks ties on indices when distances are identical.
# The best might be using a std::stable_sort and a Comparator which might need
# an Array of Structures (AoS) instead of the Structure of Arrays (SoA)
# currently used.

# Utilities functions

# TODO: Factor code also present in `tree._utils` or use `libc.math.log2` directly
cdef inline double log(double x) nogil:
    return ln(x) / ln(2.0)

cdef inline void _simultaneous_swap(
    floating* values,
    cnp.intp_t* indices,
    cnp.intp_t i,
    cnp.intp_t j,
) nogil:
    # Helper for sort
    values[i], values[j] = values[j], values[i]
    indices[i], indices[j] = indices[j], indices[i]

cdef inline floating _median3(
    floating* values,
    cnp.intp_t size,
) nogil:
    # Median of three pivot selection, after Bentley and McIlroy (1993).
    # Engineering a sort function. SP&E. Requires 8/3 comparisons on average.
    cdef floating a = values[0], b = values[size / 2], c = values[size - 1]
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

cdef inline void _sift_down(
    floating* values,
    cnp.intp_t* indices,
    cnp.intp_t start,
    cnp.intp_t end,
) nogil:
    # Restore heap order in values[start:end] by moving the max element to start.
    cdef cnp.intp_t child, maxind, root

    root = start
    while True:
        child = root * 2 + 1

        # find max of root, left child, right child
        maxind = root
        if child < end and values[maxind] < values[child]:
            maxind = child
        if child + 1 < end and values[maxind] < values[child + 1]:
            maxind = child + 1

        if maxind == root:
            break
        else:
            _simultaneous_swap(values, indices, root, maxind)
            root = maxind


# Sorting functions

cdef inline void simultaneous_introsort(
    floating* values,
    cnp.intp_t* indices,
    cnp.intp_t size,
) nogil:
    # Sort a Structure of Arrays pointed consisting of arrays of values and indices,
    # simultaneously, based on the values. Algorithm: Introsort (Musser, SP&E, 1997).
    if size == 0:
      return
    cdef int maxd = 2 * <int>log(size)
    _simultaneous_introsort(values, indices, size, maxd)


cdef int simultaneous_quick_sort(
    floating* values,
    cnp.intp_t* indices,
    cnp.intp_t size,
) nogil:
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
    cdef:
        cnp.intp_t pivot_idx, i, store_idx
        floating pivot_val

    # in the small-array case, do things efficiently
    if size <= 1:
        pass
    elif size == 2:
        if values[0] > values[1]:
            _simultaneous_swap(values, indices, 0, 1)
    elif size == 3:
        if values[0] > values[1]:
            _simultaneous_swap(values, indices, 0, 1)
        if values[1] > values[2]:
            _simultaneous_swap(values, indices, 1, 2)
            if values[0] > values[1]:
                _simultaneous_swap(values, indices, 0, 1)
    else:
        # Determine the pivot using the median-of-three rule.
        # The smallest of the three is moved to the beginning of the array,
        # the middle (the pivot value) is moved to the end, and the largest
        # is moved to the pivot index.
        pivot_idx = size // 2
        if values[0] > values[size - 1]:
            _simultaneous_swap(values, indices, 0, size - 1)
        if values[size - 1] > values[pivot_idx]:
            _simultaneous_swap(values, indices, size - 1, pivot_idx)
            if values[0] > values[size - 1]:
                _simultaneous_swap(values, indices, 0, size - 1)
        pivot_val = values[size - 1]

        # Partition indices about pivot.  At the end of this operation,
        # pivot_idx will contain the pivot value, everything to the left
        # will be smaller, and everything to the right will be larger.
        store_idx = 0
        for i in range(size - 1):
            if values[i] < pivot_val:
                _simultaneous_swap(values, indices, i, store_idx)
                store_idx += 1
        _simultaneous_swap(values, indices, store_idx, size - 1)
        pivot_idx = store_idx

        # Recursively sort each side of the pivot
        if pivot_idx > 1:
            simultaneous_quick_sort(values, indices, pivot_idx)
        if pivot_idx + 2 < size:
            simultaneous_quick_sort(values + pivot_idx + 1,
                                    indices + pivot_idx + 1,
                                    size - pivot_idx - 1)
    return 0

# Introsort with median of 3 pivot selection and 3-way partition function
# (robust to repeated elements, e.g. lots of zero features).
cdef void _simultaneous_introsort(
    floating* values,
    cnp.intp_t* indices,
    cnp.intp_t size,
    int maxd,
) nogil:
    cdef floating pivot
    cdef cnp.intp_t i, l, r

    while size > 1:
        if maxd <= 0:   # max depth limit exceeded ("gone quadratic")
            simultaneous_heapsort(values, indices, size)
            return
        maxd -= 1

        pivot = _median3(values, size)

        # Three-way partition.
        i = l = 0
        r = size
        while i < r:
            if values[i] < pivot:
                _simultaneous_swap(values, indices, i, l)
                i += 1
                l += 1
            elif values[i] > pivot:
                r -= 1
                _simultaneous_swap(values, indices, i, r)
            else:
                i += 1

        _simultaneous_introsort(values, indices, l, maxd)
        values += r
        indices += r
        size -= r

cdef void simultaneous_heapsort(
    floating* values,
    cnp.intp_t* indices,
    cnp.intp_t size,
) nogil:
    cdef cnp.intp_t start, end

    # heapify
    start = (size - 2) / 2
    end = size
    while True:
        _sift_down(values, indices, start, end)
        if start == 0:
            break
        start -= 1

    # sort by shrinking the heap, putting the max element immediately after it
    end = size - 1
    while end > 0:
        _simultaneous_swap(values, indices, 0, end)
        _sift_down(values, indices, 0, end)
        end = end - 1
