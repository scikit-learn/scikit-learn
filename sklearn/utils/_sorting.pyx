from cython cimport floating

cdef inline void dual_swap(
    floating* darr,
    intp_t *iarr,
    intp_t a,
    intp_t b,
) noexcept nogil:
    """Swap the values at index a and b of both darr and iarr"""
    cdef floating dtmp = darr[a]
    darr[a] = darr[b]
    darr[b] = dtmp

    cdef intp_t itmp = iarr[a]
    iarr[a] = iarr[b]
    iarr[b] = itmp


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
        intp_t pivot_idx, i, j
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
        i = 1  # the median-of-three step ensures values[0] <= pivot
        j = size - 2  # the median-of-three step ensures values[-1] >= pivot
        while True:
            # Find element >= pivot from left
            while i <= j and values[i] < pivot_val:
                i += 1
            # Find element <= pivot from right
            while i <= j and values[j] > pivot_val:
                j -= 1
            if i >= j:
                break
            dual_swap(values, indices, i, j)
            i += 1
            j -= 1

        # Put pivot at pivot_idx
        pivot_idx = i
        dual_swap(values, indices, pivot_idx, size - 1)

        # Recursively sort each side of the pivot
        if pivot_idx > 1:
            simultaneous_sort(values, indices, pivot_idx)
        if pivot_idx + 2 < size:
            simultaneous_sort(values + pivot_idx + 1,
                              indices + pivot_idx + 1,
                              size - pivot_idx - 1)
    return 0


def _py_simultaneous_sort(
    floating[::1] values,
    intp_t[::1] indices,
    intp_t size,
):
    """Python wrapper used for testing."""
    simultaneous_sort(&values[0], &indices[0], size)
