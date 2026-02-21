from libc.math cimport log2

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


cdef inline floating median3(floating* values, intp_t n) noexcept nogil:
    # Median of three pivot selection, after Bentley and McIlroy (1993).
    # Engineering a sort function. SP&E. Requires 8/3 comparisons on average.
    cdef floating a = values[0], b = values[n // 2], c = values[n - 1]
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


cdef inline void sift_down(
    floating* values,
    intp_t* indices,
    intp_t start,
    intp_t end,
) noexcept nogil:
    # Restore heap order in values[start:end] by sifting down the element at start.
    cdef intp_t child, maxind, root

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
            dual_swap(values, indices, root, maxind)
            root = maxind


cdef void heapsort(
    floating* values,
    intp_t* indices,
    intp_t n,
) noexcept nogil:
    cdef intp_t start, end

    # heapify
    start = (n - 2) // 2
    end = n
    while True:
        sift_down(values, indices, start, end)
        if start == 0:
            break
        start -= 1

    # sort by shrinking the heap, putting the max element immediately after it
    end = n - 1
    while end > 0:
        dual_swap(values, indices, 0, end)
        sift_down(values, indices, 0, end)
        end -= 1


# Introsort with median-of-3 pivot selection and 3-way partition function
# (robust to repeated elements, e.g. lots of duplicate distances).
cdef void introsort(
    floating* values,
    intp_t* indices,
    intp_t n,
    intp_t maxd,
) noexcept nogil:
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
                dual_swap(values, indices, i, l)
                i += 1
                l += 1
            elif values[i] > pivot:
                r -= 1
                dual_swap(values, indices, i, r)
            else:
                i += 1

        introsort(values, indices, l, maxd)
        values += r
        indices += r
        n -= r


cdef int simultaneous_sort(
    floating* values,
    intp_t* indices,
    intp_t size,
) noexcept nogil:
    """
    Perform an introsort on the values array as to sort them ascendingly.
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
    if size <= 1:
        return 0
    introsort(values, indices, size, 2 * <intp_t>log2(size))
    return 0
