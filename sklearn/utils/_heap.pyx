#!python
#cython: boundscheck=False
#cython: wraparound=False
#cython: cdivision=True


from cython cimport floating, integral, numeric

from ._typedefs cimport ITYPE_t

cdef inline void dual_swap(floating* darr, ITYPE_t* iarr,
                           ITYPE_t i1, ITYPE_t i2) nogil:
    """swap the values at inex i1 and i2 of both darr and iarr"""
    cdef floating dtmp = darr[i1]
    darr[i1] = darr[i2]
    darr[i2] = dtmp

    cdef ITYPE_t itmp = iarr[i1]
    iarr[i1] = iarr[i2]
    iarr[i2] = itmp

cdef int _simultaneous_sort(
    floating* dist,
    ITYPE_t* idx,
    ITYPE_t size
) nogil except -1:
    """
    Perform a recursive quicksort on the dist array, simultaneously
    performing the same swaps on the idx array.

    TODO: test if the following algorithms are better:
      - introselect via std::nth_element
      - heap-sort-like
    """
    cdef:
        ITYPE_t pivot_idx, i, store_idx
        floating pivot_val

    # in the small-array case, do things efficiently
    if size <= 1:
        pass
    elif size == 2:
        if dist[0] > dist[1]:
            dual_swap(dist, idx, 0, 1)
    elif size == 3:
        if dist[0] > dist[1]:
            dual_swap(dist, idx, 0, 1)
        if dist[1] > dist[2]:
            dual_swap(dist, idx, 1, 2)
            if dist[0] > dist[1]:
                dual_swap(dist, idx, 0, 1)
    else:
        # Determine the pivot using the median-of-three rule.
        # The smallest of the three is moved to the beginning of the array,
        # the middle (the pivot value) is moved to the end, and the largest
        # is moved to the pivot index.
        pivot_idx = size // 2
        if dist[0] > dist[size - 1]:
            dual_swap(dist, idx, 0, size - 1)
        if dist[size - 1] > dist[pivot_idx]:
            dual_swap(dist, idx, size - 1, pivot_idx)
            if dist[0] > dist[size - 1]:
                dual_swap(dist, idx, 0, size - 1)
        pivot_val = dist[size - 1]

        # partition indices about pivot.  At the end of this operation,
        # pivot_idx will contain the pivot value, everything to the left
        # will be smaller, and everything to the right will be larger.
        store_idx = 0
        for i in range(size - 1):
            if dist[i] < pivot_val:
                dual_swap(dist, idx, i, store_idx)
                store_idx += 1
        dual_swap(dist, idx, store_idx, size - 1)
        pivot_idx = store_idx

        # recursively sort each side of the pivot
        if pivot_idx > 1:
            _simultaneous_sort(dist, idx, pivot_idx)
        if pivot_idx + 2 < size:
            _simultaneous_sort(dist + pivot_idx + 1,
                               idx + pivot_idx + 1,
                               size - pivot_idx - 1)
    return 0


cdef inline int _push(
    floating* dist,
    ITYPE_t* idx,
    ITYPE_t size,
    floating val,
    ITYPE_t i_val,
) nogil except -1:
    """push (val, i_val) into the heap (dist, idx) of the given size"""
    cdef:
        ITYPE_t current_idx, left_child_idx, right_child_idx, swap_idx

    # check if val should be in heap
    if val > dist[0]:
        return 0

    # insert val at position zero
    dist[0] = val
    idx[0] = i_val

    # descend the heap, swapping values until the max heap criterion is met
    current_idx = 0
    while True:
        left_child_idx = 2 * current_idx + 1
        right_child_idx = left_child_idx + 1

        if left_child_idx >= size:
            break
        elif right_child_idx >= size:
            if dist[left_child_idx] > val:
                swap_idx = left_child_idx
            else:
                break
        elif dist[left_child_idx] >= dist[right_child_idx]:
            if val < dist[left_child_idx]:
                swap_idx = left_child_idx
            else:
                break
        else:
            if val < dist[right_child_idx]:
                swap_idx = right_child_idx
            else:
                break

        dist[current_idx] = dist[swap_idx]
        idx[current_idx] = idx[swap_idx]

        current_idx = swap_idx

    dist[current_idx] = val
    idx[current_idx] = i_val

    return 0
