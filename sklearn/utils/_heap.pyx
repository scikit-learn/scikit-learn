from cython cimport floating, integral, numeric

from ._typedefs cimport ITYPE_t

cdef inline void dual_swap(floating* darr, ITYPE_t* iarr,
                           ITYPE_t a, ITYPE_t b) nogil:
    """Swap the values at index a and b of both darr and iarr"""
    cdef floating dtmp = darr[a]
    darr[a] = darr[b]
    darr[b] = dtmp

    cdef ITYPE_t itmp = iarr[a]
    iarr[a] = iarr[b]
    iarr[b] = itmp

cdef int simultaneous_sort(
    floating* values,
    ITYPE_t* indices,
    ITYPE_t size
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
    # TODO: In order to support discrete distance metrics, we need to have a
    # simultaneous sort which breaks ties on indices when distances are identical.
    # The best might be using a std::stable_sort and a Comparator which might need
    # an Array of Structures (AoS) instead of the Structure of Arrays (SoA)
    # currently used.
    cdef:
        ITYPE_t pivot_idx, i, store_idx
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


cdef inline int heap_push(
    floating* values,
    ITYPE_t* indices,
    ITYPE_t size,
    floating val,
    ITYPE_t val_idx,
) nogil:
    """Push a tuple (val, val_idx) onto a fixed-size max-heap.

    The max-heap is represented as a Structure of Arrays where:
     - values is the array containing the data to construct the heap with
     - indices is the array containing the indices (meta-data) of each value

    Notes
    -----
    Arrays are manipulated via a pointer to there first element and their size
    as to ease the processing of dynamically allocated buffers.

    For instance, in pseudo-code:

        values = [1.2, 0.4, 0.1],
        indices = [42, 1, 5],
        heap_push(
            values=values,
            indices=indices,
            size=3,
            val=0.2,
            val_idx=4,
        )

    will modify values and indices inplace, giving at the end of the call:

        values  == [0.4, 0.2, 0.1]
        indices == [1, 4, 5]

    """
    cdef:
        ITYPE_t current_idx, left_child_idx, right_child_idx, swap_idx

    # Check if val should be in heap
    if val >= values[0]:
        return 0

    # Insert val at position zero
    values[0] = val
    indices[0] = val_idx

    # Descend the heap, swapping values until the max heap criterion is met
    current_idx = 0
    while True:
        left_child_idx = 2 * current_idx + 1
        right_child_idx = left_child_idx + 1

        if left_child_idx >= size:
            break
        elif right_child_idx >= size:
            if values[left_child_idx] > val:
                swap_idx = left_child_idx
            else:
                break
        elif values[left_child_idx] >= values[right_child_idx]:
            if val < values[left_child_idx]:
                swap_idx = left_child_idx
            else:
                break
        else:
            if val < values[right_child_idx]:
                swap_idx = right_child_idx
            else:
                break

        values[current_idx] = values[swap_idx]
        indices[current_idx] = indices[swap_idx]

        current_idx = swap_idx

    values[current_idx] = val
    indices[current_idx] = val_idx

    return 0
