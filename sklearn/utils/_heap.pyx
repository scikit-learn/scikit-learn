from cython cimport floating

from sklearn.utils._typedefs cimport intp_t


cdef inline int heap_push(
    floating* values,
    intp_t* indices,
    intp_t size,
    floating val,
    intp_t val_idx,
) noexcept nogil:
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
        intp_t current_idx, left_child_idx, right_child_idx, swap_idx

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
