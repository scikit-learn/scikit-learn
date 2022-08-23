from cython cimport floating


cdef inline void dual_swap(
    floating* darr,
    ITYPE_t* iarr,
    ITYPE_t a,
    ITYPE_t b,
) nogil:
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
    ITYPE_t size,
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
    Arrays are manipulated via a pointer to their first element and their size
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


cdef int simultaneous_radix_sort(
    ITYPE_t* values,
    ITYPE_t* indices,
    ITYPE_t size,
    ITYPE_t* value_copies,
    ITYPE_t* index_copies,
) nogil:
    ''' 
    Perform a radix sort on the values array as to sort them ascendingly.
    This simultaneously performs the swaps on both the values and the indices
    arrays.

    The numpy equivalent is:

        def simultaneous_sort(values, indices):
             i = np.argsort(values)
             return values[i], indices[i]
    
    The algorithm consists of repeatedly applying Counting Sort to smaller sets
    of bits of the values. For example, when sorting 32-bit integers, all
    values are first sorted with respect to their 8 left-most bits, then are
    resorted with respect to their 9-16 bits, then bits 17-24 and finally 25-32.

    Auxiliary arrays value_copies and index_copies are passed as arguments to
    avoid repeatedly allocating memory in multiple subsequent sorts.

    Notes
    -----
    Arrays are manipulated via a pointer to their first element and their size
    as to ease the processing of dynamically allocated buffers.

    Modified from https://pt.wikipedia.org/wiki/Radix_sort#C%C3%B3digo_em_C
    (in portuguese).
    '''
    cdef ITYPE_t block_size = 8   # Number of bits to use at each iteration
    cdef ITYPE_t n_buckets = 256  # 1 << n_bits
    cdef ITYPE_t mask = 255       # n_buckets - 1
    cdef ITYPE_t i, j, key
    cdef ITYPE_t max_value
    cdef ITYPE_t[256] counter

    # One extra pass to get the maximum value and enable early stopping
    max_value = values[0]
    for j in range(1, size):
        if values[j] > max_value:
            max_value = values[j]

    # The following is used multiple times bellow to get the ith set of bits
    # from value:
    #
    #     (value >> (i*block_size)) & mask
    #
    # i thus will indicate the bit partition we currently using to sort.
    # For instance, the number 5834326 will be partitioned as:
    #
    #    (0000000)(00101100)(10000011)(001010110)
    #       i=3      i=2        i=1      i=0

    i = 0
    while (max_value >> (i*block_size)) & mask:
        # Initialize buckets for counting digits (values in sets of bits).
        for j in range(n_buckets):
            counter[j] = 0
    
        # Count ocurrences of each bit combination in the current (ith) bit
        # ppartition.
        for j in range(size):
            counter[(values[j] >> (i*block_size)) & mask] += 1

        # Compute cumulative sum of counts, since they will indicate the
        # number's position in the sorted array.
        for j in range(1, n_buckets):
            counter[j] += counter[j-1]

        # Build sorted arrays with respect to the ith bit partition, storing
        # them in value_copies and index_copies
        for j in range(size-1, -1, -1):
            key = (values[j] >> (i*block_size)) & mask
            counter[key] -= 1
            value_copies[counter[key]] = values[j]
            index_copies[counter[key]] = indices[j]

        # Copy sorted arrays back to values and indices.
        for j in range(size):
            values[j] = value_copies[j]
            indices[j] = index_copies[j]

        i += 1

    return 0
