import cython

import numpy as np
cimport numpy as np
ctypedef np.float64_t DOUBLE  # This is simply `double`

@cython.boundscheck(False)
@cython.wraparound(False)
def binsearch(contcolumn, zero_interval, DOUBLE[::1] search_points):
    # TODO: Add sparse array support

    if search_points.min() < 0:
        search_points = np.insert(search_points, -np.inf)

    cdef:
        np.ndarray[DOUBLE, ndim=1, mode='c'] values = contcolumn.astype(float)
        np.ndarray[DOUBLE, ndim=1, mode="c"] output = np.zeros(len(contcolumn))
        DOUBLE lower = <DOUBLE> zero_interval[0]
        DOUBLE upper = <DOUBLE> zero_interval[1]
        int i = 0, keyindex
        int[::1] keys = np.arange(1, len(search_points))
        DOUBLE value

    with nogil:
        for value_index in range(values.shape[0]):
            value = values[value_index]
            if lower <= value < upper:
                output[i] = 0
            else:
                keyindex = binary_search(value, search_points)
                output[i] = keys[keyindex]
               # if keyindex == -1:
               #     raise ValueError() # Will never happen
            i += 1

@cython.boundscheck(False)
@cython.wraparound(False)
cdef int binary_search(DOUBLE value, DOUBLE[::1] search_points) nogil:
    cdef:
        int lower = 0, upper = search_points.shape[0]
        int mid
    #with nogil:
    while lower < upper:
        mid = (lower + upper) / 2
        if search_points[mid] <= value < search_points[mid + 1]:
            return mid
        elif search_points[mid] < value:
            lower = mid
        else:
            upper = mid
    return -1
