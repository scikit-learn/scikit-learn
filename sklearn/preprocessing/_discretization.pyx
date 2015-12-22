import cython

import numpy as np
cimport numpy as np
ctypedef np.float64_t DOUBLE  # This is simply `double`

@cython.boundscheck(False)
@cython.wraparound(False)
def binsearch(contcolumn, zero_interval, search_points):
    # TODO: Add sparse array support

    cdef DOUBLE[::1] _search_points = search_points.astype(float)

    _search_points = np.insert(_search_points, 0, -np.inf)
    _search_points = np.append(_search_points, np.inf)

    cdef:
        np.ndarray[DOUBLE, ndim=1, mode='c'] values = contcolumn.astype(float)
        np.ndarray[DOUBLE, ndim=1, mode="c"] output = np.zeros(contcolumn.shape[0])
        DOUBLE lower = <DOUBLE> zero_interval[0]
        DOUBLE upper = <DOUBLE> zero_interval[1]
        unsigned int i = 0, keyindex
        long[::1] keys = np.arange(1, search_points.shape[0])
        DOUBLE value
    print _search_points

    #with nogil:
    for value_index in range(values.shape[0]):
        print "hello world"
        value = values[value_index]
        if lower <= value < upper:
            output[i] = 0
        else:
            keyindex = binary_search(value, _search_points)
            output[i] = keys[keyindex]
        i += 1
    print output
    return output

@cython.boundscheck(False)
@cython.wraparound(False)
cdef int binary_search(DOUBLE value, DOUBLE[::1] search_points): # except -1: #nogil except -1:
    cdef:
        int lower = 0, upper = search_points.shape[0]
        int mid
    print search_points

    while lower < upper:
        mid = (lower + upper) / 2
        if search_points[mid] <= value < search_points[mid + 1]:
            return mid
        elif search_points[mid] < value:
            lower = mid + 1
        else:
            upper = mid
    return -1
