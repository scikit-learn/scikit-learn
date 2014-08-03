import numpy as np
cimport numpy as np
cimport cython

DTYPE = np.double
ctypedef np.double_t DTYPE_t

@cython.boundscheck(False)
def fht(np.ndarray[DTYPE_t] array_):
    cdef int bit, length, _, i, j
    cdef double temp
    bit = length = len(array_)
    for _ in xrange(int(np.log2(length))):
        bit >>= 1
        for i in xrange(length):
            if i & bit == 0:
                j = i | bit
                temp = array_[i]
                array_[i] += array_[j]
                array_[j] = temp - array_[j]

