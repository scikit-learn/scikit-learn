import numpy as np
cimport numpy as np
cimport cython
from libc.math cimport log2


DTYPE = np.double
ctypedef np.double_t DTYPE_t

@cython.boundscheck(False)
def fht(np.ndarray[DTYPE_t] array_):
    cdef unsigned int bit, length, _, i, j
    cdef double temp
    bit = length = len(array_)
    for _ in xrange(int(log2(length))):
        bit >>= 1
        for i in xrange(length):
            if i & bit == 0:
                j = i | bit
                temp = array_[i]
                array_[i] += array_[j]
                array_[j] = temp - array_[j]

@cython.boundscheck(False)
def fht2(np.ndarray[DTYPE_t, ndim=2] array_):
    cdef unsigned int bit, length, _, i, j, n
    cdef double temp
    n = array_.shape[0]
    for x in xrange(n):
        bit = length = array_.shape[1]
        for _ in xrange(int(log2(length))):
            bit >>= 1
            for i in xrange(length):
                if i & bit == 0:
                    j = i | bit
                    temp = array_[x, i]
                    array_[x, i] += array_[x, j]
                    array_[x, j] = temp - array_[x, j]

