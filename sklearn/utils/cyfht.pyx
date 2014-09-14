cimport numpy as np
cimport cython
from libc.math cimport log2

def is_power_of_two(input_integer):
    """Test if an integer is a power of two"""
    if input_integer == 1:
        return False
    return input_integer != 0 and ((input_integer & (input_integer - 1)) == 0)

#DTYPE = np.double
ctypedef np.double_t DTYPE_t



def fht(np.ndarray[DTYPE_t] array_):
    """ Single dimensional fht. """
    if not is_power_of_two(array_.shape[0]):
        raise ValueError('Length of input for fht must be a power of two')
    else:
        _fht(array_)

@cython.boundscheck(False)
cdef _fht(np.ndarray[DTYPE_t] array_):
    cdef unsigned int bit, length, _, i, j
    cdef double temp
    bit = length = array_.shape[0]
    for _ in xrange(<unsigned int>(log2(length))):
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
        for _ in xrange(<unsigned int>(log2(length))):
            bit >>= 1
            for i in xrange(length):
                if i & bit == 0:
                    j = i | bit
                    temp = array_[x, i]
                    array_[x, i] += array_[x, j]
                    array_[x, j] = temp - array_[x, j]

