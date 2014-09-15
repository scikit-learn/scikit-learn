""" Implementation of the Fast Hadamard Transform.

https://en.wikipedia.org/wiki/Hadamard_transform

This module supplies a single dimensional and two-dimensional row-wise
implementation. Both are non-normalized, operate in-place and can only handle
the double/float64 type.

Inspired by a Python-C-API implementation at:

https://github.com/nbarbey/fht

"""

import numpy as np
cimport numpy as np
cimport cython
from libc.math cimport log2


def is_power_of_two(input_integer):
    """ Test if an integer is a power of two. """
    if input_integer == 1:
        return False
    return input_integer != 0 and ((input_integer & (input_integer - 1)) == 0)

#DTYPE = np.double
ctypedef np.double_t DTYPE_t

def pure_python_fht(array_):
    """ Pure Python implementation for educational purposes. """
    bit = length = len(array_)
    for _ in xrange(int(np.log2(length))):
        bit >>= 1
        for i in xrange(length):
            if i & bit == 0:
                j = i | bit
                temp = array_[i]
                array_[i] += array_[j]
                array_[j] = temp - array_[j]


def fht(np.ndarray[DTYPE_t] array_):
    """ Single dimensional FHT. """
    if not is_power_of_two(array_.shape[0]):
        raise ValueError('Length of input for fht must be a power of two')
    else:
        _fht(array_)

@cython.boundscheck(False)
cdef _fht(np.ndarray[DTYPE_t, ndim=1] array_):
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



def fht2(np.ndarray[DTYPE_t, ndim=2] array_):
    """ Two dimensional row-wise FHT. """
    if not is_power_of_two(array_.shape[1]):
        raise ValueError('Length of rows for fht2 must be a power of two')
    else:
        _fht2(array_)

@cython.boundscheck(False)
cdef _fht2(np.ndarray[DTYPE_t, ndim=2] array_):
    cdef unsigned int bit, length, _, i, j, n
    cdef double temp
    n = array_.shape[0]
    for x in xrange(n):
        # TODO: This call still shows up as yellow in cython -a presumably due
        # to the [] access, but the array_ is already typed...
        _fht(array_[x])
