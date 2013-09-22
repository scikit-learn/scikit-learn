"""
Uses C++ map containers for fast dict-like behavior with keys being
integers, and values float.
"""
# Author: Gael Varoquaux
# License: BSD

cimport cython

# C++
from cython.operator cimport dereference as deref, preincrement as inc, \
    predecrement as dec
from libcpp.utility cimport pair
from libcpp.map cimport map as cpp_map

import numpy as np

# Import the C-level symbols of numpy
cimport numpy as np

# Numpy must be initialized. When using numpy from C or Cython you must
# _always_ do that, or you will have segfaults
np.import_array()

#DTYPE = np.float64
#ctypedef np.float64_t DTYPE_t

#ITYPE = np.intp
#ctypedef np.intp_t ITYPE_t

###############################################################################
# An object to be used in Python

# Lookup is faster than dict (up to 10 times), and so is full traversal
# (up to 50 times), and assignment (up to 6 times), but creation is
# slower (up to 3 times). Also, a large benefit is that memory
# consumption is reduced a lot compared to a Python dict

cdef class IntFloatDict:
    #cdef cpp_map[ITYPE_t, DTYPE_t] my_map

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def __init__(self, np.ndarray[ITYPE_t, ndim=1] keys,
                       np.ndarray[DTYPE_t, ndim=1] values):
        cdef int i
        cdef int size = values.size
        # Should check that sizes for keys and values are equal, and
        # after should boundcheck(False)
        for i in range(size):
            self.my_map[keys[i]] = values[i]

    def __len__(self):
        return self.my_map.size()

    def __getitem__(self, int key):
        cdef cpp_map[ITYPE_t, DTYPE_t].iterator it = self.my_map.find(key)
        if it == self.my_map.end():
            # The key is not in the dict
            raise KeyError('%i' % key)
        return deref(it).second

    def __setitem__(self, int key, float value):
        self.my_map[key] = value

    def __iter__(self):
        cdef cpp_map[ITYPE_t, DTYPE_t].iterator it = self.my_map.begin()
        cdef cpp_map[ITYPE_t, DTYPE_t].iterator end = self.my_map.end()
        while it != end:
            yield deref(it).first, deref(it).second
            inc(it)

    def update(self, IntFloatDict other):
        cdef cpp_map[ITYPE_t, DTYPE_t].iterator it = other.my_map.begin()
        cdef cpp_map[ITYPE_t, DTYPE_t].iterator end = other.my_map.end()
        while it != end:
            self.my_map[deref(it).first] = deref(it).second
            inc(it)

    def copy(self):
        cdef IntFloatDict out_obj = IntFloatDict.__new__(IntFloatDict)
        # The '=' operator is a copy operator for C++ maps
        out_obj.my_map = self.my_map
        return out_obj

    def append(self, ITYPE_t key, DTYPE_t value):
        cdef cpp_map[ITYPE_t, DTYPE_t].iterator end = self.my_map.end()
        # Decrement the iterator
        dec(end)
        # Construct our arguments
        cdef pair[ITYPE_t, DTYPE_t] args
        args.first = key
        args.second = value
        self.my_map.insert(end, args)


###############################################################################
# operation on dict

def argmin(IntFloatDict d):
    cdef cpp_map[ITYPE_t, DTYPE_t].iterator it = d.my_map.begin()
    cdef cpp_map[ITYPE_t, DTYPE_t].iterator end = d.my_map.end()
    cdef ITYPE_t min_key
    cdef DTYPE_t min_value = np.inf
    while it != end:
        if deref(it).second < min_value:
            min_value = deref(it).second
            min_key = deref(it).first
        inc(it)
    return min_key, min_value

