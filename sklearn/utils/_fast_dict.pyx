"""
Uses C++ map containers for fast dict-like behavior with keys being
integers, and values float.
"""
# Author: Gael Varoquaux
# License: BSD

# C++
from cython.operator cimport dereference as deref, preincrement as inc
from libcpp.utility cimport pair
from libcpp.map cimport map as cpp_map

import numpy as np

# Import the C-level symbols of numpy
cimport numpy as cnp

# Numpy must be initialized. When using numpy from C or Cython you must
# _always_ do that, or you will have segfaults
cnp.import_array()

#DTYPE = np.float64
#ctypedef cnp.float64_t DTYPE_t

#ITYPE = np.intp
#ctypedef cnp.intp_t ITYPE_t

###############################################################################
# An object to be used in Python

# Lookup is faster than dict (up to 10 times), and so is full traversal
# (up to 50 times), and assignment (up to 6 times), but creation is
# slower (up to 3 times). Also, a large benefit is that memory
# consumption is reduced a lot compared to a Python dict

cdef class IntFloatDict:

    def __init__(self, cnp.ndarray[ITYPE_t, ndim=1] keys,
                       cnp.ndarray[DTYPE_t, ndim=1] values):
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

    # Cython 0.20 generates buggy code below. Commenting this out for now
    # and relying on the to_arrays method
    #def __iter__(self):
    #    cdef cpp_map[ITYPE_t, DTYPE_t].iterator it = self.my_map.begin()
    #    cdef cpp_map[ITYPE_t, DTYPE_t].iterator end = self.my_map.end()
    #    while it != end:
    #        yield deref(it).first, deref(it).second
    #        inc(it)

    def __iter__(self):
        cdef int size = self.my_map.size()
        cdef ITYPE_t [:] keys = np.empty(size, dtype=np.intp)
        cdef DTYPE_t [:] values = np.empty(size, dtype=np.float64)
        self._to_arrays(keys, values)
        cdef int idx
        cdef ITYPE_t key
        cdef DTYPE_t value
        for idx in range(size):
            key = keys[idx]
            value = values[idx]
            yield key, value

    def to_arrays(self):
        """Return the key, value representation of the IntFloatDict
           object.

           Returns
           =======
           keys : ndarray, shape (n_items, ), dtype=int
                The indices of the data points
           values : ndarray, shape (n_items, ), dtype=float
                The values of the data points
        """
        cdef int size = self.my_map.size()
        cdef cnp.ndarray[ITYPE_t, ndim=1] keys = np.empty(size,
                                                         dtype=np.intp)
        cdef cnp.ndarray[DTYPE_t, ndim=1] values = np.empty(size,
                                                           dtype=np.float64)
        self._to_arrays(keys, values)
        return keys, values

    cdef _to_arrays(self, ITYPE_t [:] keys, DTYPE_t [:] values):
        # Internal version of to_arrays that takes already-initialized arrays
        cdef cpp_map[ITYPE_t, DTYPE_t].iterator it = self.my_map.begin()
        cdef cpp_map[ITYPE_t, DTYPE_t].iterator end = self.my_map.end()
        cdef int index = 0
        while it != end:
            keys[index] = deref(it).first
            values[index] = deref(it).second
            inc(it)
            index += 1

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
        # Construct our arguments
        cdef pair[ITYPE_t, DTYPE_t] args
        args.first = key
        args.second = value
        self.my_map.insert(args)


###############################################################################
# operation on dict

def argmin(IntFloatDict d):
    cdef cpp_map[ITYPE_t, DTYPE_t].iterator it = d.my_map.begin()
    cdef cpp_map[ITYPE_t, DTYPE_t].iterator end = d.my_map.end()
    cdef ITYPE_t min_key = -1
    cdef DTYPE_t min_value = np.inf
    while it != end:
        if deref(it).second < min_value:
            min_value = deref(it).second
            min_key = deref(it).first
        inc(it)
    return min_key, min_value
