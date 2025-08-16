cimport cython

cdef long maxint

@cython.final
cdef class TransitionMap:
    cdef list map
    cdef dict special

    cpdef add(self, event, new_state)
    cpdef add_set(self, event, new_set)
    cpdef iteritems(self)
    cdef Py_ssize_t split(self, long code)
    cdef set get_special(self, event)
