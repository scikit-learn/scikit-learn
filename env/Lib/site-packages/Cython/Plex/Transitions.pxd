cimport cython

cdef long maxint

@cython.final
cdef class TransitionMap:
    cdef list map
    cdef dict special

    @cython.locals(i=cython.Py_ssize_t, j=cython.Py_ssize_t)
    cpdef add(self, event, new_state)

    @cython.locals(i=cython.Py_ssize_t, j=cython.Py_ssize_t)
    cpdef add_set(self, event, new_set)

    @cython.locals(i=cython.Py_ssize_t, n=cython.Py_ssize_t, else_set=cython.bint)
    cpdef iteritems(self)

    @cython.locals(map=list, lo=cython.Py_ssize_t, mid=cython.Py_ssize_t, hi=cython.Py_ssize_t)
    cdef split(self, long code)

    cdef get_special(self, event)
