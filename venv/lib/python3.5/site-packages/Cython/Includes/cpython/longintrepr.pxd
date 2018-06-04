# Internals of the "long" type (Python 2) or "int" type (Python 3).
# This is not part of Python's published API.

cdef extern from "longintrepr.h":
    # Add explicit cast to avoid compiler warnings
    cdef _PyLong_New "(PyObject*)_PyLong_New"(Py_ssize_t s)

    ctypedef unsigned int digit
    ctypedef int sdigit  # Python >= 2.7 only

    ctypedef struct PyLongObject:
        digit* ob_digit

    cdef long PyLong_SHIFT
    cdef digit PyLong_BASE
    cdef digit PyLong_MASK
