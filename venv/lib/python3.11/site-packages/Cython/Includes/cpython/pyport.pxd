cdef extern from "Python.h":
    ctypedef int int32_t
    ctypedef int int64_t
    ctypedef unsigned int uint32_t
    ctypedef unsigned int uint64_t

    const Py_ssize_t PY_SSIZE_T_MIN
    const Py_ssize_t PY_SSIZE_T_MAX
