"""
Cython implementation of (parts of) the standard library time module.
"""

from libc.stdint cimport int64_t
from cpython.exc cimport PyErr_SetFromErrno

cdef extern from "Python.h":
    ctypedef int64_t _PyTime_t
    _PyTime_t _PyTime_GetSystemClock() nogil
    double _PyTime_AsSecondsDouble(_PyTime_t t) nogil

from libc.time cimport (
    tm,
    time_t,
    localtime as libc_localtime,
)


cdef inline double time() noexcept nogil:
    cdef:
        _PyTime_t tic

    tic = _PyTime_GetSystemClock()
    return _PyTime_AsSecondsDouble(tic)


cdef inline int _raise_from_errno() except -1 with gil:
    PyErr_SetFromErrno(RuntimeError)
    return <int> -1  # Let the C compiler know that this function always raises.


cdef inline tm localtime() except * nogil:
    """
    Analogue to the stdlib time.localtime.  The returned struct
    has some entries that the stdlib version does not: tm_gmtoff, tm_zone
    """
    cdef:
        time_t tic = <time_t>time()
        tm* result

    result = libc_localtime(&tic)
    if result is NULL:
        _raise_from_errno()
    # Fix 0-based date values (and the 1900-based year).
    # See tmtotuple() in https://github.com/python/cpython/blob/master/Modules/timemodule.c
    result.tm_year += 1900
    result.tm_mon += 1
    result.tm_wday = (result.tm_wday + 6) % 7
    result.tm_yday += 1
    return result[0]
