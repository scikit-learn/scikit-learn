"""
Cython implementation of (parts of) the standard library time module.
"""

from libc.stdint cimport int64_t
from cpython.exc cimport PyErr_SetFromErrno

cdef extern from *:
    """
    #if PY_VERSION_HEX >= 0x030d00b1 || defined(PyTime_t)
        #define __Pyx_PyTime_t PyTime_t
    #else
        #define __Pyx_PyTime_t _PyTime_t
    #endif

    #if PY_VERSION_HEX >= 0x030d00b1 || defined(PyTime_TimeRaw)
        static CYTHON_INLINE __Pyx_PyTime_t __Pyx_PyTime_TimeUnchecked(void) {
            __Pyx_PyTime_t tic;
            (void) PyTime_TimeRaw(&tic);
            return tic;
        }
    #else
        #define __Pyx_PyTime_TimeUnchecked()  _PyTime_GetSystemClock()
    #endif

    #if PY_VERSION_HEX >= 0x030d00b1 || defined(PyTime_AsSecondsDouble)
        #define __Pyx_PyTime_AsSecondsDouble(t)  PyTime_AsSecondsDouble(t)
    #else
        #define __Pyx_PyTime_AsSecondsDouble(t)  _PyTime_AsSecondsDouble(t)
    #endif
    """
    ctypedef int64_t PyTime_t "__Pyx_PyTime_t"
    ctypedef int64_t _PyTime_t "__Pyx_PyTime_t"

    _PyTime_t _PyTime_GetSystemClock "__Pyx_PyTime_TimeUnchecked" () nogil
    _PyTime_t PyTime_TimeUnchecked "__Pyx_PyTime_TimeUnchecked" () nogil

    double _PyTime_AsSecondsDouble "__Pyx_PyTime_AsSecondsDouble" (_PyTime_t t) nogil
    double PyTime_AsSecondsDouble "__Pyx_PyTime_AsSecondsDouble" (_PyTime_t t) nogil

from libc.time cimport (
    tm,
    time_t,
    localtime as libc_localtime,
)


cdef inline double time() noexcept nogil:
    cdef:
        _PyTime_t tic

    tic = PyTime_TimeUnchecked()
    return PyTime_AsSecondsDouble(tic)


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
    result.tm_wday = <int> ((result.tm_wday + 6) % 7)
    result.tm_yday += 1
    return result[0]
