"""
Cython implementation of (parts of) the standard library time module.

Note: On implementations that lack a C-API for monotonic/perfcounter clocks
(like PyPy), the fallback code uses the system clock which may return absolute
time values from a different value range, differing from those provided by
Python's "time" module.
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
        static CYTHON_INLINE __Pyx_PyTime_t __Pyx_PyTime_TimeRaw(void) {
            __Pyx_PyTime_t tic;
            (void) PyTime_TimeRaw(&tic);
            return tic;
        }
    #else
        #define __Pyx_PyTime_TimeRaw()  _PyTime_GetSystemClock()
    #endif

    #if PY_VERSION_HEX >= 0x030d00b1 || defined(PyTime_MonotonicRaw)
        static CYTHON_INLINE __Pyx_PyTime_t __Pyx_PyTime_MonotonicRaw(void) {
            __Pyx_PyTime_t tic;
            (void) PyTime_MonotonicRaw(&tic);
            return tic;
        }
    #elif CYTHON_COMPILING_IN_PYPY && !defined(_PyTime_GetMonotonicClock)
        #define __Pyx_PyTime_MonotonicRaw()  _PyTime_GetSystemClock()
    #else
        #define __Pyx_PyTime_MonotonicRaw()  _PyTime_GetMonotonicClock()
    #endif

    #if PY_VERSION_HEX >= 0x030d00b1 || defined(PyTime_PerfCounterRaw)
        static CYTHON_INLINE __Pyx_PyTime_t __Pyx_PyTime_PerfCounterRaw(void) {
            __Pyx_PyTime_t tic;
            (void) PyTime_PerfCounterRaw(&tic);
            return tic;
        }
    #elif CYTHON_COMPILING_IN_PYPY && !defined(_PyTime_GetPerfCounter)
        #define __Pyx_PyTime_PerfCounterRaw()  __Pyx_PyTime_MonotonicRaw()
    #else
        #define __Pyx_PyTime_PerfCounterRaw()  _PyTime_GetPerfCounter()
    #endif

    #if PY_VERSION_HEX >= 0x030d00b1 || defined(PyTime_AsSecondsDouble)
        #define __Pyx_PyTime_AsSecondsDouble(t)  PyTime_AsSecondsDouble(t)
    #else
        #define __Pyx_PyTime_AsSecondsDouble(t)  _PyTime_AsSecondsDouble(t)
    #endif
    """
    ctypedef int64_t PyTime_t "__Pyx_PyTime_t"

    ctypedef PyTime_t _PyTime_t "__Pyx_PyTime_t"  # legacy, use "PyTime_t" instead
    PyTime_t PyTime_TimeUnchecked "__Pyx_PyTime_TimeRaw" () noexcept nogil  # legacy, use "PyTime_TimeRaw" instead

    PyTime_t PyTime_TimeRaw "__Pyx_PyTime_TimeRaw" () noexcept nogil
    PyTime_t PyTime_MonotonicRaw "__Pyx_PyTime_MonotonicRaw" () noexcept nogil
    PyTime_t PyTime_PerfCounterRaw "__Pyx_PyTime_PerfCounterRaw" () noexcept nogil
    double PyTime_AsSecondsDouble "__Pyx_PyTime_AsSecondsDouble" (PyTime_t t) noexcept nogil


from libc.time cimport (
    tm,
    time_t,
    localtime as libc_localtime,
)


cdef inline double time() noexcept nogil:
    cdef PyTime_t tic = PyTime_TimeRaw()
    return PyTime_AsSecondsDouble(tic)


cdef inline int64_t time_ns() noexcept nogil:
    return <int64_t> PyTime_TimeRaw()


cdef inline double perf_counter() noexcept nogil:
    cdef PyTime_t tic = PyTime_PerfCounterRaw()
    return PyTime_AsSecondsDouble(tic)


cdef inline int64_t perf_counter_ns() noexcept nogil:
    return <int64_t> PyTime_PerfCounterRaw()


cdef inline double monotonic() noexcept nogil:
    cdef PyTime_t tic = PyTime_MonotonicRaw()
    return PyTime_AsSecondsDouble(tic)


cdef inline int64_t monotonic_ns() noexcept nogil:
    return <int64_t> PyTime_MonotonicRaw()


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
