from cpython.object cimport PyObject
from cpython.version cimport PY_VERSION_HEX

cdef extern from *:
    """
    #if CYTHON_COMPILING_IN_LIMITED_API
    #ifdef _MSC_VER
    #pragma message ("This module uses CPython specific internals of 'datetime.datetime', which are not available in the limited API.")
    #else
    #warning This module uses CPython specific internals of 'datetime.datetime', which are not available in the limited API.
    #endif
    #endif
    """

cdef extern from "Python.h":
    ctypedef struct PyTypeObject:
        pass

cdef extern from "datetime.h":
    """
    #define __Pyx_DateTime_DateTimeWithFold(year, month, day, hour, minute, second, microsecond, tz, fold) \
        PyDateTimeAPI->DateTime_FromDateAndTimeAndFold(year, month, day, hour, minute, second, \
            microsecond, tz, fold, PyDateTimeAPI->DateTimeType)
    #define __Pyx_DateTime_TimeWithFold(hour, minute, second, microsecond, tz, fold) \
        PyDateTimeAPI->Time_FromTimeAndFold(hour, minute, second, microsecond, tz, fold, PyDateTimeAPI->TimeType)

    #define __Pyx_TimeZone_UTC PyDateTime_TimeZone_UTC
    #define __Pyx_TimeZone_FromOffsetAndName(offset, name) PyTimeZone_FromOffsetAndName(offset, name)

    /* Backport for Python < 3.10 */
    #if PY_VERSION_HEX < 0x030a00a1
        #ifndef PyDateTime_TIME_GET_TZINFO
            #define PyDateTime_TIME_GET_TZINFO(o) \
                ((((PyDateTime_Time*)o)->hastzinfo) ? ((PyDateTime_Time*)o)->tzinfo : Py_None)
        #endif
        #ifndef PyDateTime_DATE_GET_TZINFO
            #define PyDateTime_DATE_GET_TZINFO(o) \
                ((((PyDateTime_DateTime*)o)->hastzinfo) ? ((PyDateTime_DateTime*)o)->tzinfo : Py_None)
        #endif
    #endif
    """

    ctypedef extern class datetime.date[object PyDateTime_Date]:
        @property
        cdef inline int year(self) noexcept:
            return PyDateTime_GET_YEAR(self)

        @property
        cdef inline int month(self) noexcept:
            return PyDateTime_GET_MONTH(self)

        @property
        cdef inline int day(self) noexcept:
            return PyDateTime_GET_DAY(self)

    ctypedef extern class datetime.time[object PyDateTime_Time]:
        @property
        cdef inline int hour(self) noexcept:
            return PyDateTime_TIME_GET_HOUR(self)

        @property
        cdef inline int minute(self) noexcept:
            return PyDateTime_TIME_GET_MINUTE(self)

        @property
        cdef inline int second(self) noexcept:
            return PyDateTime_TIME_GET_SECOND(self)

        @property
        cdef inline int microsecond(self) noexcept:
            return PyDateTime_TIME_GET_MICROSECOND(self)

        @property
        cdef inline object tzinfo(self):
            return <object>PyDateTime_TIME_GET_TZINFO(self)

        @property
        cdef inline int fold(self) noexcept:
            # For Python < 3.6 this returns 0 no matter what
            return PyDateTime_TIME_GET_FOLD(self)

    ctypedef extern class datetime.datetime[object PyDateTime_DateTime]:
        @property
        cdef inline int year(self) noexcept:
            return PyDateTime_GET_YEAR(self)

        @property
        cdef inline int month(self) noexcept:
            return PyDateTime_GET_MONTH(self)

        @property
        cdef inline int day(self) noexcept:
            return PyDateTime_GET_DAY(self)

        @property
        cdef inline int hour(self) noexcept:
            return PyDateTime_DATE_GET_HOUR(self)

        @property
        cdef inline int minute(self) noexcept:
            return PyDateTime_DATE_GET_MINUTE(self)

        @property
        cdef inline int second(self) noexcept:
            return PyDateTime_DATE_GET_SECOND(self)

        @property
        cdef inline int microsecond(self) noexcept:
            return PyDateTime_DATE_GET_MICROSECOND(self)

        @property
        cdef inline object tzinfo(self):
            return <object>PyDateTime_DATE_GET_TZINFO(self)

        @property
        cdef inline int fold(self) noexcept:
            # For Python < 3.6 this returns 0 no matter what
            return PyDateTime_DATE_GET_FOLD(self)

    ctypedef extern class datetime.timedelta[object PyDateTime_Delta]:
        @property
        cdef inline int day(self) noexcept:
            return PyDateTime_DELTA_GET_DAYS(self)

        @property
        cdef inline int second(self) noexcept:
            return PyDateTime_DELTA_GET_SECONDS(self)

        @property
        cdef inline int microsecond(self) noexcept:
            return PyDateTime_DELTA_GET_MICROSECONDS(self)

    ctypedef extern class datetime.tzinfo[object PyDateTime_TZInfo]:
        pass

    ctypedef struct PyDateTime_Date:
        pass

    ctypedef struct PyDateTime_Time:
        unsigned char fold
        char hastzinfo
        PyObject *tzinfo

    ctypedef struct PyDateTime_DateTime:
        unsigned char fold
        char hastzinfo
        PyObject *tzinfo

    ctypedef struct PyDateTime_Delta:
        int days
        int seconds
        int microseconds

    # Define structure for C API.
    ctypedef struct PyDateTime_CAPI:
        # type objects
        PyTypeObject *DateType
        PyTypeObject *DateTimeType
        PyTypeObject *TimeType
        PyTypeObject *DeltaType
        PyTypeObject *TZInfoType

        # constructors
        date (*Date_FromDate)(int, int, int, PyTypeObject*)
        datetime (*DateTime_FromDateAndTime)(int, int, int, int, int, int, int, object, PyTypeObject*)
        time (*Time_FromTime)(int, int, int, int, object, PyTypeObject*)
        timedelta (*Delta_FromDelta)(int, int, int, int, PyTypeObject*)

        # constructors for the DB API
        datetime (*DateTime_FromTimestamp)(PyObject*, object, PyObject*)
        date (*Date_FromTimestamp)(PyObject*, object)

        # We cannot use the following because they do not compile in older Python versions.
        # Instead, we use datetime.h's macros here that we can backport in C.

        # Python 3.7+ constructors
        object (*TimeZone_FromTimeZone)(object offset, PyObject *name)

        # Python 3.7+ singletons
        PyObject *TimeZone_UTC

        # Python 3.6+ PEP 495 constructors
        datetime (*DateTime_FromDateAndTimeAndFold)(int, int, int, int, int, int, int, object, int, PyTypeObject*)
        time (*Time_FromTimeAndFold)(int, int, int ,int, object, int, PyTypeObject*)

    # Check type of the object.
    bint PyDate_Check(object op)
    bint PyDate_CheckExact(object op)

    bint PyDateTime_Check(object op)
    bint PyDateTime_CheckExact(object op)

    bint PyTime_Check(object op)
    bint PyTime_CheckExact(object op)

    bint PyDelta_Check(object op)
    bint PyDelta_CheckExact(object op)

    bint PyTZInfo_Check(object op)
    bint PyTZInfo_CheckExact(object op)

    # Getters for date and datetime (C macros).
    int PyDateTime_GET_YEAR(object o)
    int PyDateTime_GET_MONTH(object o)
    int PyDateTime_GET_DAY(object o)

    # Getters for datetime (C macros).
    int PyDateTime_DATE_GET_HOUR(object o)
    int PyDateTime_DATE_GET_MINUTE(object o)
    int PyDateTime_DATE_GET_SECOND(object o)
    int PyDateTime_DATE_GET_MICROSECOND(object o)
    int PyDateTime_DATE_GET_FOLD(object o)
    PyObject* PyDateTime_DATE_GET_TZINFO(object o)  # returns a borrowed reference

    # Getters for time (C macros).
    int PyDateTime_TIME_GET_HOUR(object o)
    int PyDateTime_TIME_GET_MINUTE(object o)
    int PyDateTime_TIME_GET_SECOND(object o)
    int PyDateTime_TIME_GET_MICROSECOND(object o)
    int PyDateTime_TIME_GET_FOLD(object o)
    PyObject* PyDateTime_TIME_GET_TZINFO(object o)  # returns a borrowed reference

    # Getters for timedelta (C macros).
    int PyDateTime_DELTA_GET_DAYS(object o)
    int PyDateTime_DELTA_GET_SECONDS(object o)
    int PyDateTime_DELTA_GET_MICROSECONDS(object o)

    # Constructors
    object PyTimeZone_FromOffset(object offset)
    object PyTimeZone_FromOffsetAndName(object offset, object name)

    # The above macros is Python 3.7+ so we use these instead
    object __Pyx_TimeZone_FromOffsetAndName(object offset, PyObject* name)

    # Constructors for the DB API
    datetime PyDateTime_FromTimeStamp(object args)
    date PyDate_FromTimeStamp(object args)

    # PEP 495 constructors but patched above to allow passing tz
    datetime __Pyx_DateTime_DateTimeWithFold(int, int, int, int, int, int, int, object, int)
    datetime __Pyx_DateTime_TimeWithFold(int, int, int ,int, object, int)

    # PyDateTime CAPI object.
    PyDateTime_CAPI *PyDateTimeAPI

    PyObject* PyDateTime_TimeZone_UTC

    # PyDateTime_TimeZone_UTC is Python 3.7+ so instead we use the following macro
    PyObject* __Pyx_TimeZone_UTC

    void PyDateTime_IMPORT()

# Datetime C API initialization function.
# You have to call it before any usage of DateTime CAPI functions.
cdef inline void import_datetime() noexcept:
    PyDateTime_IMPORT

# Create date object using DateTime CAPI factory function.
# Note, there are no range checks for any of the arguments.
cdef inline date date_new(int year, int month, int day):
    return PyDateTimeAPI.Date_FromDate(year, month, day, PyDateTimeAPI.DateType)

# Create time object using DateTime CAPI factory function
# Note, there are no range checks for any of the arguments.
cdef inline time time_new(int hour, int minute, int second, int microsecond, object tz, int fold=0):
    return __Pyx_DateTime_TimeWithFold(hour, minute, second, microsecond, tz, fold)

# Create datetime object using DateTime CAPI factory function.
# Note, there are no range checks for any of the arguments.
cdef inline datetime datetime_new(int year, int month, int day, int hour, int minute, int second, int microsecond, object tz, int fold=0):
    return __Pyx_DateTime_DateTimeWithFold(year, month, day, hour, minute, second, microsecond, tz, fold)

# Create timedelta object using DateTime CAPI factory function.
# Note, there are no range checks for any of the arguments.
cdef inline timedelta timedelta_new(int days, int seconds, int useconds):
    return PyDateTimeAPI.Delta_FromDelta(days, seconds, useconds, 1, PyDateTimeAPI.DeltaType)

# Create timedelta object using DateTime CAPI factory function.
cdef inline object timezone_new(object offset, object name=None):
    return __Pyx_TimeZone_FromOffsetAndName(offset, <PyObject*>name if name is not None else NULL)

# Create datetime object using DB API constructor.
cdef inline datetime datetime_from_timestamp(timestamp, tz=None):
    return PyDateTimeAPI.DateTime_FromTimestamp(
        <PyObject*>PyDateTimeAPI.DateTimeType, (timestamp, tz) if tz is not None else (timestamp,), NULL)

# Create date object using DB API constructor.
cdef inline date date_from_timestamp(timestamp):
    return PyDateTimeAPI.Date_FromTimestamp(<PyObject*>PyDateTimeAPI.DateType, (timestamp,))

# More recognizable getters for date/time/datetime/timedelta.
# There are no setters because datetime.h hasn't them.
# This is because of immutable nature of these objects by design.
# If you would change time/date/datetime/timedelta object you need to recreate.

# Get UTC singleton
cdef inline object get_utc():
    return <object>__Pyx_TimeZone_UTC

# Get tzinfo of time
cdef inline object time_tzinfo(object o):
    return <object>PyDateTime_TIME_GET_TZINFO(o)

# Get tzinfo of datetime
cdef inline object datetime_tzinfo(object o):
    return <object>PyDateTime_DATE_GET_TZINFO(o)

# Get year of date
cdef inline int date_year(object o) noexcept:
    return PyDateTime_GET_YEAR(o)

# Get month of date
cdef inline int date_month(object o) noexcept:
    return PyDateTime_GET_MONTH(o)

# Get day of date
cdef inline int date_day(object o) noexcept:
    return PyDateTime_GET_DAY(o)

# Get year of datetime
cdef inline int datetime_year(object o) noexcept:
    return PyDateTime_GET_YEAR(o)

# Get month of datetime
cdef inline int datetime_month(object o) noexcept:
    return PyDateTime_GET_MONTH(o)

# Get day of datetime
cdef inline int datetime_day(object o) noexcept:
    return PyDateTime_GET_DAY(o)

# Get hour of time
cdef inline int time_hour(object o) noexcept:
    return PyDateTime_TIME_GET_HOUR(o)

# Get minute of time
cdef inline int time_minute(object o) noexcept:
    return PyDateTime_TIME_GET_MINUTE(o)

# Get second of time
cdef inline int time_second(object o) noexcept:
    return PyDateTime_TIME_GET_SECOND(o)

# Get microsecond of time
cdef inline int time_microsecond(object o) noexcept:
    return PyDateTime_TIME_GET_MICROSECOND(o)

# Get fold of time
cdef inline int time_fold(object o) noexcept:
    # For Python < 3.6 this returns 0 no matter what
    return PyDateTime_TIME_GET_FOLD(o)

# Get hour of datetime
cdef inline int datetime_hour(object o) noexcept:
    return PyDateTime_DATE_GET_HOUR(o)

# Get minute of datetime
cdef inline int datetime_minute(object o) noexcept:
    return PyDateTime_DATE_GET_MINUTE(o)

# Get second of datetime
cdef inline int datetime_second(object o) noexcept:
    return PyDateTime_DATE_GET_SECOND(o)

# Get microsecond of datetime
cdef inline int datetime_microsecond(object o) noexcept:
    return PyDateTime_DATE_GET_MICROSECOND(o)

# Get fold of datetime
cdef inline int datetime_fold(object o) noexcept:
    # For Python < 3.6 this returns 0 no matter what
    return PyDateTime_DATE_GET_FOLD(o)

# Get days of timedelta
cdef inline int timedelta_days(object o) noexcept:
    return PyDateTime_DELTA_GET_DAYS(o)

# Get seconds of timedelta
cdef inline int timedelta_seconds(object o) noexcept:
    return PyDateTime_DELTA_GET_SECONDS(o)

# Get microseconds of timedelta
cdef inline int timedelta_microseconds(object o) noexcept:
    return PyDateTime_DELTA_GET_MICROSECONDS(o)

cdef inline double total_seconds(timedelta obj) noexcept:
    # Mirrors the "timedelta.total_seconds()" method.
    # Note that this implementation is not guaranteed to give *exactly* the same
    # result as the original method, due to potential differences in floating point rounding.
    cdef:
        double days, seconds, micros
    days = <double>PyDateTime_DELTA_GET_DAYS(obj)
    seconds = <double>PyDateTime_DELTA_GET_SECONDS(obj)
    micros = <double>PyDateTime_DELTA_GET_MICROSECONDS(obj)
    return days * 24 * 3600 + seconds + micros / 1_000_000
