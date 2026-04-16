/*
    pybind11/chrono.h: Transparent conversion between std::chrono and python's datetime

    Copyright (c) 2016 Trent Houliston <trent@houliston.me> and
                       Wenzel Jakob <wenzel.jakob@epfl.ch>

    All rights reserved. Use of this source code is governed by a
    BSD-style license that can be found in the LICENSE file.
*/

#pragma once

#include "pybind11.h"

#include <chrono>
#include <cmath>
#include <ctime>
#include <datetime.h>
#include <mutex>

PYBIND11_NAMESPACE_BEGIN(PYBIND11_NAMESPACE)
PYBIND11_NAMESPACE_BEGIN(detail)

template <typename type>
class duration_caster {
public:
    using rep = typename type::rep;
    using period = typename type::period;

    // signed 25 bits required by the standard.
    using days = std::chrono::duration<int_least32_t, std::ratio<86400>>;

    bool load(handle src, bool) {
        using namespace std::chrono;

        // Lazy initialise the PyDateTime import
        if (!PyDateTimeAPI) {
            PyDateTime_IMPORT;
        }

        if (!src) {
            return false;
        }
        // If invoked with datetime.delta object
        if (PyDelta_Check(src.ptr())) {
            value = type(duration_cast<duration<rep, period>>(
                days(PyDateTime_DELTA_GET_DAYS(src.ptr()))
                + seconds(PyDateTime_DELTA_GET_SECONDS(src.ptr()))
                + microseconds(PyDateTime_DELTA_GET_MICROSECONDS(src.ptr()))));
            return true;
        }
        // If invoked with a float we assume it is seconds and convert
        if (PyFloat_Check(src.ptr())) {
            value = type(duration_cast<duration<rep, period>>(
                duration<double>(PyFloat_AsDouble(src.ptr()))));
            return true;
        }
        return false;
    }

    // If this is a duration just return it back
    static const std::chrono::duration<rep, period> &
    get_duration(const std::chrono::duration<rep, period> &src) {
        return src;
    }
    static const std::chrono::duration<rep, period> &
    get_duration(const std::chrono::duration<rep, period> &&) = delete;

    // If this is a time_point get the time_since_epoch
    template <typename Clock>
    static std::chrono::duration<rep, period>
    get_duration(const std::chrono::time_point<Clock, std::chrono::duration<rep, period>> &src) {
        return src.time_since_epoch();
    }

    static handle cast(const type &src, return_value_policy /* policy */, handle /* parent */) {
        using namespace std::chrono;

        // Use overloaded function to get our duration from our source
        // Works out if it is a duration or time_point and get the duration
        auto d = get_duration(src);

        // Lazy initialise the PyDateTime import
        if (!PyDateTimeAPI) {
            PyDateTime_IMPORT;
        }

        // Declare these special duration types so the conversions happen with the correct
        // primitive types (int)
        using dd_t = duration<int, std::ratio<86400>>;
        using ss_t = duration<int, std::ratio<1>>;
        using us_t = duration<int, std::micro>;

        auto dd = duration_cast<dd_t>(d);
        auto subd = d - dd;
        auto ss = duration_cast<ss_t>(subd);
        auto us = duration_cast<us_t>(subd - ss);
        return PyDelta_FromDSU(dd.count(), ss.count(), us.count());
    }

    PYBIND11_TYPE_CASTER(type, const_name("datetime.timedelta"));
};

inline std::tm *localtime_thread_safe(const std::time_t *time, std::tm *buf) {
#if (defined(__STDC_LIB_EXT1__) && defined(__STDC_WANT_LIB_EXT1__)) || defined(_MSC_VER)
    if (localtime_s(buf, time))
        return nullptr;
    return buf;
#else
    static std::mutex mtx;
    std::lock_guard<std::mutex> lock(mtx);
    std::tm *tm_ptr = std::localtime(time);
    if (tm_ptr != nullptr) {
        *buf = *tm_ptr;
    }
    return tm_ptr;
#endif
}

// This is for casting times on the system clock into datetime.datetime instances
template <typename Duration>
class type_caster<std::chrono::time_point<std::chrono::system_clock, Duration>> {
public:
    using type = std::chrono::time_point<std::chrono::system_clock, Duration>;
    bool load(handle src, bool) {
        using namespace std::chrono;

        // Lazy initialise the PyDateTime import
        if (!PyDateTimeAPI) {
            PyDateTime_IMPORT;
        }

        if (!src) {
            return false;
        }

        std::tm cal;
        microseconds msecs;

        if (PyDateTime_Check(src.ptr())) {
            cal.tm_sec = PyDateTime_DATE_GET_SECOND(src.ptr());
            cal.tm_min = PyDateTime_DATE_GET_MINUTE(src.ptr());
            cal.tm_hour = PyDateTime_DATE_GET_HOUR(src.ptr());
            cal.tm_mday = PyDateTime_GET_DAY(src.ptr());
            cal.tm_mon = PyDateTime_GET_MONTH(src.ptr()) - 1;
            cal.tm_year = PyDateTime_GET_YEAR(src.ptr()) - 1900;
            cal.tm_isdst = -1;
            msecs = microseconds(PyDateTime_DATE_GET_MICROSECOND(src.ptr()));
        } else if (PyDate_Check(src.ptr())) {
            cal.tm_sec = 0;
            cal.tm_min = 0;
            cal.tm_hour = 0;
            cal.tm_mday = PyDateTime_GET_DAY(src.ptr());
            cal.tm_mon = PyDateTime_GET_MONTH(src.ptr()) - 1;
            cal.tm_year = PyDateTime_GET_YEAR(src.ptr()) - 1900;
            cal.tm_isdst = -1;
            msecs = microseconds(0);
        } else if (PyTime_Check(src.ptr())) {
            cal.tm_sec = PyDateTime_TIME_GET_SECOND(src.ptr());
            cal.tm_min = PyDateTime_TIME_GET_MINUTE(src.ptr());
            cal.tm_hour = PyDateTime_TIME_GET_HOUR(src.ptr());
            cal.tm_mday = 1;  // This date (day, month, year) = (1, 0, 70)
            cal.tm_mon = 0;   // represents 1-Jan-1970, which is the first
            cal.tm_year = 70; // earliest available date for Python's datetime
            cal.tm_isdst = -1;
            msecs = microseconds(PyDateTime_TIME_GET_MICROSECOND(src.ptr()));
        } else {
            return false;
        }

        value = time_point_cast<Duration>(system_clock::from_time_t(std::mktime(&cal)) + msecs);
        return true;
    }

    static handle cast(const std::chrono::time_point<std::chrono::system_clock, Duration> &src,
                       return_value_policy /* policy */,
                       handle /* parent */) {
        using namespace std::chrono;

        // Lazy initialise the PyDateTime import
        if (!PyDateTimeAPI) {
            PyDateTime_IMPORT;
        }

        // Get out microseconds, and make sure they are positive, to avoid bug in eastern
        // hemisphere time zones (cfr. https://github.com/pybind/pybind11/issues/2417)
        using us_t = duration<int, std::micro>;
        auto us = duration_cast<us_t>(src.time_since_epoch() % seconds(1));
        if (us.count() < 0) {
            us += duration_cast<us_t>(seconds(1));
        }

        // Subtract microseconds BEFORE `system_clock::to_time_t`, because:
        // > If std::time_t has lower precision, it is implementation-defined whether the value is
        // rounded or truncated. (https://en.cppreference.com/w/cpp/chrono/system_clock/to_time_t)
        std::time_t tt
            = system_clock::to_time_t(time_point_cast<system_clock::duration>(src - us));

        std::tm localtime;
        std::tm *localtime_ptr = localtime_thread_safe(&tt, &localtime);
        if (!localtime_ptr) {
            throw cast_error("Unable to represent system_clock in local time");
        }
        return PyDateTime_FromDateAndTime(localtime.tm_year + 1900,
                                          localtime.tm_mon + 1,
                                          localtime.tm_mday,
                                          localtime.tm_hour,
                                          localtime.tm_min,
                                          localtime.tm_sec,
                                          us.count());
    }
    PYBIND11_TYPE_CASTER(type, const_name("datetime.datetime"));
};

// Other clocks that are not the system clock are not measured as datetime.datetime objects
// since they are not measured on calendar time. So instead we just make them timedeltas
// Or if they have passed us a time as a float we convert that
template <typename Clock, typename Duration>
class type_caster<std::chrono::time_point<Clock, Duration>>
    : public duration_caster<std::chrono::time_point<Clock, Duration>> {};

template <typename Rep, typename Period>
class type_caster<std::chrono::duration<Rep, Period>>
    : public duration_caster<std::chrono::duration<Rep, Period>> {};

PYBIND11_NAMESPACE_END(detail)
PYBIND11_NAMESPACE_END(PYBIND11_NAMESPACE)
