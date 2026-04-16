/*
Copyright Rene Rivera 2015
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE_1_0.txt or copy at
http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_PREDEF_OS_ANDROID_H
#define BOOST_PREDEF_OS_ANDROID_H

#include <boost/predef/version_number.h>
#include <boost/predef/make.h>

/*`
[heading `BOOST_OS_ANDROID`]

NOTE: `BOOST_OS_ANDROID` is deprecated, and will be removed in a following release.
Please use `BOOST_PLAT_ANDROID` instead.

[@http://en.wikipedia.org/wiki/Android_%28operating_system%29 Android] operating system.

[table
    [[__predef_symbol__] [__predef_version__]]

    [[`__ANDROID__`] [__predef_detection__]]
    ]
 */

#define BOOST_OS_ANDROID BOOST_VERSION_NUMBER_NOT_AVAILABLE

#if !defined(BOOST_PREDEF_DETAIL_OS_DETECTED) && ( \
    defined(__ANDROID__) \
    )
#   undef BOOST_OS_ANDROID
#   define BOOST_OS_ANDROID BOOST_VERSION_NUMBER_AVAILABLE
#endif

#if BOOST_OS_ANDROID
#   define BOOST_OS_ANDROID_AVAILABLE
#   include <boost/predef/detail/os_detected.h>
#endif

#define BOOST_OS_ANDROID_NAME "Android"

#endif

#include <boost/predef/detail/test.h>
BOOST_PREDEF_DECLARE_TEST(BOOST_OS_ANDROID,BOOST_OS_ANDROID_NAME)
