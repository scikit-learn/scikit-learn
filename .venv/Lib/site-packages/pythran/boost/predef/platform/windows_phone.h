/*
Copyright (c) Microsoft Corporation 2014
Copyright Rene Rivera 2015
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE_1_0.txt or copy at
http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_PREDEF_PLAT_WINDOWS_PHONE_H
#define BOOST_PREDEF_PLAT_WINDOWS_PHONE_H

#include <boost/predef/make.h>
#include <boost/predef/os/windows.h>
#include <boost/predef/platform/windows_uwp.h>
#include <boost/predef/version_number.h>

/*`
[heading `BOOST_PLAT_WINDOWS_PHONE`]

[@https://docs.microsoft.com/en-us/windows/uwp/get-started/universal-application-platform-guide UWP]
for Windows Phone development.

[table
    [[__predef_symbol__] [__predef_version__]]

    [[`WINAPI_FAMILY == WINAPI_FAMILY_PHONE_APP`] [__predef_detection__]]
    ]
 */

#define BOOST_PLAT_WINDOWS_PHONE BOOST_VERSION_NUMBER_NOT_AVAILABLE

#if BOOST_OS_WINDOWS && \
    defined(WINAPI_FAMILY_PHONE_APP) && WINAPI_FAMILY == WINAPI_FAMILY_PHONE_APP
#   undef BOOST_PLAT_WINDOWS_PHONE
#   define BOOST_PLAT_WINDOWS_PHONE BOOST_VERSION_NUMBER_AVAILABLE
#endif
 
#if BOOST_PLAT_WINDOWS_PHONE
#   define BOOST_PLAT_WINDOWS_PHONE_AVAILABLE
#   include <boost/predef/detail/platform_detected.h>
#endif

#define BOOST_PLAT_WINDOWS_PHONE_NAME "Windows Phone"

#endif

#include <boost/predef/detail/test.h>
BOOST_PREDEF_DECLARE_TEST(BOOST_PLAT_WINDOWS_PHONE,BOOST_PLAT_WINDOWS_PHONE_NAME)
