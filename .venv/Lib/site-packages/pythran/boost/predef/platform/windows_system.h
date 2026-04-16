/*
Copyright James E. King III, 2017
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE_1_0.txt or copy at
http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_PREDEF_PLAT_WINDOWS_SYSTEM_H
#define BOOST_PREDEF_PLAT_WINDOWS_SYSTEM_H

#include <boost/predef/make.h>
#include <boost/predef/os/windows.h>
#include <boost/predef/platform/windows_uwp.h>
#include <boost/predef/version_number.h>

/*`
[heading `BOOST_PLAT_WINDOWS_SYSTEM`]

[@https://docs.microsoft.com/en-us/windows/uwp/get-started/universal-application-platform-guide UWP]
for Windows System development.

[table
    [[__predef_symbol__] [__predef_version__]]

    [[`WINAPI_FAMILY == WINAPI_FAMILY_SYSTEM`] [__predef_detection__]]
    ]
 */

#define BOOST_PLAT_WINDOWS_SYSTEM BOOST_VERSION_NUMBER_NOT_AVAILABLE

#if BOOST_OS_WINDOWS && \
    defined(WINAPI_FAMILY_SYSTEM) && WINAPI_FAMILY == WINAPI_FAMILY_SYSTEM
#   undef BOOST_PLAT_WINDOWS_SYSTEM
#   define BOOST_PLAT_WINDOWS_SYSTEM BOOST_VERSION_NUMBER_AVAILABLE
#endif
 
#if BOOST_PLAT_WINDOWS_SYSTEM
#   define BOOST_PLAT_WINDOWS_SYSTEM_AVAILABLE
#   include <boost/predef/detail/platform_detected.h>
#endif

#define BOOST_PLAT_WINDOWS_SYSTEM_NAME "Windows Drivers and Tools"

#endif

#include <boost/predef/detail/test.h>
BOOST_PREDEF_DECLARE_TEST(BOOST_PLAT_WINDOWS_SYSTEM,BOOST_PLAT_WINDOWS_SYSTEM_NAME)
