/*
Copyright Rene Rivera 2008-2015
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE_1_0.txt or copy at
http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_PREDEF_LIBRARY_C_VMS_H
#define BOOST_PREDEF_LIBRARY_C_VMS_H

#include <boost/predef/library/c/_prefix.h>

#include <boost/predef/version_number.h>
#include <boost/predef/make.h>

/*`
[heading `BOOST_LIB_C_VMS`]

VMS libc Standard C library.
Version number available as major, minor, and patch.

[table
    [[__predef_symbol__] [__predef_version__]]

    [[`__CRTL_VER`] [__predef_detection__]]

    [[`__CRTL_VER`] [V.R.P]]
    ]
 */

#define BOOST_LIB_C_VMS BOOST_VERSION_NUMBER_NOT_AVAILABLE

#if defined(__CRTL_VER)
#   undef BOOST_LIB_C_VMS
#   define BOOST_LIB_C_VMS BOOST_PREDEF_MAKE_10_VVRR0PP00(__CRTL_VER)
#endif

#if BOOST_LIB_C_VMS
#   define BOOST_LIB_C_VMS_AVAILABLE
#endif

#define BOOST_LIB_C_VMS_NAME "VMS"

#endif

#include <boost/predef/detail/test.h>
BOOST_PREDEF_DECLARE_TEST(BOOST_LIB_C_VMS,BOOST_LIB_C_VMS_NAME)
