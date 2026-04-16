/*
Copyright Rene Rivera 2008-2015
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE_1_0.txt or copy at
http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_PREDEF_LIBRARY_STD_DINKUMWARE_H
#define BOOST_PREDEF_LIBRARY_STD_DINKUMWARE_H

#include <boost/predef/library/std/_prefix.h>

#include <boost/predef/version_number.h>
#include <boost/predef/make.h>

/*`
[heading `BOOST_LIB_STD_DINKUMWARE`]

[@http://en.wikipedia.org/wiki/Dinkumware Dinkumware] Standard C++ Library.
If available version number as major, minor, and patch.

[table
    [[__predef_symbol__] [__predef_version__]]

    [[`_YVALS`, `__IBMCPP__`] [__predef_detection__]]
    [[`_CPPLIB_VER`] [__predef_detection__]]

    [[`_CPPLIB_VER`] [V.R.0]]
    ]
 */

#define BOOST_LIB_STD_DINKUMWARE BOOST_VERSION_NUMBER_NOT_AVAILABLE

#if (defined(_YVALS) && !defined(__IBMCPP__)) || defined(_CPPLIB_VER)
#   undef BOOST_LIB_STD_DINKUMWARE
#   if defined(_CPPLIB_VER)
#       define BOOST_LIB_STD_DINKUMWARE BOOST_PREDEF_MAKE_10_VVRR(_CPPLIB_VER)
#   else
#       define BOOST_LIB_STD_DINKUMWARE BOOST_VERSION_NUMBER_AVAILABLE
#   endif
#endif

#if BOOST_LIB_STD_DINKUMWARE
#   define BOOST_LIB_STD_DINKUMWARE_AVAILABLE
#endif

#define BOOST_LIB_STD_DINKUMWARE_NAME "Dinkumware"

#endif

#include <boost/predef/detail/test.h>
BOOST_PREDEF_DECLARE_TEST(BOOST_LIB_STD_DINKUMWARE,BOOST_LIB_STD_DINKUMWARE_NAME)
