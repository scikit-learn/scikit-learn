/*
Copyright Rene Rivera 2008-2015
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE_1_0.txt or copy at
http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_PREDEF_LIBRARY_STD_MSL_H
#define BOOST_PREDEF_LIBRARY_STD_MSL_H

#include <boost/predef/library/std/_prefix.h>

#include <boost/predef/version_number.h>
#include <boost/predef/make.h>

/*`
[heading `BOOST_LIB_STD_MSL`]

[@http://www.freescale.com/ Metrowerks] Standard C++ Library.
Version number available as major, minor, and patch.

[table
    [[__predef_symbol__] [__predef_version__]]

    [[`__MSL_CPP__`] [__predef_detection__]]
    [[`__MSL__`] [__predef_detection__]]

    [[`__MSL_CPP__`] [V.R.P]]
    [[`__MSL__`] [V.R.P]]
    ]
 */

#define BOOST_LIB_STD_MSL BOOST_VERSION_NUMBER_NOT_AVAILABLE

#if defined(__MSL_CPP__) || defined(__MSL__)
#   undef BOOST_LIB_STD_MSL
#   if defined(__MSL_CPP__)
#       define BOOST_LIB_STD_MSL BOOST_PREDEF_MAKE_0X_VRPP(__MSL_CPP__)
#   else
#       define BOOST_LIB_STD_MSL BOOST_PREDEF_MAKE_0X_VRPP(__MSL__)
#   endif
#endif

#if BOOST_LIB_STD_MSL
#   define BOOST_LIB_STD_MSL_AVAILABLE
#endif

#define BOOST_LIB_STD_MSL_NAME "Metrowerks"

#endif

#include <boost/predef/detail/test.h>
BOOST_PREDEF_DECLARE_TEST(BOOST_LIB_STD_MSL,BOOST_LIB_STD_MSL_NAME)
