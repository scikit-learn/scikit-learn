/*
Copyright Rene Rivera 2008-2015
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE_1_0.txt or copy at
http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_PREDEF_COMPILER_WATCOM_H
#define BOOST_PREDEF_COMPILER_WATCOM_H

#include <boost/predef/version_number.h>
#include <boost/predef/make.h>

/*`
[heading `BOOST_COMP_WATCOM`]

[@http://en.wikipedia.org/wiki/Watcom Watcom C++] compiler.
Version number available as major, and minor.

[table
    [[__predef_symbol__] [__predef_version__]]

    [[`__WATCOMC__`] [__predef_detection__]]

    [[`__WATCOMC__`] [V.R.P]]
    ]
 */

#define BOOST_COMP_WATCOM BOOST_VERSION_NUMBER_NOT_AVAILABLE

#if defined(__WATCOMC__)
#   define BOOST_COMP_WATCOM_DETECTION BOOST_PREDEF_MAKE_10_VVRR(__WATCOMC__)
#endif

#ifdef BOOST_COMP_WATCOM_DETECTION
#   if defined(BOOST_PREDEF_DETAIL_COMP_DETECTED)
#       define BOOST_COMP_WATCOM_EMULATED BOOST_COMP_WATCOM_DETECTION
#   else
#       undef BOOST_COMP_WATCOM
#       define BOOST_COMP_WATCOM BOOST_COMP_WATCOM_DETECTION
#   endif
#   define BOOST_COMP_WATCOM_AVAILABLE
#   include <boost/predef/detail/comp_detected.h>
#endif

#define BOOST_COMP_WATCOM_NAME "Watcom C++"

#endif

#include <boost/predef/detail/test.h>
BOOST_PREDEF_DECLARE_TEST(BOOST_COMP_WATCOM,BOOST_COMP_WATCOM_NAME)

#ifdef BOOST_COMP_WATCOM_EMULATED
#include <boost/predef/detail/test.h>
BOOST_PREDEF_DECLARE_TEST(BOOST_COMP_WATCOM_EMULATED,BOOST_COMP_WATCOM_NAME)
#endif
