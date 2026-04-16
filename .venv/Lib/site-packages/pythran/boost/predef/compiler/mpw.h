/*
Copyright Rene Rivera 2008-2015
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE_1_0.txt or copy at
http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_PREDEF_COMPILER_MPW_H
#define BOOST_PREDEF_COMPILER_MPW_H

#include <boost/predef/version_number.h>
#include <boost/predef/make.h>

/*`
[heading `BOOST_COMP_MPW`]

[@http://en.wikipedia.org/wiki/Macintosh_Programmer%27s_Workshop MPW C++] compiler.
Version number available as major, and minor.

[table
    [[__predef_symbol__] [__predef_version__]]

    [[`__MRC__`] [__predef_detection__]]
    [[`MPW_C`] [__predef_detection__]]
    [[`MPW_CPLUS`] [__predef_detection__]]

    [[`__MRC__`] [V.R.0]]
    ]
 */

#define BOOST_COMP_MPW BOOST_VERSION_NUMBER_NOT_AVAILABLE

#if defined(__MRC__) || defined(MPW_C) || defined(MPW_CPLUS)
#   if !defined(BOOST_COMP_MPW_DETECTION) && defined(__MRC__)
#       define BOOST_COMP_MPW_DETECTION BOOST_PREDEF_MAKE_0X_VVRR(__MRC__)
#   endif
#   if !defined(BOOST_COMP_MPW_DETECTION)
#       define BOOST_COMP_MPW_DETECTION BOOST_VERSION_NUMBER_AVAILABLE
#   endif
#endif

#ifdef BOOST_COMP_MPW_DETECTION
#   if defined(BOOST_PREDEF_DETAIL_COMP_DETECTED)
#       define BOOST_COMP_MPW_EMULATED BOOST_COMP_MPW_DETECTION
#   else
#       undef BOOST_COMP_MPW
#       define BOOST_COMP_MPW BOOST_COMP_MPW_DETECTION
#   endif
#   define BOOST_COMP_MPW_AVAILABLE
#   include <boost/predef/detail/comp_detected.h>
#endif

#define BOOST_COMP_MPW_NAME "MPW C++"

#endif

#include <boost/predef/detail/test.h>
BOOST_PREDEF_DECLARE_TEST(BOOST_COMP_MPW,BOOST_COMP_MPW_NAME)

#ifdef BOOST_COMP_MPW_EMULATED
#include <boost/predef/detail/test.h>
BOOST_PREDEF_DECLARE_TEST(BOOST_COMP_MPW_EMULATED,BOOST_COMP_MPW_NAME)
#endif
