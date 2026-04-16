/*
Copyright Rene Rivera 2008-2015
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE_1_0.txt or copy at
http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_PREDEF_COMPILER_DIAB_H
#define BOOST_PREDEF_COMPILER_DIAB_H

#include <boost/predef/version_number.h>
#include <boost/predef/make.h>

/*`
[heading `BOOST_COMP_DIAB`]

[@http://www.windriver.com/products/development_suite/wind_river_compiler/ Diab C/C++] compiler.
Version number available as major, minor, and patch.

[table
    [[__predef_symbol__] [__predef_version__]]

    [[`__DCC__`] [__predef_detection__]]

    [[`__VERSION_NUMBER__`] [V.R.P]]
    ]
 */

#define BOOST_COMP_DIAB BOOST_VERSION_NUMBER_NOT_AVAILABLE

#if defined(__DCC__)
#   define BOOST_COMP_DIAB_DETECTION BOOST_PREDEF_MAKE_10_VRPP(__VERSION_NUMBER__)
#endif

#ifdef BOOST_COMP_DIAB_DETECTION
#   if defined(BOOST_PREDEF_DETAIL_COMP_DETECTED)
#       define BOOST_COMP_DIAB_EMULATED BOOST_COMP_DIAB_DETECTION
#   else
#       undef BOOST_COMP_DIAB
#       define BOOST_COMP_DIAB BOOST_COMP_DIAB_DETECTION
#   endif
#   define BOOST_COMP_DIAB_AVAILABLE
#   include <boost/predef/detail/comp_detected.h>
#endif

#define BOOST_COMP_DIAB_NAME "Diab C/C++"

#endif

#include <boost/predef/detail/test.h>
BOOST_PREDEF_DECLARE_TEST(BOOST_COMP_DIAB,BOOST_COMP_DIAB_NAME)

#ifdef BOOST_COMP_DIAB_EMULATED
#include <boost/predef/detail/test.h>
BOOST_PREDEF_DECLARE_TEST(BOOST_COMP_DIAB_EMULATED,BOOST_COMP_DIAB_NAME)
#endif
