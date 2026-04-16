/*
Copyright Rene Rivera 2008-2015
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE_1_0.txt or copy at
http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_PREDEF_LIBRARY_STD_STDCPP3_H
#define BOOST_PREDEF_LIBRARY_STD_STDCPP3_H

#include <boost/predef/library/std/_prefix.h>

#include <boost/predef/version_number.h>
#include <boost/predef/make.h>

/*`
[heading `BOOST_LIB_STD_GNU`]

[@http://gcc.gnu.org/libstdc++/ GNU libstdc++] Standard C++ library.
Version number available as year (from 1970), month, and day.

[table
    [[__predef_symbol__] [__predef_version__]]

    [[`__GLIBCXX__`] [__predef_detection__]]
    [[`__GLIBCPP__`] [__predef_detection__]]

    [[`__GLIBCXX__`] [V.R.P]]
    [[`__GLIBCPP__`] [V.R.P]]
    ]
 */

#define BOOST_LIB_STD_GNU BOOST_VERSION_NUMBER_NOT_AVAILABLE

#if defined(__GLIBCPP__) || defined(__GLIBCXX__)
#   undef BOOST_LIB_STD_GNU
#   if defined(__GLIBCXX__)
#       define BOOST_LIB_STD_GNU BOOST_PREDEF_MAKE_YYYYMMDD(__GLIBCXX__)
#   else
#       define BOOST_LIB_STD_GNU BOOST_PREDEF_MAKE_YYYYMMDD(__GLIBCPP__)
#   endif
#endif

#if BOOST_LIB_STD_GNU
#   define BOOST_LIB_STD_GNU_AVAILABLE
#endif

#define BOOST_LIB_STD_GNU_NAME "GNU"

#endif

#include <boost/predef/detail/test.h>
BOOST_PREDEF_DECLARE_TEST(BOOST_LIB_STD_GNU,BOOST_LIB_STD_GNU_NAME)
