/*
Copyright Rene Rivera 2008-2015
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE_1_0.txt or copy at
http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_PREDEF_LIBRARY_STD_LIBCOMO_H
#define BOOST_PREDEF_LIBRARY_STD_LIBCOMO_H

#include <boost/predef/library/std/_prefix.h>

#include <boost/predef/version_number.h>
#include <boost/predef/make.h>

/*`
[heading `BOOST_LIB_STD_COMO`]

[@http://www.comeaucomputing.com/libcomo/ Comeau Computing] Standard C++ Library.
Version number available as major.

[table
    [[__predef_symbol__] [__predef_version__]]

    [[`__LIBCOMO__`] [__predef_detection__]]

    [[`__LIBCOMO_VERSION__`] [V.0.0]]
    ]
 */

#define BOOST_LIB_STD_COMO BOOST_VERSION_NUMBER_NOT_AVAILABLE

#if defined(__LIBCOMO__)
#   undef BOOST_LIB_STD_COMO
#   define BOOST_LIB_STD_COMO BOOST_VERSION_NUMBER(__LIBCOMO_VERSION__,0,0)
#endif

#if BOOST_LIB_STD_COMO
#   define BOOST_LIB_STD_COMO_AVAILABLE
#endif

#define BOOST_LIB_STD_COMO_NAME "Comeau Computing"

#endif

#include <boost/predef/detail/test.h>
BOOST_PREDEF_DECLARE_TEST(BOOST_LIB_STD_COMO,BOOST_LIB_STD_COMO_NAME)
