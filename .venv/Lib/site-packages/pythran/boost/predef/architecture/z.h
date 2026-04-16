/*
Copyright Rene Rivera 2008-2015
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE_1_0.txt or copy at
http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_PREDEF_ARCHITECTURE_Z_H
#define BOOST_PREDEF_ARCHITECTURE_Z_H

#include <boost/predef/version_number.h>
#include <boost/predef/make.h>

/*`
[heading `BOOST_ARCH_Z`]

[@http://en.wikipedia.org/wiki/Z/Architecture z/Architecture] architecture.

[table
    [[__predef_symbol__] [__predef_version__]]

    [[`__SYSC_ZARCH__`] [__predef_detection__]]
    ]
 */

#define BOOST_ARCH_Z BOOST_VERSION_NUMBER_NOT_AVAILABLE

#if defined(__SYSC_ZARCH__)
#   undef BOOST_ARCH_Z
#   define BOOST_ARCH_Z BOOST_VERSION_NUMBER_AVAILABLE
#endif

#if BOOST_ARCH_Z
#   define BOOST_ARCH_Z_AVAILABLE
#endif

#define BOOST_ARCH_Z_NAME "z/Architecture"

#endif

#include <boost/predef/detail/test.h>
BOOST_PREDEF_DECLARE_TEST(BOOST_ARCH_Z,BOOST_ARCH_Z_NAME)
