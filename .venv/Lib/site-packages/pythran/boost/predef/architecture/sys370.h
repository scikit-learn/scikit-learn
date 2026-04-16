/*
Copyright Rene Rivera 2008-2015
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE_1_0.txt or copy at
http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_PREDEF_ARCHITECTURE_SYS370_H
#define BOOST_PREDEF_ARCHITECTURE_SYS370_H

#include <boost/predef/version_number.h>
#include <boost/predef/make.h>

/*`
[heading `BOOST_ARCH_SYS370`]

[@http://en.wikipedia.org/wiki/System/370 System/370] architecture.

[table
    [[__predef_symbol__] [__predef_version__]]

    [[`__370__`] [__predef_detection__]]
    [[`__THW_370__`] [__predef_detection__]]
    ]
 */

#define BOOST_ARCH_SYS370 BOOST_VERSION_NUMBER_NOT_AVAILABLE

#if defined(__370__) || defined(__THW_370__)
#   undef BOOST_ARCH_SYS370
#   define BOOST_ARCH_SYS370 BOOST_VERSION_NUMBER_AVAILABLE
#endif

#if BOOST_ARCH_SYS370
#   define BOOST_ARCH_SYS370_AVAILABLE
#endif

#define BOOST_ARCH_SYS370_NAME "System/370"

#endif

#include <boost/predef/detail/test.h>
BOOST_PREDEF_DECLARE_TEST(BOOST_ARCH_SYS370,BOOST_ARCH_SYS370_NAME)
