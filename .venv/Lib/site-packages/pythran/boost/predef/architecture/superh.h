/*
Copyright Rene Rivera 2008-2015
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE_1_0.txt or copy at
http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_PREDEF_ARCHITECTURE_SUPERH_H
#define BOOST_PREDEF_ARCHITECTURE_SUPERH_H

#include <boost/predef/version_number.h>
#include <boost/predef/make.h>

/*`
[heading `BOOST_ARCH_SH`]

[@http://en.wikipedia.org/wiki/SuperH SuperH] architecture:
If available versions \[1-5\] are specifically detected.

[table
    [[__predef_symbol__] [__predef_version__]]

    [[`__sh__`] [__predef_detection__]]

    [[`__SH5__`] [5.0.0]]
    [[`__SH4__`] [4.0.0]]
    [[`__sh3__`] [3.0.0]]
    [[`__SH3__`] [3.0.0]]
    [[`__sh2__`] [2.0.0]]
    [[`__sh1__`] [1.0.0]]
    ]
 */

#define BOOST_ARCH_SH BOOST_VERSION_NUMBER_NOT_AVAILABLE

#if defined(__sh__)
#   undef BOOST_ARCH_SH
#   if !defined(BOOST_ARCH_SH) && (defined(__SH5__))
#       define BOOST_ARCH_SH BOOST_VERSION_NUMBER(5,0,0)
#   endif
#   if !defined(BOOST_ARCH_SH) && (defined(__SH4__))
#       define BOOST_ARCH_SH BOOST_VERSION_NUMBER(4,0,0)
#   endif
#   if !defined(BOOST_ARCH_SH) && (defined(__sh3__) || defined(__SH3__))
#       define BOOST_ARCH_SH BOOST_VERSION_NUMBER(3,0,0)
#   endif
#   if !defined(BOOST_ARCH_SH) && (defined(__sh2__))
#       define BOOST_ARCH_SH BOOST_VERSION_NUMBER(2,0,0)
#   endif
#   if !defined(BOOST_ARCH_SH) && (defined(__sh1__))
#       define BOOST_ARCH_SH BOOST_VERSION_NUMBER(1,0,0)
#   endif
#   if !defined(BOOST_ARCH_SH)
#       define BOOST_ARCH_SH BOOST_VERSION_NUMBER_AVAILABLE
#   endif
#endif

#if BOOST_ARCH_SH
#   define BOOST_ARCH_SH_AVAILABLE
#endif

#define BOOST_ARCH_SH_NAME "SuperH"

#endif

#include <boost/predef/detail/test.h>
BOOST_PREDEF_DECLARE_TEST(BOOST_ARCH_SH,BOOST_ARCH_SH_NAME)
