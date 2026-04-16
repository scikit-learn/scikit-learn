/*
Copyright Rene Rivera 2008-2015
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE_1_0.txt or copy at
http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_PREDEF_ARCHITECTURE_ALPHA_H
#define BOOST_PREDEF_ARCHITECTURE_ALPHA_H

#include <boost/predef/version_number.h>
#include <boost/predef/make.h>

/*`
[heading `BOOST_ARCH_ALPHA`]

[@http://en.wikipedia.org/wiki/DEC_Alpha DEC Alpha] architecture.

[table
    [[__predef_symbol__] [__predef_version__]]
    [[`__alpha__`] [__predef_detection__]]
    [[`__alpha`] [__predef_detection__]]
    [[`_M_ALPHA`] [__predef_detection__]]

    [[`__alpha_ev4__`] [4.0.0]]
    [[`__alpha_ev5__`] [5.0.0]]
    [[`__alpha_ev6__`] [6.0.0]]
    ]
 */

#define BOOST_ARCH_ALPHA BOOST_VERSION_NUMBER_NOT_AVAILABLE

#if defined(__alpha__) || defined(__alpha) || \
    defined(_M_ALPHA)
#   undef BOOST_ARCH_ALPHA
#   if !defined(BOOST_ARCH_ALPHA) && defined(__alpha_ev4__)
#       define BOOST_ARCH_ALPHA BOOST_VERSION_NUMBER(4,0,0)
#   endif
#   if !defined(BOOST_ARCH_ALPHA) && defined(__alpha_ev5__)
#       define BOOST_ARCH_ALPHA BOOST_VERSION_NUMBER(5,0,0)
#   endif
#   if !defined(BOOST_ARCH_ALPHA) && defined(__alpha_ev6__)
#       define BOOST_ARCH_ALPHA BOOST_VERSION_NUMBER(6,0,0)
#   endif
#   if !defined(BOOST_ARCH_ALPHA)
#       define BOOST_ARCH_ALPHA BOOST_VERSION_NUMBER_AVAILABLE
#   endif
#endif

#if BOOST_ARCH_ALPHA
#   define BOOST_ARCH_ALPHA_AVAILABLE
#endif

#define BOOST_ARCH_ALPHA_NAME "DEC Alpha"

#endif

#include <boost/predef/detail/test.h>
BOOST_PREDEF_DECLARE_TEST(BOOST_ARCH_ALPHA,BOOST_ARCH_ALPHA_NAME)
