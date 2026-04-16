/*
Copyright Rene Rivera 2008-2019
Copyright Franz Detro 2014
Copyright (c) Microsoft Corporation 2014
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE_1_0.txt or copy at
http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_PREDEF_ARCHITECTURE_ARM_H
#define BOOST_PREDEF_ARCHITECTURE_ARM_H

#include <boost/predef/version_number.h>
#include <boost/predef/make.h>

/*`
[heading `BOOST_ARCH_ARM`]

[@http://en.wikipedia.org/wiki/ARM_architecture ARM] architecture.

[table
    [[__predef_symbol__] [__predef_version__]]

    [[`__ARM_ARCH`] [__predef_detection__]]
    [[`__TARGET_ARCH_ARM`] [__predef_detection__]]
    [[`__TARGET_ARCH_THUMB`] [__predef_detection__]]
    [[`_M_ARM`] [__predef_detection__]]
    [[`__arm__`] [__predef_detection__]]
    [[`__arm64`] [__predef_detection__]]
    [[`__thumb__`] [__predef_detection__]]
    [[`_M_ARM64`] [__predef_detection__]]
    [[`__aarch64__`] [__predef_detection__]]
    [[`__AARCH64EL__`] [__predef_detection__]]
    [[`__ARM_ARCH_7__`] [__predef_detection__]]
    [[`__ARM_ARCH_7A__`] [__predef_detection__]]
    [[`__ARM_ARCH_7R__`] [__predef_detection__]]
    [[`__ARM_ARCH_7M__`] [__predef_detection__]]
    [[`__ARM_ARCH_6K__`] [__predef_detection__]]
    [[`__ARM_ARCH_6Z__`] [__predef_detection__]]
    [[`__ARM_ARCH_6KZ__`] [__predef_detection__]]
    [[`__ARM_ARCH_6T2__`] [__predef_detection__]]
    [[`__ARM_ARCH_5TE__`] [__predef_detection__]]
    [[`__ARM_ARCH_5TEJ__`] [__predef_detection__]]
    [[`__ARM_ARCH_4T__`] [__predef_detection__]]
    [[`__ARM_ARCH_4__`] [__predef_detection__]]

    [[`__ARM_ARCH`] [V.0.0]]
    [[`__TARGET_ARCH_ARM`] [V.0.0]]
    [[`__TARGET_ARCH_THUMB`] [V.0.0]]
    [[`_M_ARM`] [V.0.0]]
    [[`__arm64`] [8.0.0]]
    [[`_M_ARM64`] [8.0.0]]
    [[`__aarch64__`] [8.0.0]]
    [[`__AARCH64EL__`] [8.0.0]]
    [[`__ARM_ARCH_7__`] [7.0.0]]
    [[`__ARM_ARCH_7A__`] [7.0.0]]
    [[`__ARM_ARCH_7R__`] [7.0.0]]
    [[`__ARM_ARCH_7M__`] [7.0.0]]
    [[`__ARM_ARCH_6K__`] [6.0.0]]
    [[`__ARM_ARCH_6Z__`] [6.0.0]]
    [[`__ARM_ARCH_6KZ__`] [6.0.0]]
    [[`__ARM_ARCH_6T2__`] [6.0.0]]
    [[`__ARM_ARCH_5TE__`] [5.0.0]]
    [[`__ARM_ARCH_5TEJ__`] [5.0.0]]
    [[`__ARM_ARCH_4T__`] [4.0.0]]
    [[`__ARM_ARCH_4__`] [4.0.0]]
    ]
 */

#define BOOST_ARCH_ARM BOOST_VERSION_NUMBER_NOT_AVAILABLE

#if \
    defined(__ARM_ARCH) || defined(__TARGET_ARCH_ARM) || \
    defined(__TARGET_ARCH_THUMB) || defined(_M_ARM) || \
    defined(__arm__) || defined(__arm64) || defined(__thumb__) || \
    defined(_M_ARM64) || defined(__aarch64__) || defined(__AARCH64EL__) || \
    defined(__ARM_ARCH_7__) || defined(__ARM_ARCH_7A__) || \
    defined(__ARM_ARCH_7R__) || defined(__ARM_ARCH_7M__) || \
    defined(__ARM_ARCH_6K__) || defined(__ARM_ARCH_6Z__) || \
    defined(__ARM_ARCH_6KZ__) || defined(__ARM_ARCH_6T2__) || \
    defined(__ARM_ARCH_5TE__) || defined(__ARM_ARCH_5TEJ__) || \
    defined(__ARM_ARCH_4T__) || defined(__ARM_ARCH_4__)
#   undef BOOST_ARCH_ARM
#   if !defined(BOOST_ARCH_ARM) && defined(__ARM_ARCH)
#       define BOOST_ARCH_ARM BOOST_VERSION_NUMBER(__ARM_ARCH,0,0)
#   endif
#   if !defined(BOOST_ARCH_ARM) && defined(__TARGET_ARCH_ARM)
#       define BOOST_ARCH_ARM BOOST_VERSION_NUMBER(__TARGET_ARCH_ARM,0,0)
#   endif
#   if !defined(BOOST_ARCH_ARM) && defined(__TARGET_ARCH_THUMB)
#       define BOOST_ARCH_ARM BOOST_VERSION_NUMBER(__TARGET_ARCH_THUMB,0,0)
#   endif
#   if !defined(BOOST_ARCH_ARM) && defined(_M_ARM)
#       define BOOST_ARCH_ARM BOOST_VERSION_NUMBER(_M_ARM,0,0)
#   endif
#   if !defined(BOOST_ARCH_ARM) && ( \
        defined(__arm64) || defined(_M_ARM64) || defined(__aarch64__) || \
        defined(__AARCH64EL__) )
#       define BOOST_ARCH_ARM BOOST_VERSION_NUMBER(8,0,0)
#   endif
#   if !defined(BOOST_ARCH_ARM) && ( \
    defined(__ARM_ARCH_7__) || defined(__ARM_ARCH_7A__) || \
    defined(__ARM_ARCH_7R__) || defined(__ARM_ARCH_7M__) )
#       define BOOST_ARCH_ARM BOOST_VERSION_NUMBER(7,0,0)
#   endif
#   if !defined(BOOST_ARCH_ARM) && ( \
    defined(__ARM_ARCH_6K__) || defined(__ARM_ARCH_6Z__) || \
    defined(__ARM_ARCH_6KZ__) || defined(__ARM_ARCH_6T2__) )
#       define BOOST_ARCH_ARM BOOST_VERSION_NUMBER(6,0,0)
#   endif
#   if !defined(BOOST_ARCH_ARM) && ( \
    defined(__ARM_ARCH_5TE__) || defined(__ARM_ARCH_5TEJ__) )
#       define BOOST_ARCH_ARM BOOST_VERSION_NUMBER(5,0,0)
#   endif
#   if !defined(BOOST_ARCH_ARM) && ( \
    defined(__ARM_ARCH_4T__) || defined(__ARM_ARCH_4__) )
#       define BOOST_ARCH_ARM BOOST_VERSION_NUMBER(4,0,0)
#   endif
#   if !defined(BOOST_ARCH_ARM)
#       define BOOST_ARCH_ARM BOOST_VERSION_NUMBER_AVAILABLE
#   endif
#endif

#if BOOST_ARCH_ARM
#   define BOOST_ARCH_ARM_AVAILABLE
#endif

#define BOOST_ARCH_ARM_NAME "ARM"

#endif

#include <boost/predef/detail/test.h>
BOOST_PREDEF_DECLARE_TEST(BOOST_ARCH_ARM,BOOST_ARCH_ARM_NAME)
