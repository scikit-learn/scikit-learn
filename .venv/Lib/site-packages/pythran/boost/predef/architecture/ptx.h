/*
Copyright Benjamin Worpitz 2018
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE_1_0.txt or copy at
http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_PREDEF_ARCHITECTURE_PTX_H
#define BOOST_PREDEF_ARCHITECTURE_PTX_H

#include <boost/predef/version_number.h>
#include <boost/predef/make.h>

/*`
[heading `BOOST_ARCH_PTX`]

[@https://en.wikipedia.org/wiki/Parallel_Thread_Execution PTX] architecture.

[table
    [[__predef_symbol__] [__predef_version__]]

    [[`__CUDA_ARCH__`] [__predef_detection__]]

    [[`__CUDA_ARCH__`] [V.R.0]]
    ]
 */

#define BOOST_ARCH_PTX BOOST_VERSION_NUMBER_NOT_AVAILABLE

#if defined(__CUDA_ARCH__)
#   undef BOOST_ARCH_PTX
#   define BOOST_ARCH_PTX BOOST_PREDEF_MAKE_10_VR0(__CUDA_ARCH__)
#endif

#if BOOST_ARCH_PTX
#   define BOOST_ARCH_PTX_AVAILABLE
#endif

#define BOOST_ARCH_PTX_NAME "PTX"

#endif

#include <boost/predef/detail/test.h>
BOOST_PREDEF_DECLARE_TEST(BOOST_ARCH_PTX,BOOST_ARCH_PTX_NAME)
