/*
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 * Copyright (c) 2017 Andrey Semashev
 */
/*!
 * \file   atomic/detail/hwcaps_gcc_ppc.hpp
 *
 * This header defines hardware capabilities macros for PowerPC
 */

#ifndef BOOST_ATOMIC_DETAIL_HWCAPS_GCC_PPC_HPP_INCLUDED_
#define BOOST_ATOMIC_DETAIL_HWCAPS_GCC_PPC_HPP_INCLUDED_

#include <boost/atomic/detail/config.hpp>

#ifdef BOOST_HAS_PRAGMA_ONCE
#pragma once
#endif

#if defined(__POWERPC__) || defined(__PPC__)

#if defined(_ARCH_PWR8)
// Power8 and later architectures have 8 and 16-bit instructions
#define BOOST_ATOMIC_DETAIL_PPC_HAS_LBARX_STBCX
#define BOOST_ATOMIC_DETAIL_PPC_HAS_LHARX_STHCX
#endif

#if defined(__powerpc64__) || defined(__PPC64__)
// Power7 and later architectures in 64-bit mode have 64-bit instructions
#define BOOST_ATOMIC_DETAIL_PPC_HAS_LDARX_STDCX
#if defined(_ARCH_PWR8)
// Power8 also has 128-bit instructions
#define BOOST_ATOMIC_DETAIL_PPC_HAS_LQARX_STQCX
#endif
#endif

#endif // defined(__POWERPC__) || defined(__PPC__)

#endif // BOOST_ATOMIC_DETAIL_HWCAPS_GCC_PPC_HPP_INCLUDED_
