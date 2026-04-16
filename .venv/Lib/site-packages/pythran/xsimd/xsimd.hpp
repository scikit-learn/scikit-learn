/***************************************************************************
 * Copyright (c) Johan Mabille, Sylvain Corlay, Wolf Vollprecht and         *
 * Martin Renou                                                             *
 * Copyright (c) QuantStack                                                 *
 * Copyright (c) Serge Guelton                                              *
 *                                                                          *
 * Distributed under the terms of the BSD 3-Clause License.                 *
 *                                                                          *
 * The full license is in the file LICENSE, distributed with this software. *
 ****************************************************************************/

#ifndef XSIMD_HPP
#define XSIMD_HPP

#if defined(__has_cpp_attribute)
// if this check passes, then the compiler supports feature test macros
#if __has_cpp_attribute(nodiscard) >= 201603L
// if this check passes, then the compiler supports [[nodiscard]] without a message
#define XSIMD_NO_DISCARD [[nodiscard]]
#endif
#endif

#if !defined(XSIMD_NO_DISCARD) && __cplusplus >= 201703L
// this means that the previous tests failed, but we are using C++17 or higher
#define XSIMD_NO_DISCARD [[nodiscard]]
#endif

#if !defined(XSIMD_NO_DISCARD) && (defined(__GNUC__) || defined(__clang__))
// this means that the previous checks failed, but we are using GCC or Clang
#define XSIMD_NO_DISCARD __attribute__((warn_unused_result))
#endif

#if !defined(XSIMD_NO_DISCARD)
// this means that all the previous checks failed, so we fallback to doing nothing
#define XSIMD_NO_DISCARD
#endif

#ifdef __cpp_if_constexpr
// this means that the compiler supports the `if constexpr` construct
#define XSIMD_IF_CONSTEXPR if constexpr
#endif

#if !defined(XSIMD_IF_CONSTEXPR) && __cplusplus >= 201703L
// this means that the previous test failed, but we are using C++17 or higher
#define XSIMD_IF_CONSTEXPR if constexpr
#endif

#if !defined(XSIMD_IF_CONSTEXPR)
// this means that all the previous checks failed, so we fallback to a normal `if`
#define XSIMD_IF_CONSTEXPR if
#endif

#include "config/xsimd_config.hpp"
#include "config/xsimd_inline.hpp"

#include "arch/xsimd_scalar.hpp"
#include "memory/xsimd_aligned_allocator.hpp"

#if defined(XSIMD_NO_SUPPORTED_ARCHITECTURE)
// to type definition or anything appart from scalar definition and aligned allocator
#else
#include "types/xsimd_batch.hpp"
#include "types/xsimd_batch_constant.hpp"
#include "types/xsimd_traits.hpp"

// This include must come last
#include "types/xsimd_api.hpp"
#endif
#endif
