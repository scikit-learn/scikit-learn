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

#ifndef XSIMD_SSE2_REGISTER_HPP
#define XSIMD_SSE2_REGISTER_HPP

#include "./xsimd_generic_arch.hpp"
#include "./xsimd_register.hpp"

#if XSIMD_WITH_SSE2
#include <emmintrin.h>
#include <xmmintrin.h>
#endif

namespace xsimd
{
    /**
     * @ingroup architectures
     *
     * SSE2 instructions
     */
    struct sse2 : generic
    {
        static constexpr bool supported() noexcept { return XSIMD_WITH_SSE2; }
        static constexpr bool available() noexcept { return true; }
        static constexpr bool requires_alignment() noexcept { return true; }
        static constexpr std::size_t alignment() noexcept { return 16; }
        static constexpr char const* name() noexcept { return "sse2"; }
    };

#if XSIMD_WITH_SSE2
    namespace types
    {
        XSIMD_DECLARE_SIMD_REGISTER(signed char, sse2, __m128i);
        XSIMD_DECLARE_SIMD_REGISTER(unsigned char, sse2, __m128i);
        XSIMD_DECLARE_SIMD_REGISTER(char, sse2, __m128i);
        XSIMD_DECLARE_SIMD_REGISTER(unsigned short, sse2, __m128i);
        XSIMD_DECLARE_SIMD_REGISTER(short, sse2, __m128i);
        XSIMD_DECLARE_SIMD_REGISTER(unsigned int, sse2, __m128i);
        XSIMD_DECLARE_SIMD_REGISTER(int, sse2, __m128i);
        XSIMD_DECLARE_SIMD_REGISTER(unsigned long int, sse2, __m128i);
        XSIMD_DECLARE_SIMD_REGISTER(long int, sse2, __m128i);
        XSIMD_DECLARE_SIMD_REGISTER(unsigned long long int, sse2, __m128i);
        XSIMD_DECLARE_SIMD_REGISTER(long long int, sse2, __m128i);
        XSIMD_DECLARE_SIMD_REGISTER(float, sse2, __m128);
        XSIMD_DECLARE_SIMD_REGISTER(double, sse2, __m128d);
    }
#endif
}

#endif
