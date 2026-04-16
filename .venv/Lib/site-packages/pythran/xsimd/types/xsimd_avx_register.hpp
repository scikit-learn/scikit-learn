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

#ifndef XSIMD_AVX_REGISTER_HPP
#define XSIMD_AVX_REGISTER_HPP

#include "./xsimd_generic_arch.hpp"

namespace xsimd
{

    /**
     * @ingroup architectures
     *
     * AVX instructions
     */
    struct avx : generic
    {
        static constexpr bool supported() noexcept { return XSIMD_WITH_AVX; }
        static constexpr bool available() noexcept { return true; }
        static constexpr std::size_t alignment() noexcept { return 32; }
        static constexpr bool requires_alignment() noexcept { return true; }
        static constexpr char const* name() noexcept { return "avx"; }
    };
}

#if XSIMD_WITH_AVX

#include <immintrin.h>

namespace xsimd
{
    namespace types
    {

        XSIMD_DECLARE_SIMD_REGISTER(signed char, avx, __m256i);
        XSIMD_DECLARE_SIMD_REGISTER(unsigned char, avx, __m256i);
        XSIMD_DECLARE_SIMD_REGISTER(char, avx, __m256i);
        XSIMD_DECLARE_SIMD_REGISTER(unsigned short, avx, __m256i);
        XSIMD_DECLARE_SIMD_REGISTER(short, avx, __m256i);
        XSIMD_DECLARE_SIMD_REGISTER(unsigned int, avx, __m256i);
        XSIMD_DECLARE_SIMD_REGISTER(int, avx, __m256i);
        XSIMD_DECLARE_SIMD_REGISTER(unsigned long int, avx, __m256i);
        XSIMD_DECLARE_SIMD_REGISTER(long int, avx, __m256i);
        XSIMD_DECLARE_SIMD_REGISTER(unsigned long long int, avx, __m256i);
        XSIMD_DECLARE_SIMD_REGISTER(long long int, avx, __m256i);
        XSIMD_DECLARE_SIMD_REGISTER(float, avx, __m256);
        XSIMD_DECLARE_SIMD_REGISTER(double, avx, __m256d);
    }
}
#endif
#endif
