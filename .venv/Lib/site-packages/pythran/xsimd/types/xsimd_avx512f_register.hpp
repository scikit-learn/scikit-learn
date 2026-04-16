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

#ifndef XSIMD_AVX512F_REGISTER_HPP
#define XSIMD_AVX512F_REGISTER_HPP

#include "./xsimd_generic_arch.hpp"

namespace xsimd
{

    /**
     * @ingroup architectures
     *
     * AVX512F instructions
     */
    struct avx512f : generic
    {
        static constexpr bool supported() noexcept { return XSIMD_WITH_AVX512F; }
        static constexpr bool available() noexcept { return true; }
        static constexpr std::size_t alignment() noexcept { return 64; }
        static constexpr bool requires_alignment() noexcept { return true; }
        static constexpr char const* name() noexcept { return "avx512f"; }
    };

#if XSIMD_WITH_AVX512F

    namespace types
    {
        template <class T>
        struct simd_avx512_bool_register
        {
            using register_type = typename std::conditional<
                (sizeof(T) < 4), std::conditional<(sizeof(T) == 1), __mmask64, __mmask32>,
                std::conditional<(sizeof(T) == 4), __mmask16, __mmask8>>::type::type;
            register_type data;
            simd_avx512_bool_register() = default;
            simd_avx512_bool_register(register_type r) { data = r; }
            operator register_type() const noexcept { return data; }
        };
        template <class T>
        struct get_bool_simd_register<T, avx512f>
        {
            using type = simd_avx512_bool_register<T>;
        };

        XSIMD_DECLARE_SIMD_REGISTER(signed char, avx512f, __m512i);
        XSIMD_DECLARE_SIMD_REGISTER(unsigned char, avx512f, __m512i);
        XSIMD_DECLARE_SIMD_REGISTER(char, avx512f, __m512i);
        XSIMD_DECLARE_SIMD_REGISTER(unsigned short, avx512f, __m512i);
        XSIMD_DECLARE_SIMD_REGISTER(short, avx512f, __m512i);
        XSIMD_DECLARE_SIMD_REGISTER(unsigned int, avx512f, __m512i);
        XSIMD_DECLARE_SIMD_REGISTER(int, avx512f, __m512i);
        XSIMD_DECLARE_SIMD_REGISTER(unsigned long int, avx512f, __m512i);
        XSIMD_DECLARE_SIMD_REGISTER(long int, avx512f, __m512i);
        XSIMD_DECLARE_SIMD_REGISTER(unsigned long long int, avx512f, __m512i);
        XSIMD_DECLARE_SIMD_REGISTER(long long int, avx512f, __m512i);
        XSIMD_DECLARE_SIMD_REGISTER(float, avx512f, __m512);
        XSIMD_DECLARE_SIMD_REGISTER(double, avx512f, __m512d);

    }
#endif
}

#endif
