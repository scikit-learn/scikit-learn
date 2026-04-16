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

#ifndef XSIMD_NEON_REGISTER_HPP
#define XSIMD_NEON_REGISTER_HPP

#include "xsimd_generic_arch.hpp"
#include "xsimd_register.hpp"

#if XSIMD_WITH_NEON
#include <arm_neon.h>
#endif

namespace xsimd
{
    /**
     * @ingroup architectures
     *
     * NEON instructions for arm32
     */
    struct neon : generic
    {
        static constexpr bool supported() noexcept { return XSIMD_WITH_NEON; }
        static constexpr bool available() noexcept { return true; }
        static constexpr bool requires_alignment() noexcept { return true; }
        static constexpr std::size_t alignment() noexcept { return 16; }
        static constexpr char const* name() noexcept { return "arm32+neon"; }
    };

#if XSIMD_WITH_NEON
    namespace types
    {
        namespace detail
        {
            template <size_t S>
            struct neon_vector_type_impl;

            template <>
            struct neon_vector_type_impl<8>
            {
                using signed_type = int8x16_t;
                using unsigned_type = uint8x16_t;
            };

            template <>
            struct neon_vector_type_impl<16>
            {
                using signed_type = int16x8_t;
                using unsigned_type = uint16x8_t;
            };

            template <>
            struct neon_vector_type_impl<32>
            {
                using signed_type = int32x4_t;
                using unsigned_type = uint32x4_t;
            };

            template <>
            struct neon_vector_type_impl<64>
            {
                using signed_type = int64x2_t;
                using unsigned_type = uint64x2_t;
            };

            template <class T>
            using signed_neon_vector_type = typename neon_vector_type_impl<8 * sizeof(T)>::signed_type;

            template <class T>
            using unsigned_neon_vector_type = typename neon_vector_type_impl<8 * sizeof(T)>::unsigned_type;

            template <class T>
            using neon_vector_type = typename std::conditional<std::is_signed<T>::value,
                                                               signed_neon_vector_type<T>,
                                                               unsigned_neon_vector_type<T>>::type;

            using char_neon_vector_type = typename std::conditional<std::is_signed<char>::value,
                                                                    signed_neon_vector_type<char>,
                                                                    unsigned_neon_vector_type<char>>::type;
        }

        XSIMD_DECLARE_SIMD_REGISTER(signed char, neon, detail::neon_vector_type<signed char>);
        XSIMD_DECLARE_SIMD_REGISTER(unsigned char, neon, detail::neon_vector_type<unsigned char>);
        XSIMD_DECLARE_SIMD_REGISTER(char, neon, detail::char_neon_vector_type);
        XSIMD_DECLARE_SIMD_REGISTER(short, neon, detail::neon_vector_type<short>);
        XSIMD_DECLARE_SIMD_REGISTER(unsigned short, neon, detail::neon_vector_type<unsigned short>);
        XSIMD_DECLARE_SIMD_REGISTER(int, neon, detail::neon_vector_type<int>);
        XSIMD_DECLARE_SIMD_REGISTER(unsigned int, neon, detail::neon_vector_type<unsigned int>);
        XSIMD_DECLARE_SIMD_REGISTER(long int, neon, detail::neon_vector_type<long int>);
        XSIMD_DECLARE_SIMD_REGISTER(unsigned long int, neon, detail::neon_vector_type<unsigned long int>);
        XSIMD_DECLARE_SIMD_REGISTER(long long int, neon, detail::neon_vector_type<long long int>);
        XSIMD_DECLARE_SIMD_REGISTER(unsigned long long int, neon, detail::neon_vector_type<unsigned long long int>);
        XSIMD_DECLARE_SIMD_REGISTER(float, neon, float32x4_t);
        XSIMD_DECLARE_INVALID_SIMD_REGISTER(double, neon);

        namespace detail
        {
            template <size_t S>
            struct get_unsigned_type;

            template <>
            struct get_unsigned_type<1>
            {
                using type = uint8_t;
            };

            template <>
            struct get_unsigned_type<2>
            {
                using type = uint16_t;
            };

            template <>
            struct get_unsigned_type<4>
            {
                using type = uint32_t;
            };

            template <>
            struct get_unsigned_type<8>
            {
                using type = uint64_t;
            };

            template <size_t S>
            using get_unsigned_type_t = typename get_unsigned_type<S>::type;

            template <class T, class A>
            struct neon_bool_simd_register
            {
                using type = simd_register<get_unsigned_type_t<sizeof(T)>, A>;
            };
        }

        template <class T>
        struct get_bool_simd_register<T, neon>
            : detail::neon_bool_simd_register<T, neon>
        {
        };

    }
#endif

}

#endif
