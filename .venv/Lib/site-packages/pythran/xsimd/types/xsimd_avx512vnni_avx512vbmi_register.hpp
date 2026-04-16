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

#ifndef XSIMD_AVX512VNNI_AVX512VBMI_REGISTER_HPP
#define XSIMD_AVX512VNNI_AVX512VBMI_REGISTER_HPP

#include "./xsimd_avx512vbmi_register.hpp"

namespace xsimd
{
    template <typename arch>
    struct avx512vnni;

    /**
     * @ingroup architectures
     *
     * AVX512VNNI instructions
     */
    template <>
    struct avx512vnni<avx512vbmi> : avx512vbmi
    {
        static constexpr bool supported() noexcept { return XSIMD_WITH_AVX512VNNI_AVX512VBMI; }
        static constexpr bool available() noexcept { return true; }
        static constexpr char const* name() noexcept { return "avx512vnni+avx512vbmi"; }
    };

#if XSIMD_WITH_AVX512VNNI_AVX512VBMI

    namespace types
    {
        template <class T>
        struct get_bool_simd_register<T, avx512vnni<avx512vbmi>>
        {
            using type = simd_avx512_bool_register<T>;
        };

        XSIMD_DECLARE_SIMD_REGISTER_ALIAS(avx512vnni<avx512vbmi>, avx512vbmi);

    }
#endif
}
#endif
