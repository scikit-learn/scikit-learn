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

#ifndef XSIMD_AVX512IFMA_REGISTER_HPP
#define XSIMD_AVX512IFMA_REGISTER_HPP

#include "./xsimd_avx512bw_register.hpp"

namespace xsimd
{

    /**
     * @ingroup architectures
     *
     * AVX512IFMA instructions
     */
    struct avx512ifma : avx512bw
    {
        static constexpr bool supported() noexcept { return XSIMD_WITH_AVX512IFMA; }
        static constexpr bool available() noexcept { return true; }
        static constexpr char const* name() noexcept { return "avx512ifma"; }
    };

#if XSIMD_WITH_AVX512IFMA

    namespace types
    {
        template <class T>
        struct get_bool_simd_register<T, avx512ifma>
        {
            using type = simd_avx512_bool_register<T>;
        };

        XSIMD_DECLARE_SIMD_REGISTER_ALIAS(avx512ifma, avx512bw);

    }
#endif
}
#endif
