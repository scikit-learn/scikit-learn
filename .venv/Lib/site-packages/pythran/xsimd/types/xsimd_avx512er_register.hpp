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

#ifndef XSIMD_AVX512ER_REGISTER_HPP
#define XSIMD_AVX512ER_REGISTER_HPP

#include "./xsimd_avx512dq_register.hpp"

namespace xsimd
{

    /**
     * @ingroup architectures
     *
     * AVX512ER instructions
     */
    struct avx512er : avx512cd
    {
        static constexpr bool supported() noexcept { return XSIMD_WITH_AVX512ER; }
        static constexpr bool available() noexcept { return true; }
        static constexpr char const* name() noexcept { return "avx512er"; }
    };

#if XSIMD_WITH_AVX512ER

    namespace types
    {
        template <class T>
        struct get_bool_simd_register<T, avx512er>
        {
            using type = simd_avx512_bool_register<T>;
        };

        XSIMD_DECLARE_SIMD_REGISTER_ALIAS(avx512er, avx512cd);

    }
#endif
}
#endif
