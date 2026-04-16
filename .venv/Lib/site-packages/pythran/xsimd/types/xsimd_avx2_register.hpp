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

#ifndef XSIMD_AVX2_REGISTER_HPP
#define XSIMD_AVX2_REGISTER_HPP

#include "./xsimd_avx_register.hpp"

namespace xsimd
{
    /**
     * @ingroup architectures
     *
     * AVX2 instructions
     */
    struct avx2 : avx
    {
        static constexpr bool supported() noexcept { return XSIMD_WITH_AVX2; }
        static constexpr bool available() noexcept { return true; }
        static constexpr char const* name() noexcept { return "avx2"; }
    };

#if XSIMD_WITH_AVX2
    namespace types
    {
        XSIMD_DECLARE_SIMD_REGISTER_ALIAS(avx2, avx);
    }
#endif
}

#endif
