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

#ifndef XSIMD_SSE4_2_REGISTER_HPP
#define XSIMD_SSE4_2_REGISTER_HPP

#include "./xsimd_sse4_1_register.hpp"

#if XSIMD_WITH_SSE4_2
#include <nmmintrin.h>
#endif

namespace xsimd
{
    /**
     * @ingroup architectures
     *
     * SSE4.2 instructions
     */
    struct sse4_2 : sse4_1
    {
        static constexpr bool supported() noexcept { return XSIMD_WITH_SSE4_2; }
        static constexpr bool available() noexcept { return true; }
        static constexpr char const* name() noexcept { return "sse4.2"; }
    };

#if XSIMD_WITH_SSE4_2
    namespace types
    {
        XSIMD_DECLARE_SIMD_REGISTER_ALIAS(sse4_2, sse4_1);
    }
#endif
}

#endif
