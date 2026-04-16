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

#ifndef XSIMD_SSE3_REGISTER_HPP
#define XSIMD_SSE3_REGISTER_HPP

#include "./xsimd_sse2_register.hpp"

#if XSIMD_WITH_SSE3
#include <pmmintrin.h>
#endif

namespace xsimd
{
    /**
     * @ingroup architectures
     *
     * SSE3 instructions
     */
    struct sse3 : sse2
    {
        static constexpr bool supported() noexcept { return XSIMD_WITH_SSE3; }
        static constexpr bool available() noexcept { return true; }
        static constexpr char const* name() noexcept { return "sse3"; }
    };

#if XSIMD_WITH_SSE3
    namespace types
    {

        XSIMD_DECLARE_SIMD_REGISTER_ALIAS(sse3, sse2);
    }
#endif
}

#endif
