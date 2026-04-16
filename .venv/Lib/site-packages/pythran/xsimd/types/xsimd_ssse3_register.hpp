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

#ifndef XSIMD_SSSE3_REGISTER_HPP
#define XSIMD_SSSE3_REGISTER_HPP

#include "./xsimd_sse3_register.hpp"

#if XSIMD_WITH_SSSE3
#include <tmmintrin.h>
#endif

namespace xsimd
{
    /**
     * @ingroup architectures
     *
     * SSSE3 instructions
     */
    struct ssse3 : sse3
    {
        static constexpr bool supported() noexcept { return XSIMD_WITH_SSSE3; }
        static constexpr bool available() noexcept { return true; }
        static constexpr char const* name() noexcept { return "ssse3"; }
    };

#if XSIMD_WITH_SSSE3
    namespace types
    {
        XSIMD_DECLARE_SIMD_REGISTER_ALIAS(ssse3, sse3);
    }
#endif
}

#endif
