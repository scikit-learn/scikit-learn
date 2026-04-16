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

#ifndef XSIMD_FMA4_REGISTER_HPP
#define XSIMD_FMA4_REGISTER_HPP

#include "./xsimd_sse4_2_register.hpp"

namespace xsimd
{
    /**
     * @ingroup architectures
     *
     * SSE4.2 + FMA4 instructions
     */
    struct fma4 : sse4_2
    {
        static constexpr bool supported() noexcept { return XSIMD_WITH_FMA4; }
        static constexpr bool available() noexcept { return true; }
        static constexpr char const* name() noexcept { return "fma4"; }
    };

#if XSIMD_WITH_FMA4
    namespace types
    {

        XSIMD_DECLARE_SIMD_REGISTER_ALIAS(fma4, sse4_2);

    }
#endif

}
#endif
