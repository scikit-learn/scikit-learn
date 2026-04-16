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

#ifndef XSIMD_NEON64_REGISTER_HPP
#define XSIMD_NEON64_REGISTER_HPP

#include "xsimd_neon_register.hpp"

namespace xsimd
{
    /**
     * @ingroup architectures
     *
     * NEON instructions for arm64
     */
    struct neon64 : neon
    {
        static constexpr bool supported() noexcept { return XSIMD_WITH_NEON64; }
        static constexpr bool available() noexcept { return true; }
        static constexpr bool requires_alignment() noexcept { return true; }
        static constexpr std::size_t alignment() noexcept { return 16; }
        static constexpr char const* name() noexcept { return "arm64+neon"; }
    };

#if XSIMD_WITH_NEON64

    namespace types
    {
        XSIMD_DECLARE_SIMD_REGISTER_ALIAS(neon64, neon);
        XSIMD_DECLARE_SIMD_REGISTER(double, neon64, float64x2_t);

        template <class T>
        struct get_bool_simd_register<T, neon64>
            : detail::neon_bool_simd_register<T, neon64>
        {
        };
    }

#endif

}

#endif
