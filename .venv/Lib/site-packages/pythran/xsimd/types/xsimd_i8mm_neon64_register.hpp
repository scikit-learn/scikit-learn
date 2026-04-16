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

#ifndef XSIMD_I8MM_NEON64_REGISTER_HPP
#define XSIMD_I8MM_NEON64_REGISTER_HPP

#include "./xsimd_neon64_register.hpp"

namespace xsimd
{
    template <typename arch>
    struct i8mm;

    /**
     * @ingroup architectures
     *
     * Neon64 + i8mm instructions
     */
    template <>
    struct i8mm<neon64> : neon64
    {
        static constexpr bool supported() noexcept { return XSIMD_WITH_I8MM_NEON64; }
        static constexpr bool available() noexcept { return true; }
        static constexpr char const* name() noexcept { return "i8mm+neon64"; }
    };

#if XSIMD_WITH_I8MM_NEON64
    namespace types
    {

        XSIMD_DECLARE_SIMD_REGISTER_ALIAS(i8mm<neon64>, neon64);

        template <class T>
        struct get_bool_simd_register<T, i8mm<neon64>>
            : detail::neon_bool_simd_register<T, i8mm<neon64>>
        {
        };
    }
#endif

}
#endif
