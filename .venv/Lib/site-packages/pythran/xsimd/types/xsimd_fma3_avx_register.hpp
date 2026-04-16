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

#ifndef XSIMD_FMA3_AVX_REGISTER_HPP
#define XSIMD_FMA3_AVX_REGISTER_HPP

#include "./xsimd_avx_register.hpp"

namespace xsimd
{
    template <typename arch>
    struct fma3;

    /**
     * @ingroup architectures
     *
     * AVX + FMA instructions
     */
    template <>
    struct fma3<avx> : avx
    {
        static constexpr bool supported() noexcept { return XSIMD_WITH_FMA3_AVX; }
        static constexpr bool available() noexcept { return true; }
        static constexpr char const* name() noexcept { return "fma3+avx"; }
    };

#if XSIMD_WITH_FMA3_AVX
    namespace types
    {

        XSIMD_DECLARE_SIMD_REGISTER_ALIAS(fma3<avx>, avx);

    }
#endif

}
#endif
