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

#ifndef XSIMD_FMA3_AVX2_REGISTER_HPP
#define XSIMD_FMA3_AVX2_REGISTER_HPP

#include "./xsimd_avx2_register.hpp"

namespace xsimd
{
    template <typename arch>
    struct fma3;

    /**
     * @ingroup architectures
     *
     * AVX2 + FMA instructions
     */
    template <>
    struct fma3<avx2> : avx2
    {
        static constexpr bool supported() noexcept { return XSIMD_WITH_FMA3_AVX2; }
        static constexpr bool available() noexcept { return true; }
        static constexpr char const* name() noexcept { return "fma3+avx2"; }
    };

#if XSIMD_WITH_FMA3_AVX2
    namespace types
    {

        XSIMD_DECLARE_SIMD_REGISTER_ALIAS(fma3<avx2>, avx2);

    }
#endif

}
#endif
