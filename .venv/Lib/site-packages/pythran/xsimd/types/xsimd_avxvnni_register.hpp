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

#ifndef XSIMD_AVXVNNI_REGISTER_HPP
#define XSIMD_AVXVNNI_REGISTER_HPP

#include "./xsimd_avx2_register.hpp"

namespace xsimd
{
    /**
     * @ingroup architectures
     *
     * AVXVNNI instructions
     */
    struct avxvnni : avx2
    {
        static constexpr bool supported() noexcept { return XSIMD_WITH_AVXVNNI; }
        static constexpr bool available() noexcept { return true; }
        static constexpr char const* name() noexcept { return "avxvnni"; }
    };

#if XSIMD_WITH_AVXVNNI
    namespace types
    {
        XSIMD_DECLARE_SIMD_REGISTER_ALIAS(avxvnni, avx2);
    }
#endif
}

#endif
