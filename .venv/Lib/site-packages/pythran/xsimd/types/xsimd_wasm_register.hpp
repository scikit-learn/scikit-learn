/***************************************************************************
 * Copyright (c) Johan Mabille, Sylvain Corlay, Wolf Vollprecht and         *
 * Martin Renou                                                             *
 * Copyright (c) QuantStack                                                 *
 * Copyright (c) Serge Guelton                                              *
 * Copyright (c) Anutosh Bhat                                               *
 *                                                                          *
 * Distributed under the terms of the BSD 3-Clause License.                 *
 *                                                                          *
 * The full license is in the file LICENSE, distributed with this software. *
 ****************************************************************************/

#ifndef XSIMD_WASM_REGISTER_HPP
#define XSIMD_WASM_REGISTER_HPP

#include "xsimd_generic_arch.hpp"
#include "xsimd_register.hpp"

#if XSIMD_WITH_WASM
#include <wasm_simd128.h>
#endif

namespace xsimd
{
    /**
     * @ingroup architectures
     *
     * WASM instructions
     */
    struct wasm : generic
    {
        static constexpr bool supported() noexcept { return XSIMD_WITH_WASM; }
        static constexpr bool available() noexcept { return true; }
        static constexpr bool requires_alignment() noexcept { return true; }
        static constexpr std::size_t alignment() noexcept { return 16; }
        static constexpr char const* name() noexcept { return "wasm"; }
    };

#if XSIMD_WITH_WASM
    namespace types
    {
        XSIMD_DECLARE_SIMD_REGISTER(signed char, wasm, v128_t);
        XSIMD_DECLARE_SIMD_REGISTER(unsigned char, wasm, v128_t);
        XSIMD_DECLARE_SIMD_REGISTER(char, wasm, v128_t);
        XSIMD_DECLARE_SIMD_REGISTER(unsigned short, wasm, v128_t);
        XSIMD_DECLARE_SIMD_REGISTER(short, wasm, v128_t);
        XSIMD_DECLARE_SIMD_REGISTER(unsigned int, wasm, v128_t);
        XSIMD_DECLARE_SIMD_REGISTER(int, wasm, v128_t);
        XSIMD_DECLARE_SIMD_REGISTER(unsigned long int, wasm, v128_t);
        XSIMD_DECLARE_SIMD_REGISTER(long int, wasm, v128_t);
        XSIMD_DECLARE_SIMD_REGISTER(unsigned long long int, wasm, v128_t);
        XSIMD_DECLARE_SIMD_REGISTER(long long int, wasm, v128_t);
        XSIMD_DECLARE_SIMD_REGISTER(float, wasm, v128_t);
        XSIMD_DECLARE_SIMD_REGISTER(double, wasm, v128_t);
    }
#endif
}

#endif
