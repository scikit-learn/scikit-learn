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

#ifndef XSIMD_FMA3_AVX2_HPP
#define XSIMD_FMA3_AVX2_HPP

#include "../types/xsimd_fma3_avx2_register.hpp"

// Allow inclusion of xsimd_fma3_avx.hpp
#ifdef XSIMD_FMA3_AVX_HPP
#undef XSIMD_FMA3_AVX_HPP
#define XSIMD_FORCE_FMA3_AVX_HPP
#endif

// Disallow inclusion of ./xsimd_fma3_avx_register.hpp
#ifndef XSIMD_FMA3_AVX_REGISTER_HPP
#define XSIMD_FMA3_AVX_REGISTER_HPP
#define XSIMD_FORCE_FMA3_AVX_REGISTER_HPP
#endif

// Include ./xsimd_fma3_avx.hpp but s/avx/avx2
#define avx avx2
#include "./xsimd_fma3_avx.hpp"
#undef avx
#undef XSIMD_FMA3_AVX_HPP

// Carefully restore guards
#ifdef XSIMD_FORCE_FMA3_AVX_HPP
#define XSIMD_FMA3_AVX_HPP
#undef XSIMD_FORCE_FMA3_AVX_HPP
#endif

#ifdef XSIMD_FORCE_FMA3_AVX_REGISTER_HPP
#undef XSIMD_FMA3_AVX_REGISTER_HPP
#undef XSIMD_FORCE_FMA3_AVX_REGISTER_HPP
#endif

#endif
