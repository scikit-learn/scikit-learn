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

#ifndef XSIMD_ISA_HPP
#define XSIMD_ISA_HPP

#include "../config/xsimd_arch.hpp"

#include "./xsimd_generic_fwd.hpp"

#if XSIMD_WITH_EMULATED
#include "./xsimd_emulated.hpp"
#endif

#if XSIMD_WITH_SSE2
#include "./xsimd_sse2.hpp"
#endif

#if XSIMD_WITH_SSE3
#include "./xsimd_sse3.hpp"
#endif

#if XSIMD_WITH_SSSE3
#include "./xsimd_ssse3.hpp"
#endif

#if XSIMD_WITH_SSE4_1
#include "./xsimd_sse4_1.hpp"
#endif

#if XSIMD_WITH_SSE4_2
#include "./xsimd_sse4_2.hpp"
#endif

#if XSIMD_WITH_FMA3_SSE
#include "./xsimd_fma3_sse.hpp"
#endif

#if XSIMD_WITH_FMA4
#include "./xsimd_fma4.hpp"
#endif

#if XSIMD_WITH_AVX
#include "./xsimd_avx.hpp"
#endif

#if XSIMD_WITH_FMA3_AVX
#include "./xsimd_fma3_avx.hpp"
#endif

#if XSIMD_WITH_AVXVNNI
#include "./xsimd_avxvnni.hpp"
#endif

#if XSIMD_WITH_AVX2
#include "./xsimd_avx2.hpp"
#endif

#if XSIMD_WITH_FMA3_AVX2
#include "./xsimd_fma3_avx2.hpp"
#endif

#if XSIMD_WITH_AVX512F
#include "./xsimd_avx512f.hpp"
#endif

#if XSIMD_WITH_AVX512BW
#include "./xsimd_avx512bw.hpp"
#endif

#if XSIMD_WITH_AVX512ER
#include "./xsimd_avx512er.hpp"
#endif

#if XSIMD_WITH_AVX512PF
#include "./xsimd_avx512pf.hpp"
#endif

#if XSIMD_WITH_AVX512IFMA
#include "./xsimd_avx512ifma.hpp"
#endif

#if XSIMD_WITH_AVX512VBMI
#include "./xsimd_avx512vbmi.hpp"
#endif

#if XSIMD_WITH_AVX512VNNI_AVX512BW
#include "./xsimd_avx512vnni_avx512bw.hpp"
#endif

#if XSIMD_WITH_AVX512VNNI_AVX512VBMI
#include "./xsimd_avx512vnni_avx512vbmi.hpp"
#endif

#if XSIMD_WITH_NEON
#include "./xsimd_neon.hpp"
#endif

#if XSIMD_WITH_NEON64
#include "./xsimd_neon64.hpp"
#endif

#if XSIMD_WITH_I8MM_NEON64
#include "./xsimd_i8mm_neon64.hpp"
#endif

#if XSIMD_WITH_SVE
#include "./xsimd_sve.hpp"
#endif

#if XSIMD_WITH_RVV
#include "./xsimd_rvv.hpp"
#endif

#if XSIMD_WITH_WASM
#include "./xsimd_wasm.hpp"
#endif

// Must come last to have access to all conversion specializations.
#include "./xsimd_generic.hpp"

#endif
