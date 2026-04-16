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

#include "xsimd_fma3_sse_register.hpp"
#include "xsimd_fma4_register.hpp"
#include "xsimd_sse2_register.hpp"
#include "xsimd_sse3_register.hpp"
#include "xsimd_sse4_1_register.hpp"
#include "xsimd_sse4_2_register.hpp"

#include "xsimd_avx2_register.hpp"
#include "xsimd_avx_register.hpp"
#include "xsimd_avxvnni_register.hpp"
#include "xsimd_fma3_avx2_register.hpp"
#include "xsimd_fma3_avx_register.hpp"

#include "xsimd_avx512vnni_avx512bw_register.hpp"
#include "xsimd_avx512vnni_avx512vbmi_register.hpp"

#include "xsimd_avx512ifma_register.hpp"
#include "xsimd_avx512vbmi_register.hpp"

#include "xsimd_avx512er_register.hpp"
#include "xsimd_avx512pf_register.hpp"

#include "xsimd_avx512bw_register.hpp"
#include "xsimd_avx512cd_register.hpp"
#include "xsimd_avx512dq_register.hpp"
#include "xsimd_avx512f_register.hpp"

#include "xsimd_i8mm_neon64_register.hpp"

#include "xsimd_neon64_register.hpp"
#include "xsimd_neon_register.hpp"

#include "xsimd_sve_register.hpp"

#include "xsimd_rvv_register.hpp"

#include "xsimd_wasm_register.hpp"

#if XSIMD_WITH_EMULATED
#include "xsimd_emulated_register.hpp"
#endif
