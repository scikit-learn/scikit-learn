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

#ifndef XSIMD_SSE3_HPP
#define XSIMD_SSE3_HPP

#include "../types/xsimd_sse3_register.hpp"
#include <type_traits>

namespace xsimd
{

    namespace kernel
    {
        using namespace types;

        // haddp
        template <class A>
        XSIMD_INLINE batch<float, A> haddp(batch<float, A> const* row, requires_arch<sse3>) noexcept
        {
            return _mm_hadd_ps(_mm_hadd_ps(row[0], row[1]),
                               _mm_hadd_ps(row[2], row[3]));
        }
        template <class A>
        XSIMD_INLINE batch<double, A> haddp(batch<double, A> const* row, requires_arch<sse3>) noexcept
        {
            return _mm_hadd_pd(row[0], row[1]);
        }

        // load_unaligned
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> load_unaligned(T const* mem, convert<T>, requires_arch<sse3>) noexcept
        {
            return _mm_lddqu_si128((__m128i const*)mem);
        }

        // reduce_add
        template <class A>
        XSIMD_INLINE float reduce_add(batch<float, A> const& self, requires_arch<sse3>) noexcept
        {
            __m128 tmp0 = _mm_hadd_ps(self, self);
            __m128 tmp1 = _mm_hadd_ps(tmp0, tmp0);
            return _mm_cvtss_f32(tmp1);
        }
        template <class A>
        XSIMD_INLINE double reduce_add(batch<double, A> const& self, requires_arch<sse3>) noexcept
        {
            __m128d tmp0 = _mm_hadd_pd(self, self);
            return _mm_cvtsd_f64(tmp0);
        }

    }

}

#endif
