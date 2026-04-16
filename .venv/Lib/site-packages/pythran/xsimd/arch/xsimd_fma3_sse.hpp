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

#ifndef XSIMD_FMA3_SSE_HPP
#define XSIMD_FMA3_SSE_HPP

#include "../types/xsimd_fma3_sse_register.hpp"

namespace xsimd
{

    namespace kernel
    {
        using namespace types;
        // fnma
        template <class A>
        XSIMD_INLINE batch<float, A> fnma(batch<float, A> const& x, batch<float, A> const& y, batch<float, A> const& z, requires_arch<fma3<sse4_2>>) noexcept
        {
            return _mm_fnmadd_ps(x, y, z);
        }

        template <class A>
        XSIMD_INLINE batch<double, A> fnma(batch<double, A> const& x, batch<double, A> const& y, batch<double, A> const& z, requires_arch<fma3<sse4_2>>) noexcept
        {
            return _mm_fnmadd_pd(x, y, z);
        }

        // fnms
        template <class A>
        XSIMD_INLINE batch<float, A> fnms(batch<float, A> const& x, batch<float, A> const& y, batch<float, A> const& z, requires_arch<fma3<sse4_2>>) noexcept
        {
            return _mm_fnmsub_ps(x, y, z);
        }

        template <class A>
        XSIMD_INLINE batch<double, A> fnms(batch<double, A> const& x, batch<double, A> const& y, batch<double, A> const& z, requires_arch<fma3<sse4_2>>) noexcept
        {
            return _mm_fnmsub_pd(x, y, z);
        }

        // fma
        template <class A>
        XSIMD_INLINE batch<float, A> fma(batch<float, A> const& x, batch<float, A> const& y, batch<float, A> const& z, requires_arch<fma3<sse4_2>>) noexcept
        {
            return _mm_fmadd_ps(x, y, z);
        }

        template <class A>
        XSIMD_INLINE batch<double, A> fma(batch<double, A> const& x, batch<double, A> const& y, batch<double, A> const& z, requires_arch<fma3<sse4_2>>) noexcept
        {
            return _mm_fmadd_pd(x, y, z);
        }

        // fms
        template <class A>
        XSIMD_INLINE batch<float, A> fms(batch<float, A> const& x, batch<float, A> const& y, batch<float, A> const& z, requires_arch<fma3<sse4_2>>) noexcept
        {
            return _mm_fmsub_ps(x, y, z);
        }

        template <class A>
        XSIMD_INLINE batch<double, A> fms(batch<double, A> const& x, batch<double, A> const& y, batch<double, A> const& z, requires_arch<fma3<sse4_2>>) noexcept
        {
            return _mm_fmsub_pd(x, y, z);
        }

    }

}

#endif
