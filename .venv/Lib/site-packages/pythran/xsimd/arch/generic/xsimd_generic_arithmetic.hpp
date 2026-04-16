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

#ifndef XSIMD_GENERIC_ARITHMETIC_HPP
#define XSIMD_GENERIC_ARITHMETIC_HPP

#include <complex>
#include <limits>
#include <type_traits>

#include "./xsimd_generic_details.hpp"

namespace xsimd
{

    namespace kernel
    {

        using namespace types;

        // bitwise_lshift
        template <class A, class T, class /*=typename std::enable_if<std::is_integral<T>::value, void>::type*/>
        XSIMD_INLINE batch<T, A> bitwise_lshift(batch<T, A> const& self, batch<T, A> const& other, requires_arch<generic>) noexcept
        {
            return detail::apply([](T x, T y) noexcept
                                 { return x << y; },
                                 self, other);
        }

        // bitwise_rshift
        template <class A, class T, class /*=typename std::enable_if<std::is_integral<T>::value, void>::type*/>
        XSIMD_INLINE batch<T, A> bitwise_rshift(batch<T, A> const& self, batch<T, A> const& other, requires_arch<generic>) noexcept
        {
            return detail::apply([](T x, T y) noexcept
                                 { return x >> y; },
                                 self, other);
        }

        // decr
        template <class A, class T>
        XSIMD_INLINE batch<T, A> decr(batch<T, A> const& self, requires_arch<generic>) noexcept
        {
            return self - T(1);
        }

        // decr_if
        template <class A, class T, class Mask>
        XSIMD_INLINE batch<T, A> decr_if(batch<T, A> const& self, Mask const& mask, requires_arch<generic>) noexcept
        {
            return select(mask, decr(self), self);
        }

        // div
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> div(batch<T, A> const& self, batch<T, A> const& other, requires_arch<generic>) noexcept
        {
            return detail::apply([](T x, T y) noexcept -> T
                                 { return x / y; },
                                 self, other);
        }

        // fma
        template <class A, class T>
        XSIMD_INLINE batch<T, A> fma(batch<T, A> const& x, batch<T, A> const& y, batch<T, A> const& z, requires_arch<generic>) noexcept
        {
            return x * y + z;
        }

        template <class A, class T>
        XSIMD_INLINE batch<std::complex<T>, A> fma(batch<std::complex<T>, A> const& x, batch<std::complex<T>, A> const& y, batch<std::complex<T>, A> const& z, requires_arch<generic>) noexcept
        {
            auto res_r = fms(x.real(), y.real(), fms(x.imag(), y.imag(), z.real()));
            auto res_i = fma(x.real(), y.imag(), fma(x.imag(), y.real(), z.imag()));
            return { res_r, res_i };
        }

        // fms
        template <class A, class T>
        XSIMD_INLINE batch<T, A> fms(batch<T, A> const& x, batch<T, A> const& y, batch<T, A> const& z, requires_arch<generic>) noexcept
        {
            return x * y - z;
        }

        template <class A, class T>
        XSIMD_INLINE batch<std::complex<T>, A> fms(batch<std::complex<T>, A> const& x, batch<std::complex<T>, A> const& y, batch<std::complex<T>, A> const& z, requires_arch<generic>) noexcept
        {
            auto res_r = fms(x.real(), y.real(), fma(x.imag(), y.imag(), z.real()));
            auto res_i = fma(x.real(), y.imag(), fms(x.imag(), y.real(), z.imag()));
            return { res_r, res_i };
        }

        // fnma
        template <class A, class T>
        XSIMD_INLINE batch<T, A> fnma(batch<T, A> const& x, batch<T, A> const& y, batch<T, A> const& z, requires_arch<generic>) noexcept
        {
            return -x * y + z;
        }

        template <class A, class T>
        XSIMD_INLINE batch<std::complex<T>, A> fnma(batch<std::complex<T>, A> const& x, batch<std::complex<T>, A> const& y, batch<std::complex<T>, A> const& z, requires_arch<generic>) noexcept
        {
            auto res_r = -fms(x.real(), y.real(), fma(x.imag(), y.imag(), z.real()));
            auto res_i = -fma(x.real(), y.imag(), fms(x.imag(), y.real(), z.imag()));
            return { res_r, res_i };
        }

        // fnms
        template <class A, class T>
        XSIMD_INLINE batch<T, A> fnms(batch<T, A> const& x, batch<T, A> const& y, batch<T, A> const& z, requires_arch<generic>) noexcept
        {
            return -x * y - z;
        }

        template <class A, class T>
        XSIMD_INLINE batch<std::complex<T>, A> fnms(batch<std::complex<T>, A> const& x, batch<std::complex<T>, A> const& y, batch<std::complex<T>, A> const& z, requires_arch<generic>) noexcept
        {
            auto res_r = -fms(x.real(), y.real(), fms(x.imag(), y.imag(), z.real()));
            auto res_i = -fma(x.real(), y.imag(), fma(x.imag(), y.real(), z.imag()));
            return { res_r, res_i };
        }

        // hadd
        template <class A, class T, class /*=typename std::enable_if<std::is_integral<T>::value, void>::type*/>
        XSIMD_INLINE T hadd(batch<T, A> const& self, requires_arch<generic>) noexcept
        {
            alignas(A::alignment()) T buffer[batch<T, A>::size];
            self.store_aligned(buffer);
            T res = 0;
            for (T val : buffer)
            {
                res += val;
            }
            return res;
        }

        // incr
        template <class A, class T>
        XSIMD_INLINE batch<T, A> incr(batch<T, A> const& self, requires_arch<generic>) noexcept
        {
            return self + T(1);
        }

        // incr_if
        template <class A, class T, class Mask>
        XSIMD_INLINE batch<T, A> incr_if(batch<T, A> const& self, Mask const& mask, requires_arch<generic>) noexcept
        {
            return select(mask, incr(self), self);
        }

        // mul
        template <class A, class T, class /*=typename std::enable_if<std::is_integral<T>::value, void>::type*/>
        XSIMD_INLINE batch<T, A> mul(batch<T, A> const& self, batch<T, A> const& other, requires_arch<generic>) noexcept
        {
            return detail::apply([](T x, T y) noexcept -> T
                                 { return x * y; },
                                 self, other);
        }

        // rotl
        template <class A, class T, class STy>
        XSIMD_INLINE batch<T, A> rotl(batch<T, A> const& self, STy other, requires_arch<generic>) noexcept
        {
            constexpr auto N = std::numeric_limits<T>::digits;
            return (self << other) | (self >> (N - other));
        }

        // rotr
        template <class A, class T, class STy>
        XSIMD_INLINE batch<T, A> rotr(batch<T, A> const& self, STy other, requires_arch<generic>) noexcept
        {
            constexpr auto N = std::numeric_limits<T>::digits;
            return (self >> other) | (self << (N - other));
        }

        // sadd
        template <class A>
        XSIMD_INLINE batch<float, A> sadd(batch<float, A> const& self, batch<float, A> const& other, requires_arch<generic>) noexcept
        {
            return add(self, other); // no saturated arithmetic on floating point numbers
        }
        template <class A, class T, class /*=typename std::enable_if<std::is_integral<T>::value, void>::type*/>
        XSIMD_INLINE batch<T, A> sadd(batch<T, A> const& self, batch<T, A> const& other, requires_arch<generic>) noexcept
        {
            if (std::is_signed<T>::value)
            {
                auto mask = (other >> (8 * sizeof(T) - 1));
                auto self_pos_branch = min(std::numeric_limits<T>::max() - other, self);
                auto self_neg_branch = max(std::numeric_limits<T>::min() - other, self);
                return other + select(batch_bool<T, A>(mask.data), self_neg_branch, self_pos_branch);
            }
            else
            {
                const auto diffmax = std::numeric_limits<T>::max() - self;
                const auto mindiff = min(diffmax, other);
                return self + mindiff;
            }
        }
        template <class A>
        XSIMD_INLINE batch<double, A> sadd(batch<double, A> const& self, batch<double, A> const& other, requires_arch<generic>) noexcept
        {
            return add(self, other); // no saturated arithmetic on floating point numbers
        }

        // ssub
        template <class A>
        XSIMD_INLINE batch<float, A> ssub(batch<float, A> const& self, batch<float, A> const& other, requires_arch<generic>) noexcept
        {
            return sub(self, other); // no saturated arithmetic on floating point numbers
        }
        template <class A, class T, class /*=typename std::enable_if<std::is_integral<T>::value, void>::type*/>
        XSIMD_INLINE batch<T, A> ssub(batch<T, A> const& self, batch<T, A> const& other, requires_arch<generic>) noexcept
        {
            if (std::is_signed<T>::value)
            {
                return sadd(self, -other);
            }
            else
            {
                const auto diff = min(self, other);
                return self - diff;
            }
        }
        template <class A>
        XSIMD_INLINE batch<double, A> ssub(batch<double, A> const& self, batch<double, A> const& other, requires_arch<generic>) noexcept
        {
            return sub(self, other); // no saturated arithmetic on floating point numbers
        }

    }

}

#endif
