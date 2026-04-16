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

#ifndef XSIMD_EMULATED_HPP
#define XSIMD_EMULATED_HPP

#include <complex>
#include <limits>
#include <numeric>
#include <type_traits>

#include "../arch/xsimd_scalar.hpp"

#include "../types/xsimd_emulated_register.hpp"
#include "../types/xsimd_utils.hpp"

namespace xsimd
{
    template <typename T, class A, bool... Values>
    struct batch_bool_constant;

    template <class T_out, class T_in, class A>
    XSIMD_INLINE batch<T_out, A> bitwise_cast(batch<T_in, A> const& x) noexcept;

    template <typename T, class A, T... Values>
    struct batch_constant;

    namespace kernel
    {
        using namespace types;

        // fwd
        template <class A, class T, size_t I>
        XSIMD_INLINE batch<T, A> insert(batch<T, A> const& self, T val, index<I>, requires_arch<generic>) noexcept;
        template <class A, typename T, typename ITy, ITy... Indices>
        XSIMD_INLINE batch<T, A> shuffle(batch<T, A> const& x, batch<T, A> const& y, batch_constant<ITy, A, Indices...>, requires_arch<generic>) noexcept;

        namespace detail
        {
            template <size_t I, class F, class... Bs>
            auto emulated_apply(F func, Bs const&... bs) -> decltype(func(bs.data[I]...))
            {
                return func(bs.data[I]...);
            }

            template <class F, class B, class... Bs, size_t... Is>
            auto emulated_apply(F func, ::xsimd::detail::index_sequence<Is...>, B const& b, Bs const&... bs) -> std::array<decltype(func(b.data[0], bs.data[0]...)), B::size>
            {
                return { emulated_apply<Is>(func, b, bs...)... };
            }

            template <class B, class F, class... Bs>
            auto emulated_apply(F func, B const& b, Bs const&... bs) -> std::array<decltype(func(b.data[0], bs.data[0]...)), B::size>
            {
                return emulated_apply(func, ::xsimd::detail::make_index_sequence<B::size>(), b, bs...);
            }
        }

        // abs
        template <class A, class T, size_t N = 8 * sizeof(T) * batch<T, A>::size>
        XSIMD_INLINE batch<T, A> abs(batch<T, A> const& self, requires_arch<emulated<N>>) noexcept
        {
            return detail::emulated_apply([](T v)
                                          { return xsimd::abs(v); },
                                          self);
        }

        // add
        template <class A, class T, size_t N = 8 * sizeof(T) * batch<T, A>::size>
        XSIMD_INLINE batch<T, A> add(batch<T, A> const& self, batch<T, A> const& other, requires_arch<emulated<N>>) noexcept
        {
            return detail::emulated_apply([](T v0, T v1)
                                          { return xsimd::add(v0, v1); },
                                          self, other);
        }

        // all
        template <class A, class T, size_t N = 8 * sizeof(T) * batch<T, A>::size>
        XSIMD_INLINE bool all(batch_bool<T, A> const& self, requires_arch<emulated<N>>) noexcept
        {
            return std::all_of(self.data.begin(), self.data.end(), [](T v)
                               { return bool(v); });
        }

        // any
        template <class A, class T, size_t N = 8 * sizeof(T) * batch<T, A>::size>
        XSIMD_INLINE bool any(batch_bool<T, A> const& self, requires_arch<emulated<N>>) noexcept
        {
            return std::any_of(self.data.begin(), self.data.end(), [](T v)
                               { return bool(v); });
        }

        // batch_bool_cast
        template <class A, class T_out, class T_in, size_t N = 8 * sizeof(T_in) * batch<T_in, A>::size>
        XSIMD_INLINE batch_bool<T_out, A> batch_bool_cast(batch_bool<T_in, A> const& self, batch_bool<T_out, A> const&, requires_arch<emulated<N>>) noexcept
        {
            return { self.data };
        }

        // bitwise_and
        template <class A, class T, size_t N = 8 * sizeof(T) * batch<T, A>::size>
        XSIMD_INLINE batch<T, A> bitwise_and(batch<T, A> const& self, batch<T, A> const& other, requires_arch<emulated<N>>) noexcept
        {
            return detail::emulated_apply([](T v0, T v1)
                                          { return xsimd::bitwise_and(v0, v1); },
                                          self, other);
        }

        template <class A, class T, size_t N = 8 * sizeof(T) * batch<T, A>::size>
        XSIMD_INLINE batch_bool<T, A> bitwise_and(batch_bool<T, A> const& self, batch_bool<T, A> const& other, requires_arch<emulated<N>>) noexcept
        {
            return detail::emulated_apply([](bool v0, bool v1)
                                          { return xsimd::bitwise_and(v0, v1); },
                                          self, other);
        }

        // bitwise_andnot
        template <class A, class T, size_t N = 8 * sizeof(T) * batch<T, A>::size>
        XSIMD_INLINE batch<T, A> bitwise_andnot(batch<T, A> const& self, batch<T, A> const& other, requires_arch<emulated<N>>) noexcept
        {
            return detail::emulated_apply([](T v0, T v1)
                                          { return xsimd::bitwise_andnot(v0, v1); },
                                          self, other);
        }

        template <class A, class T, size_t N = 8 * sizeof(T) * batch<T, A>::size>
        XSIMD_INLINE batch_bool<T, A> bitwise_andnot(batch_bool<T, A> const& self, batch_bool<T, A> const& other, requires_arch<emulated<N>>) noexcept
        {
            return detail::emulated_apply([](bool v0, bool v1)
                                          { return xsimd::bitwise_andnot(v0, v1); },
                                          self, other);
        }

        // bitwise_lshift
        template <class A, class T, size_t N = 8 * sizeof(T) * batch<T, A>::size>
        XSIMD_INLINE batch<T, A> bitwise_lshift(batch<T, A> const& self, int32_t other, requires_arch<emulated<N>>) noexcept
        {
            return detail::emulated_apply([other](T v)
                                          { return xsimd::bitwise_lshift(v, other); },
                                          self);
        }

        // bitwise_not
        template <class A, class T, size_t N = 8 * sizeof(T) * batch<T, A>::size>
        XSIMD_INLINE batch<T, A> bitwise_not(batch<T, A> const& self, requires_arch<emulated<N>>) noexcept
        {
            return detail::emulated_apply([](T v)
                                          { return xsimd::bitwise_not(v); },
                                          self);
        }

        template <class A, class T, size_t N = 8 * sizeof(T) * batch<T, A>::size>
        XSIMD_INLINE batch_bool<T, A> bitwise_not(batch_bool<T, A> const& self, requires_arch<emulated<N>>) noexcept
        {
            return detail::emulated_apply([](bool v)
                                          { return xsimd::bitwise_not(v); },
                                          self);
        }

        // bitwise_or
        template <class A, class T, size_t N = 8 * sizeof(T) * batch<T, A>::size>
        XSIMD_INLINE batch<T, A> bitwise_or(batch<T, A> const& self, batch<T, A> const& other, requires_arch<emulated<N>>) noexcept
        {
            return detail::emulated_apply([](T v0, T v1)
                                          { return xsimd::bitwise_or(v0, v1); },
                                          self, other);
        }

        template <class A, class T, size_t N = 8 * sizeof(T) * batch<T, A>::size>
        XSIMD_INLINE batch_bool<T, A> bitwise_or(batch_bool<T, A> const& self, batch_bool<T, A> const& other, requires_arch<emulated<N>>) noexcept
        {
            return detail::emulated_apply([](bool v0, bool v1)
                                          { return xsimd::bitwise_or(v0, v1); },
                                          self, other);
        }

        // bitwise_rshift
        template <class A, class T, size_t N = 8 * sizeof(T) * batch<T, A>::size>
        XSIMD_INLINE batch<T, A> bitwise_rshift(batch<T, A> const& self, int32_t other, requires_arch<emulated<N>>) noexcept
        {
            return detail::emulated_apply([other](T v)
                                          { return xsimd::bitwise_rshift(v, other); },
                                          self);
        }

        // bitwise_xor
        template <class A, class T, size_t N = 8 * sizeof(T) * batch<T, A>::size>
        XSIMD_INLINE batch<T, A> bitwise_xor(batch<T, A> const& self, batch<T, A> const& other, requires_arch<emulated<N>>) noexcept
        {
            return detail::emulated_apply([](T v0, T v1)
                                          { return xsimd::bitwise_xor(v0, v1); },
                                          self, other);
        }

        template <class A, class T, size_t N = 8 * sizeof(T) * batch<T, A>::size>
        XSIMD_INLINE batch_bool<T, A> bitwise_xor(batch_bool<T, A> const& self, batch_bool<T, A> const& other, requires_arch<emulated<N>>) noexcept
        {
            return detail::emulated_apply([](bool v0, bool v1)
                                          { return xsimd::bitwise_xor(v0, v1); },
                                          self, other);
        }

        // bitwise_cast
        template <class A, class T_in, class T_out, size_t N = 8 * sizeof(T_in) * batch<T_in, A>::size>
        XSIMD_INLINE batch<T_out, A> bitwise_cast(batch<T_in, A> const& self, batch<T_out, A> const&, requires_arch<emulated<N>>) noexcept
        {
            constexpr size_t size = batch<T_out, A>::size;
            std::array<T_out, size> result;
            char* raw_data = reinterpret_cast<char*>(result.data());
            const char* raw_input = reinterpret_cast<const char*>(self.data.data());
            memcpy(raw_data, raw_input, size * sizeof(T_out));
            return result;
        }

        // broadcast
        template <class A, class T, size_t N = 8 * sizeof(T) * batch<T, A>::size>
        batch<T, A> XSIMD_INLINE broadcast(T val, requires_arch<emulated<N>>) noexcept
        {
            constexpr size_t size = batch<T, A>::size;
            std::array<T, size> r;
            std::fill(r.begin(), r.end(), val);
            return r;
        }

        // store_complex
        namespace detail
        {
            // complex_low
            template <class A, typename T, size_t N = 8 * sizeof(T) * batch<T, A>::size>
            XSIMD_INLINE batch<T, A> complex_low(batch<std::complex<T>, A> const& self, requires_arch<emulated<N>>) noexcept
            {
                constexpr size_t size = batch<T, A>::size;
                std::array<T, size> result;
                for (size_t i = 0; i < size / 2; ++i)
                {
                    result[2 * i] = self.real().data[i];
                    result[1 + 2 * i] = self.imag().data[i];
                }
                return result;
            }
            // complex_high
            template <class A, typename T, size_t N = 8 * sizeof(T) * batch<T, A>::size>
            XSIMD_INLINE batch<T, A> complex_high(batch<std::complex<T>, A> const& self, requires_arch<emulated<N>>) noexcept
            {
                constexpr size_t size = batch<T, A>::size;
                std::array<T, size> result;
                for (size_t i = 0; i < size / 2; ++i)
                {
                    result[2 * i] = self.real().data[i + size / 2];
                    result[1 + 2 * i] = self.imag().data[i + size / 2];
                }
                return result;
            }
        }

        // decr_if
        template <class A, class T, size_t N = 8 * sizeof(T) * batch<T, A>::size>
        XSIMD_INLINE batch<T, A> decr_if(batch<T, A> const& self, batch_bool<T, A> const& mask, requires_arch<emulated<N>>) noexcept
        {
            return self - batch<T, A>(mask.data);
        }

        // div
        template <class A, class T, size_t N = 8 * sizeof(T) * batch<T, A>::size>
        XSIMD_INLINE batch<T, A> div(batch<T, A> const& self, batch<T, A> const& other, requires_arch<emulated<N>>) noexcept
        {
            return detail::emulated_apply([](T v0, T v1)
                                          { return xsimd::div(v0, v1); },
                                          self, other);
        }

        // fast_cast
        namespace detail
        {
            template <class A, size_t N = 8 * sizeof(float) * batch<float, A>::size>
            XSIMD_INLINE batch<float, A> fast_cast(batch<int32_t, A> const& self, batch<float, A> const&, requires_arch<emulated<N>>) noexcept
            {
                return detail::emulated_apply([](int32_t v)
                                              { return float(v); },
                                              self);
            }

            template <class A, size_t N = 8 * sizeof(float) * batch<float, A>::size>
            XSIMD_INLINE batch<float, A> fast_cast(batch<uint32_t, A> const& self, batch<float, A> const&, requires_arch<emulated<N>>) noexcept
            {
                return detail::emulated_apply([](uint32_t v)
                                              { return float(v); },
                                              self);
            }

            template <class A, size_t N = 8 * sizeof(double) * batch<double, A>::size>
            XSIMD_INLINE batch<double, A> fast_cast(batch<int64_t, A> const& self, batch<double, A> const&, requires_arch<emulated<N>>) noexcept
            {
                return detail::emulated_apply([](int64_t v)
                                              { return double(v); },
                                              self);
            }

            template <class A, size_t N = 8 * sizeof(double) * batch<double, A>::size>
            XSIMD_INLINE batch<double, A> fast_cast(batch<uint64_t, A> const& self, batch<double, A> const&, requires_arch<emulated<N>>) noexcept
            {
                return detail::emulated_apply([](uint64_t v)
                                              { return double(v); },
                                              self);
            }

            template <class A, size_t N = 8 * sizeof(int32_t) * batch<int32_t, A>::size>
            XSIMD_INLINE batch<int32_t, A> fast_cast(batch<float, A> const& self, batch<int32_t, A> const&, requires_arch<emulated<N>>) noexcept
            {
                return detail::emulated_apply([](float v)
                                              { return int32_t(v); },
                                              self);
            }

            template <class A, size_t N = 8 * sizeof(double) * batch<double, A>::size>
            XSIMD_INLINE batch<int64_t, A> fast_cast(batch<double, A> const& self, batch<int64_t, A> const&, requires_arch<emulated<N>>) noexcept
            {
                return detail::emulated_apply([](double v)
                                              { return int64_t(v); },
                                              self);
            }
        }

        // eq
        template <class A, class T, size_t N = 8 * sizeof(T) * batch<T, A>::size>
        XSIMD_INLINE batch_bool<T, emulated<N>> eq(batch<T, emulated<N>> const& self, batch<T, emulated<N>> const& other, requires_arch<emulated<N>>) noexcept
        {
            return detail::emulated_apply([](T v0, T v1)
                                          { return xsimd::eq(v0, v1); },
                                          self, other);
        }

        template <class A, class T, size_t N = 8 * sizeof(T) * batch_bool<T, A>::size>
        XSIMD_INLINE batch_bool<T, emulated<N>> eq(batch_bool<T, emulated<N>> const& self, batch_bool<T, emulated<N>> const& other, requires_arch<emulated<N>>) noexcept
        {
            return detail::emulated_apply([](bool v0, bool v1)
                                          { return xsimd::eq(v0, v1); },
                                          self, other);
        }

        // from_bool
        template <class A, class T, size_t N = 8 * sizeof(T) * batch<T, A>::size>
        XSIMD_INLINE batch<T, A> from_bool(batch_bool<T, A> const& self, requires_arch<emulated<N>>) noexcept
        {
            return detail::emulated_apply([](bool v)
                                          { return T(v); },
                                          self);
        }

        // from_mask
        template <class A, class T, size_t N = 8 * sizeof(T) * batch<T, A>::size>
        XSIMD_INLINE batch_bool<T, A> from_mask(batch_bool<T, A> const&, uint64_t mask, requires_arch<emulated<N>>) noexcept
        {
            constexpr size_t size = batch<T, A>::size;
            std::array<bool, size> vmask;
            for (size_t i = 0; i < size; ++i)
                vmask[i] = (mask >> i) & 1u;
            return vmask;
        }

        // ge
        template <class A, class T, size_t N = 8 * sizeof(T) * batch<T, A>::size>
        XSIMD_INLINE batch_bool<T, emulated<N>> ge(batch<T, emulated<N>> const& self, batch<T, emulated<N>> const& other, requires_arch<emulated<N>>) noexcept
        {
            return detail::emulated_apply([](T v0, T v1)
                                          { return xsimd::ge(v0, v1); },
                                          self, other);
        }

        // gt
        template <class A, class T, size_t N = 8 * sizeof(T) * batch<T, A>::size>
        XSIMD_INLINE batch_bool<T, emulated<N>> gt(batch<T, emulated<N>> const& self, batch<T, emulated<N>> const& other, requires_arch<emulated<N>>) noexcept
        {
            return detail::emulated_apply([](T v0, T v1)
                                          { return xsimd::gt(v0, v1); },
                                          self, other);
        }

        // haddp
        template <class A, typename T, size_t N = 8 * sizeof(T) * batch<T, A>::size>
        XSIMD_INLINE batch<T, A> haddp(batch<T, A> const* row, requires_arch<emulated<N>>) noexcept
        {
            constexpr size_t size = batch<T, A>::size;
            std::array<T, size> r;
            for (size_t i = 0; i < size; ++i)
                r[i] = std::accumulate(row[i].data.begin() + 1, row[i].data.end(), row[i].data.front());
            return r;
        }

        // incr_if
        template <class A, class T, size_t N = 8 * sizeof(T) * batch<T, A>::size>
        XSIMD_INLINE batch<T, A> incr_if(batch<T, A> const& self, batch_bool<T, A> const& mask, requires_arch<emulated<N>>) noexcept
        {
            return self + batch<T, A>(mask.data);
        }

        // insert
        template <class A, class T, size_t I, size_t N = 8 * sizeof(T) * batch<T, A>::size>
        XSIMD_INLINE batch<T, A> insert(batch<T, A> const& self, T val, index<I>, requires_arch<emulated<N>>) noexcept
        {
            batch<T, A> other = self;
            other.data[I] = val;
            return other;
        }

        // isnan
        template <class A, typename T, size_t N = 8 * sizeof(T) * batch<T, A>::size, class = typename std::enable_if<std::is_floating_point<T>::value, void>::type>
        XSIMD_INLINE batch_bool<T, A> isnan(batch<T, A> const& self, requires_arch<emulated<N>>) noexcept
        {
            return detail::emulated_apply([](T v)
                                          { return xsimd::isnan(v); },
                                          self);
        }

        // load_aligned
        template <class A, typename T, size_t N = 8 * sizeof(T) * batch<T, A>::size>
        XSIMD_INLINE batch<T, A> load_aligned(T const* mem, convert<T>, requires_arch<emulated<N>>) noexcept
        {
            constexpr size_t size = batch<T, A>::size;
            std::array<T, size> res;
            std::copy(mem, mem + size, res.begin());
            return res;
        }

        // load_unaligned
        template <class A, typename T, size_t N = 8 * sizeof(T) * batch<T, A>::size>
        XSIMD_INLINE batch<T, A> load_unaligned(T const* mem, convert<T>, requires_arch<emulated<N>>) noexcept
        {
            constexpr size_t size = batch<T, A>::size;
            std::array<T, size> res;
            std::copy(mem, mem + size, res.begin());
            return res;
        }

        // load_complex
        namespace detail
        {
            template <class A, typename T, size_t N = 8 * sizeof(T) * batch<T, A>::size>
            XSIMD_INLINE batch<std::complex<T>, A> load_complex(batch<T, A> const& hi, batch<T, A> const& lo, requires_arch<emulated<N>>) noexcept
            {
                constexpr size_t size = batch<T, A>::size;
                std::array<T, size> real, imag;
                for (size_t i = 0; i < size / 2; ++i)
                {
                    real[i] = hi.data[2 * i];
                    imag[i] = hi.data[1 + 2 * i];
                }
                for (size_t i = 0; i < size / 2; ++i)
                {
                    real[size / 2 + i] = lo.data[2 * i];
                    imag[size / 2 + i] = lo.data[1 + 2 * i];
                }
                return { real, imag };
            }
        }

        // le
        template <class A, class T, size_t N = 8 * sizeof(T) * batch<T, A>::size>
        XSIMD_INLINE batch_bool<T, emulated<N>> le(batch<T, emulated<N>> const& self, batch<T, emulated<N>> const& other, requires_arch<emulated<N>>) noexcept
        {
            return detail::emulated_apply([](T v0, T v1)
                                          { return xsimd::le(v0, v1); },
                                          self, other);
        }

        // lt
        template <class A, class T, size_t N = 8 * sizeof(T) * batch<T, A>::size>
        XSIMD_INLINE batch_bool<T, emulated<N>> lt(batch<T, emulated<N>> const& self, batch<T, emulated<N>> const& other, requires_arch<emulated<N>>) noexcept
        {
            return detail::emulated_apply([](T v0, T v1)
                                          { return xsimd::lt(v0, v1); },
                                          self, other);
        }

        // mask
        template <class A, class T, size_t N = 8 * sizeof(T) * batch<T, A>::size>
        XSIMD_INLINE uint64_t mask(batch_bool<T, A> const& self, requires_arch<emulated<N>>) noexcept
        {
            constexpr size_t size = batch<T, A>::size;
            uint64_t res = 0;
            for (size_t i = 0; i < size; ++i)
                res |= (self.data[i] ? 1u : 0u) << i;
            return res;
        }

        // max
        template <class A, typename T, size_t N = 8 * sizeof(T) * batch<T, A>::size>
        XSIMD_INLINE batch<T, A> max(batch<T, A> const& self, batch<T, A> const& other, requires_arch<emulated<N>>) noexcept
        {
            return detail::emulated_apply([](T v0, T v1)
                                          { return xsimd::max(v0, v1); },
                                          self, other);
        }

        // min
        template <class A, typename T, size_t N = 8 * sizeof(T) * batch<T, A>::size>
        XSIMD_INLINE batch<T, A> min(batch<T, A> const& self, batch<T, A> const& other, requires_arch<emulated<N>>) noexcept
        {
            return detail::emulated_apply([](T v0, T v1)
                                          { return xsimd::min(v0, v1); },
                                          self, other);
        }

        // mul
        template <class A, typename T, size_t N = 8 * sizeof(T) * batch<T, A>::size>
        XSIMD_INLINE batch<T, A> mul(batch<T, A> const& self, batch<T, A> const& other, requires_arch<emulated<N>>) noexcept
        {
            return detail::emulated_apply([](T v0, T v1)
                                          { return xsimd::mul(v0, v1); },
                                          self, other);
        }

        // nearbyint_as_int
        template <class A, typename T, size_t N = 8 * sizeof(T) * batch<T, A>::size>
        XSIMD_INLINE batch<as_integer_t<T>, A> nearbyint_as_int(batch<T, A> const& self,
                                                                requires_arch<emulated<N>>) noexcept
        {
            return detail::emulated_apply([](T v)
                                          { return xsimd::nearbyint_as_int(v); },
                                          self);
        }

        // neg
        template <class A, class T, size_t N = 8 * sizeof(T) * batch<T, A>::size>
        XSIMD_INLINE batch<T, A> neg(batch<T, A> const& self, requires_arch<emulated<N>>) noexcept
        {
            return detail::emulated_apply([](T v)
                                          { return xsimd::neg(v); },
                                          self);
        }

        // neq
        template <class A, class T, size_t N = 8 * sizeof(T) * batch<T, A>::size>
        XSIMD_INLINE batch_bool<T, A> neq(batch<T, A> const& self, batch<T, A> const& other, requires_arch<emulated<N>>) noexcept
        {
            return detail::emulated_apply([](T v0, T v1)
                                          { return xsimd::neq(v0, v1); },
                                          self, other);
        }

        template <class A, class T, size_t N = 8 * sizeof(T) * batch<T, A>::size>
        XSIMD_INLINE batch_bool<T, A> neq(batch_bool<T, A> const& self, batch_bool<T, A> const& other, requires_arch<emulated<N>>) noexcept
        {
            return detail::emulated_apply([](bool v0, bool v1)
                                          { return xsimd::neq(v0, v1); },
                                          self, other);
        }

        // reduce_add
        template <class A, class T, size_t N = 8 * sizeof(T) * batch<T, A>::size>
        XSIMD_INLINE T reduce_add(batch<T, A> const& self, requires_arch<emulated<N>>) noexcept
        {
            constexpr size_t size = batch<T, A>::size;
            std::array<T, size> buffer;
            self.store_unaligned(buffer.data());
            return std::accumulate(buffer.begin() + 1, buffer.end(), *buffer.begin());
        }

        // reduce_max
        template <class A, class T, size_t N = 8 * sizeof(T) * batch<T, A>::size>
        XSIMD_INLINE T reduce_max(batch<T, A> const& self, requires_arch<emulated<N>>) noexcept
        {
            return std::accumulate(self.data.begin() + 1, self.data.end(), *self.data.begin(), [](T const& x, T const& y)
                                   { return xsimd::max(x, y); });
        }

        // reduce_min
        template <class A, class T, size_t N = 8 * sizeof(T) * batch<T, A>::size>
        XSIMD_INLINE T reduce_min(batch<T, A> const& self, requires_arch<emulated<N>>) noexcept
        {
            return std::accumulate(self.data.begin() + 1, self.data.end(), *self.data.begin(), [](T const& x, T const& y)
                                   { return xsimd::min(x, y); });
        }

        // rsqrt
        template <class A, class T, size_t N = 8 * sizeof(T) * batch<T, A>::size>
        XSIMD_INLINE batch<T, A> rsqrt(batch<T, A> const& self, requires_arch<emulated<N>>) noexcept
        {
            return detail::emulated_apply([](T v)
                                          { return xsimd::rsqrt(v); },
                                          self);
        }

        // select
        template <class A, class T, size_t N = 8 * sizeof(T) * batch<T, A>::size>
        XSIMD_INLINE batch<T, A> select(batch_bool<T, A> const& cond, batch<T, A> const& true_br, batch<T, A> const& false_br, requires_arch<emulated<N>>) noexcept
        {
            return detail::emulated_apply([](bool c, T t, T f)
                                          { return xsimd::select(c, t, f); },
                                          cond, true_br, false_br);
        }

        template <class A, class T, bool... Values>
        XSIMD_INLINE batch<T, A> select(batch_bool_constant<T, A, Values...> const& cond, batch<T, A> const& true_br, batch<T, A> const& false_br, requires_arch<emulated<8 * sizeof(T) * batch<T, A>::size>>) noexcept
        {
            constexpr size_t size = batch<T, A>::size;
            static_assert(sizeof...(Values) == size, "consistent init");
            return select((batch_bool<T, A>)cond, true_br, false_br, emulated<8 * sizeof(T) * size> {});
        }

        // shuffle
        template <class A, typename T, class ITy, ITy... Is>
        XSIMD_INLINE batch<T, A> shuffle(batch<T, A> const& x, batch<float, A> const& y, batch_constant<ITy, A, Is...> mask, requires_arch<emulated<batch<T, A>::size>>) noexcept
        {
            constexpr size_t size = batch<T, A>::size;
            batch<ITy, A> bmask = mask;
            std::array<T, size> res;
            for (size_t i = 0; i < size; ++i)
                res[i] = bmask.data[i] < size ? x.data[bmask.data[i]] : y.data[bmask.data[i] - size];
            return res;
        }

        // sqrt
        template <class A, class T, size_t N = 8 * sizeof(T) * batch<T, A>::size>
        XSIMD_INLINE batch<T, A> sqrt(batch<T, A> const& self, requires_arch<emulated<N>>) noexcept
        {
            return detail::emulated_apply([](T v)
                                          { return xsimd::sqrt(v); },
                                          self);
        }

        // slide_left
        template <size_t M, class A, class T, size_t N = 8 * sizeof(T) * batch<T, A>::size>
        XSIMD_INLINE batch<T, A> slide_left(batch<T, A> const& x, requires_arch<emulated<N>>) noexcept
        {
            constexpr size_t size = batch<T, A>::size;
            std::array<T, size> result;
            char* raw_data = reinterpret_cast<char*>(result.data());
            memset(raw_data, 0, M);
            memcpy(raw_data + M, reinterpret_cast<const char*>(x.data.data()), sizeof(T) * result.size() - M);
            return result;
        }

        // slide_right
        template <size_t M, class A, class T, size_t N = 8 * sizeof(T) * batch<T, A>::size>
        XSIMD_INLINE batch<T, A> slide_right(batch<T, A> const& x, requires_arch<emulated<N>>) noexcept
        {
            constexpr size_t size = batch<T, A>::size;
            std::array<T, size> result;
            char* raw_data = reinterpret_cast<char*>(result.data());
            memcpy(raw_data, reinterpret_cast<const char*>(x.data.data()) + M, sizeof(T) * result.size() - M);
            memset(raw_data + sizeof(T) * result.size() - M, 0, M);
            return result;
        }

        // sadd
        template <class A, class T, size_t N = 8 * sizeof(T) * batch<T, A>::size>
        XSIMD_INLINE batch<T, A> sadd(batch<T, A> const& self, batch<T, A> const& other, requires_arch<emulated<N>>) noexcept
        {
            return detail::emulated_apply([](T v0, T v1)
                                          { return xsimd::sadd(v0, v1); },
                                          self, other);
        }

        // set
        template <class A, class T, size_t N, class... Values>
        XSIMD_INLINE batch<T, emulated<N>> set(batch<T, emulated<N>> const&, requires_arch<emulated<N>>, Values... values) noexcept
        {
            static_assert(sizeof...(Values) == batch<T, emulated<N>>::size, "consistent init");
            return { typename batch<T, emulated<N>>::register_type { static_cast<T>(values)... } };
        }

        template <class A, class T, size_t N, class... Values>
        XSIMD_INLINE batch_bool<T, emulated<N>> set(batch_bool<T, emulated<N>> const&, requires_arch<emulated<N>>, Values... values) noexcept
        {
            static_assert(sizeof...(Values) == batch<T, emulated<N>>::size, "consistent init");
            return { std::array<bool, sizeof...(Values)> { static_cast<bool>(values)... } };
        }

        // ssub
        template <class A, class T, size_t N = 8 * sizeof(T) * batch<T, A>::size>
        XSIMD_INLINE batch<T, A> ssub(batch<T, A> const& self, batch<T, A> const& other, requires_arch<emulated<N>>) noexcept
        {
            return detail::emulated_apply([](T v0, T v1)
                                          { return xsimd::ssub(v0, v1); },
                                          self, other);
        }

        // store_aligned
        template <class A, class T, size_t N>
        XSIMD_INLINE void store_aligned(T* mem, batch<T, emulated<N>> const& self, requires_arch<emulated<N>>) noexcept
        {
            std::copy(self.data.begin(), self.data.end(), mem);
        }

        // store_unaligned
        template <class A, class T, size_t N>
        XSIMD_INLINE void store_unaligned(T* mem, batch<T, emulated<N>> const& self, requires_arch<emulated<N>>) noexcept
        {
            std::copy(self.data.begin(), self.data.end(), mem);
        }

        // sub
        template <class A, class T, size_t N = 8 * sizeof(T) * batch<T, A>::size>
        XSIMD_INLINE batch<T, A> sub(batch<T, A> const& self, batch<T, A> const& other, requires_arch<emulated<N>>) noexcept
        {
            return detail::emulated_apply([](T v0, T v1)
                                          { return xsimd::sub(v0, v1); },
                                          self, other);
        }

        // swizzle

        template <class A, typename T, class ITy, ITy... Is>
        XSIMD_INLINE batch<T, A> swizzle(batch<T, A> const& self, batch_constant<ITy, A, Is...> mask, requires_arch<emulated<8 * sizeof(T) * batch<T, A>::size>>) noexcept
        {
            constexpr size_t size = batch<T, A>::size;
            batch<ITy, A> bmask = mask;
            std::array<T, size> res;
            for (size_t i = 0; i < size; ++i)
                res[i] = self.data[bmask.data[i]];
            return res;
        }

        // zip_hi
        template <class A, class T, size_t N = 8 * sizeof(T) * batch<T, A>::size>
        XSIMD_INLINE batch<T, A> zip_hi(batch<T, A> const& self, batch<T, A> const& other, requires_arch<emulated<N>>) noexcept
        {
            constexpr size_t size = batch<T, A>::size;
            // Note: irregular behavior for odd numbers.
            std::array<T, size> res;
            if (size % 2)
            {
                for (size_t i = 0; i < size; ++i)
                    res[i] = (i % 2 ? self : other).data[size / 2 + i / 2];
            }
            else
            {
                for (size_t i = 0; i < size; ++i)
                    res[i] = (i % 2 ? other : self).data[size / 2 + i / 2];
            }
            return res;
        }

        // zip_lo
        template <class A, class T, size_t N = 8 * sizeof(T) * batch<T, A>::size>
        XSIMD_INLINE batch<T, A> zip_lo(batch<T, A> const& self, batch<T, A> const& other, requires_arch<emulated<N>>) noexcept
        {
            constexpr size_t size = batch<T, A>::size;
            // Note: irregular behavior for odd numbers.
            std::array<T, size> res;
            for (size_t i = 0; i < size; ++i)
                res[i] = (i % 2 ? other : self).data[i / 2];
            return res;
        }
    }
}

#endif
