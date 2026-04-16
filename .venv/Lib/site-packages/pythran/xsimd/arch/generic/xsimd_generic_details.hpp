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

#ifndef XSIMD_GENERIC_DETAILS_HPP
#define XSIMD_GENERIC_DETAILS_HPP

#include <complex>

#include "../../math/xsimd_rem_pio2.hpp"
#include "../../types/xsimd_generic_arch.hpp"
#include "../../types/xsimd_utils.hpp"
#include "../xsimd_constants.hpp"

namespace xsimd
{
    // Forward declaration. Should we put them in a separate file?
    template <class T, class A>
    XSIMD_INLINE batch<T, A> abs(batch<T, A> const& self) noexcept;
    template <class T, class A>
    XSIMD_INLINE batch<T, A> abs(batch<std::complex<T>, A> const& self) noexcept;
    template <class T, class A>
    XSIMD_INLINE bool any(batch_bool<T, A> const& self) noexcept;
    template <class T, class A>
    XSIMD_INLINE batch<T, A> atan2(batch<T, A> const& self, batch<T, A> const& other) noexcept;
    template <class A, class T_out, class T_in>
    XSIMD_INLINE batch<T_out, A> batch_cast(batch<T_in, A> const&, batch<T_out, A> const& out) noexcept;
    template <class T, class A>
    XSIMD_INLINE batch<T, A> bitofsign(batch<T, A> const& self) noexcept;
    template <class T_out, class T_in, class A>
    XSIMD_INLINE batch<T_out, A> bitwise_cast(batch<T_in, A> const& self) noexcept;
    template <class T, class A>
    XSIMD_INLINE batch<T, A> cos(batch<T, A> const& self) noexcept;
    template <class T, class A>
    XSIMD_INLINE batch<T, A> cosh(batch<T, A> const& self) noexcept;
    template <class T, class A>
    XSIMD_INLINE batch<T, A> exp(batch<T, A> const& self) noexcept;
    template <class T, class A>
    XSIMD_INLINE batch<T, A> fma(batch<T, A> const& x, batch<T, A> const& y, batch<T, A> const& z) noexcept;
    template <class T, class A>
    XSIMD_INLINE batch<T, A> fms(batch<T, A> const& x, batch<T, A> const& y, batch<T, A> const& z) noexcept;
    template <class T, class A>
    XSIMD_INLINE batch<T, A> frexp(const batch<T, A>& x, const batch<as_integer_t<T>, A>& e) noexcept;
    template <class T, class A, uint64_t... Coefs>
    XSIMD_INLINE batch<T, A> horner(const batch<T, A>& self) noexcept;
    template <class T, class A>
    XSIMD_INLINE batch<T, A> hypot(const batch<T, A>& self) noexcept;
    template <class T, class A>
    XSIMD_INLINE batch_bool<T, A> is_even(batch<T, A> const& self) noexcept;
    template <class T, class A>
    XSIMD_INLINE batch_bool<T, A> is_flint(batch<T, A> const& self) noexcept;
    template <class T, class A>
    XSIMD_INLINE batch_bool<T, A> is_odd(batch<T, A> const& self) noexcept;
    template <class T, class A>
    XSIMD_INLINE typename batch<T, A>::batch_bool_type isinf(batch<T, A> const& self) noexcept;
    template <class T, class A>
    XSIMD_INLINE typename batch<T, A>::batch_bool_type isfinite(batch<T, A> const& self) noexcept;
    template <class T, class A>
    XSIMD_INLINE typename batch<T, A>::batch_bool_type isnan(batch<T, A> const& self) noexcept;
    template <class T, class A>
    XSIMD_INLINE batch<T, A> ldexp(const batch<T, A>& x, const batch<as_integer_t<T>, A>& e) noexcept;
    template <class T, class A>
    XSIMD_INLINE batch<T, A> log(batch<T, A> const& self) noexcept;
    template <class T, class A>
    XSIMD_INLINE batch<T, A> nearbyint(batch<T, A> const& self) noexcept;
    template <class T, class A>
    XSIMD_INLINE batch<as_integer_t<T>, A> nearbyint_as_int(const batch<T, A>& x) noexcept;
    template <class T, class A>
    XSIMD_INLINE T reduce_add(batch<T, A> const&) noexcept;
    template <class T, class A>
    XSIMD_INLINE batch<T, A> select(batch_bool<T, A> const&, batch<T, A> const&, batch<T, A> const&) noexcept;
    template <class T, class A>
    XSIMD_INLINE batch<std::complex<T>, A> select(batch_bool<T, A> const&, batch<std::complex<T>, A> const&, batch<std::complex<T>, A> const&) noexcept;
    template <class T, class A>
    XSIMD_INLINE batch<T, A> sign(batch<T, A> const& self) noexcept;
    template <class T, class A>
    XSIMD_INLINE batch<T, A> signnz(batch<T, A> const& self) noexcept;
    template <class T, class A>
    XSIMD_INLINE batch<T, A> sin(batch<T, A> const& self) noexcept;
    template <class T, class A>
    XSIMD_INLINE batch<T, A> sinh(batch<T, A> const& self) noexcept;
    template <class T, class A>
    XSIMD_INLINE std::pair<batch<T, A>, batch<T, A>> sincos(batch<T, A> const& self) noexcept;
    template <class T, class A>
    XSIMD_INLINE batch<T, A> sqrt(batch<T, A> const& self) noexcept;
    template <class T, class A>
    XSIMD_INLINE batch<T, A> tan(batch<T, A> const& self) noexcept;
    template <class T, class A>
    XSIMD_INLINE batch<as_float_t<T>, A> to_float(batch<T, A> const& self) noexcept;
    template <class T, class A>
    XSIMD_INLINE batch<as_integer_t<T>, A> to_int(batch<T, A> const& self) noexcept;
    template <class T, class A>
    XSIMD_INLINE batch<T, A> trunc(batch<T, A> const& self) noexcept;

    namespace kernel
    {

        namespace detail
        {
            template <class F, class A, class T, class... Batches>
            XSIMD_INLINE batch<T, A> apply(F&& func, batch<T, A> const& self, batch<T, A> const& other) noexcept
            {
                constexpr std::size_t size = batch<T, A>::size;
                alignas(A::alignment()) T self_buffer[size];
                alignas(A::alignment()) T other_buffer[size];
                self.store_aligned(&self_buffer[0]);
                other.store_aligned(&other_buffer[0]);
                for (std::size_t i = 0; i < size; ++i)
                {
                    self_buffer[i] = func(self_buffer[i], other_buffer[i]);
                }
                return batch<T, A>::load_aligned(self_buffer);
            }

            template <class U, class F, class A, class T>
            XSIMD_INLINE batch<U, A> apply_transform(F&& func, batch<T, A> const& self) noexcept
            {
                static_assert(batch<T, A>::size == batch<U, A>::size,
                              "Source and destination sizes must match");
                constexpr std::size_t src_size = batch<T, A>::size;
                constexpr std::size_t dest_size = batch<U, A>::size;
                alignas(A::alignment()) T self_buffer[src_size];
                alignas(A::alignment()) U other_buffer[dest_size];
                self.store_aligned(&self_buffer[0]);
                for (std::size_t i = 0; i < src_size; ++i)
                {
                    other_buffer[i] = func(self_buffer[i]);
                }
                return batch<U, A>::load_aligned(other_buffer);
            }
        }

        // some generic fast_cast conversion
        namespace detail
        {
            template <class A>
            XSIMD_INLINE batch<uint8_t, A> fast_cast(batch<int8_t, A> const& self, batch<uint8_t, A> const&, requires_arch<generic>) noexcept
            {
                return bitwise_cast<uint8_t>(self);
            }
            template <class A>
            XSIMD_INLINE batch<uint16_t, A> fast_cast(batch<int16_t, A> const& self, batch<uint16_t, A> const&, requires_arch<generic>) noexcept
            {
                return bitwise_cast<uint16_t>(self);
            }
            template <class A>
            XSIMD_INLINE batch<uint32_t, A> fast_cast(batch<int32_t, A> const& self, batch<uint32_t, A> const&, requires_arch<generic>) noexcept
            {
                return bitwise_cast<uint32_t>(self);
            }
            template <class A>
            XSIMD_INLINE batch<uint64_t, A> fast_cast(batch<int64_t, A> const& self, batch<uint64_t, A> const&, requires_arch<generic>) noexcept
            {
                return bitwise_cast<uint64_t>(self);
            }
            template <class A>
            XSIMD_INLINE batch<int8_t, A> fast_cast(batch<uint8_t, A> const& self, batch<int8_t, A> const&, requires_arch<generic>) noexcept
            {
                return bitwise_cast<int8_t>(self);
            }
            template <class A>
            XSIMD_INLINE batch<int16_t, A> fast_cast(batch<uint16_t, A> const& self, batch<int16_t, A> const&, requires_arch<generic>) noexcept
            {
                return bitwise_cast<int16_t>(self);
            }
            template <class A>
            XSIMD_INLINE batch<int32_t, A> fast_cast(batch<uint32_t, A> const& self, batch<int32_t, A> const&, requires_arch<generic>) noexcept
            {
                return bitwise_cast<int32_t>(self);
            }
            template <class A>
            XSIMD_INLINE batch<int64_t, A> fast_cast(batch<uint64_t, A> const& self, batch<int64_t, A> const&, requires_arch<generic>) noexcept
            {
                return bitwise_cast<int64_t>(self);
            }

            // Provide a generic uint32_t -> float cast only if we have a
            // non-generic int32_t -> float fast_cast
            template <class A, class _ = decltype(fast_cast(std::declval<batch<int32_t, A> const&>(), std::declval<batch<float, A> const&>(), A {}))>
            XSIMD_INLINE batch<float, A> fast_cast(batch<uint32_t, A> const& v, batch<float, A> const&, requires_arch<generic>) noexcept
            {
                // see https://stackoverflow.com/questions/34066228/how-to-perform-uint32-float-conversion-with-sse
                batch<uint32_t, A> msk_lo(0xFFFF);
                batch<float, A> cnst65536f(65536.0f);

                auto v_lo = batch_cast<int32_t>(v & msk_lo); /* extract the 16 lowest significant bits of self                             */
                auto v_hi = batch_cast<int32_t>(v >> 16); /* 16 most significant bits of v                                                 */
                auto v_lo_flt = batch_cast<float>(v_lo); /* No rounding                                                                */
                auto v_hi_flt = batch_cast<float>(v_hi); /* No rounding                                                                */
                v_hi_flt = cnst65536f * v_hi_flt; /* No rounding                                                            */
                return v_hi_flt + v_lo_flt; /* Rounding may occur here, mul and add may fuse to fma for haswell and newer   */
            }

            // Provide a generic float -> uint32_t cast only if we have a
            // non-generic float -> int32_t fast_cast
            template <class A, class _ = decltype(fast_cast(std::declval<batch<float, A> const&>(), std::declval<batch<int32_t, A> const&>(), A {}))>
            XSIMD_INLINE batch<uint32_t, A> fast_cast(batch<float, A> const& v, batch<uint32_t, A> const&, requires_arch<generic>) noexcept
            {
                auto is_large = v >= batch<float, A>(1u << 31);
                auto small = bitwise_cast<float>(batch_cast<int32_t>(v));
                auto large = bitwise_cast<float>(
                    batch_cast<int32_t>(v - batch<float, A>(1u << 31))
                    ^ batch<int32_t, A>(1u << 31));
                return bitwise_cast<uint32_t>(select(is_large, large, small));
            }
        }

        namespace detail
        {
            // Generic conversion handling machinery. Each architecture must define
            // conversion function when such conversions exits in the form of
            // intrinsic. Then we use that information to automatically decide whether
            // to use scalar or vector conversion when doing load / store / batch_cast
            struct with_fast_conversion
            {
            };
            struct with_slow_conversion
            {
            };

            template <class A, class From, class To, class = void>
            struct conversion_type_impl
            {
                using type = with_slow_conversion;
            };

            using xsimd::detail::void_t;

            template <class A, class From, class To>
            struct conversion_type_impl<A, From, To,
                                        void_t<decltype(fast_cast(std::declval<const batch<From, A>&>(),
                                                                  std::declval<const batch<To, A>&>(),
                                                                  std::declval<const A&>()))>>
            {
                using type = with_fast_conversion;
            };

            template <class A, class From, class To>
            using conversion_type = typename conversion_type_impl<A, From, To>::type;
        }

        namespace detail
        {
            /* origin: boost/simdfunction/horn.hpp*/
            /*
             * ====================================================
             * copyright 2016 NumScale SAS
             *
             * Distributed under the Boost Software License, Version 1.0.
             * (See copy at http://boost.org/LICENSE_1_0.txt)
             * ====================================================
             */
            template <class B, uint64_t c>
            XSIMD_INLINE B coef() noexcept
            {
                using value_type = typename B::value_type;
                return B(bit_cast<value_type>(as_unsigned_integer_t<value_type>(c)));
            }
            template <class B>
            XSIMD_INLINE B horner(const B&) noexcept
            {
                return B(typename B::value_type(0.));
            }

            template <class B, uint64_t c0>
            XSIMD_INLINE B horner(const B&) noexcept
            {
                return coef<B, c0>();
            }

            template <class B, uint64_t c0, uint64_t c1, uint64_t... args>
            XSIMD_INLINE B horner(const B& self) noexcept
            {
                return fma(self, horner<B, c1, args...>(self), coef<B, c0>());
            }

            /* origin: boost/simdfunction/horn1.hpp*/
            /*
             * ====================================================
             * copyright 2016 NumScale SAS
             *
             * Distributed under the Boost Software License, Version 1.0.
             * (See copy at http://boost.org/LICENSE_1_0.txt)
             * ====================================================
             */
            template <class B>
            XSIMD_INLINE B horner1(const B&) noexcept
            {
                return B(1.);
            }

            template <class B, uint64_t c0>
            XSIMD_INLINE B horner1(const B& x) noexcept
            {
                return x + detail::coef<B, c0>();
            }

            template <class B, uint64_t c0, uint64_t c1, uint64_t... args>
            XSIMD_INLINE B horner1(const B& x) noexcept
            {
                return fma(x, horner1<B, c1, args...>(x), detail::coef<B, c0>());
            }
        }

    }

}

#endif
