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

#ifndef XSIMD_NEON64_HPP
#define XSIMD_NEON64_HPP

#include <complex>
#include <cstddef>
#include <tuple>

#include "../types/xsimd_neon64_register.hpp"
#include "../types/xsimd_utils.hpp"

namespace xsimd
{
    template <typename T, class A, bool... Values>
    struct batch_bool_constant;

    namespace kernel
    {
        using namespace types;

        /*******
         * all *
         *******/

        template <class A, class T, detail::enable_sized_t<T, 4> = 0>
        XSIMD_INLINE bool all(batch_bool<T, A> const& arg, requires_arch<neon64>) noexcept
        {
            return vminvq_u32(arg) == ~0U;
        }

        template <class A, class T, detail::enable_sized_t<T, 1> = 0>
        XSIMD_INLINE bool all(batch_bool<T, A> const& arg, requires_arch<neon64>) noexcept
        {
            return all(batch_bool<uint32_t, A>(vreinterpretq_u32_u8(arg)), neon64 {});
        }

        template <class A, class T, detail::enable_sized_t<T, 2> = 0>
        XSIMD_INLINE bool all(batch_bool<T, A> const& arg, requires_arch<neon64>) noexcept
        {
            return all(batch_bool<uint32_t, A>(vreinterpretq_u32_u16(arg)), neon64 {});
        }

        template <class A, class T, detail::enable_sized_t<T, 8> = 0>
        XSIMD_INLINE bool all(batch_bool<T, A> const& arg, requires_arch<neon64>) noexcept
        {
            return all(batch_bool<uint32_t, A>(vreinterpretq_u32_u64(arg)), neon64 {});
        }

        /*******
         * any *
         *******/

        template <class A, class T, detail::enable_sized_t<T, 4> = 0>
        XSIMD_INLINE bool any(batch_bool<T, A> const& arg, requires_arch<neon64>) noexcept
        {
            return vmaxvq_u32(arg) != 0;
        }

        template <class A, class T, detail::enable_sized_t<T, 1> = 0>
        XSIMD_INLINE bool any(batch_bool<T, A> const& arg, requires_arch<neon64>) noexcept
        {
            return any(batch_bool<uint32_t, A>(vreinterpretq_u32_u8(arg)), neon64 {});
        }

        template <class A, class T, detail::enable_sized_t<T, 2> = 0>
        XSIMD_INLINE bool any(batch_bool<T, A> const& arg, requires_arch<neon64>) noexcept
        {
            return any(batch_bool<uint32_t, A>(vreinterpretq_u32_u16(arg)), neon64 {});
        }

        template <class A, class T, detail::enable_sized_t<T, 8> = 0>
        XSIMD_INLINE bool any(batch_bool<T, A> const& arg, requires_arch<neon64>) noexcept
        {
            return any(batch_bool<uint32_t, A>(vreinterpretq_u32_u64(arg)), neon64 {});
        }

        /*************
         * broadcast *
         *************/

        // Required to avoid ambiguous call
        template <class A, class T>
        XSIMD_INLINE batch<T, A> broadcast(T val, requires_arch<neon64>) noexcept
        {
            return broadcast<A>(val, neon {});
        }

        template <class A>
        XSIMD_INLINE batch<double, A> broadcast(double val, requires_arch<neon64>) noexcept
        {
            return vdupq_n_f64(val);
        }

        /*******
         * set *
         *******/

        template <class A>
        XSIMD_INLINE batch<double, A> set(batch<double, A> const&, requires_arch<neon64>, double d0, double d1) noexcept
        {
            return float64x2_t { d0, d1 };
        }

        template <class A>
        XSIMD_INLINE batch_bool<double, A> set(batch_bool<double, A> const&, requires_arch<neon64>, bool b0, bool b1) noexcept
        {
            using register_type = typename batch_bool<double, A>::register_type;
            using unsigned_type = as_unsigned_integer_t<double>;
            return register_type { static_cast<unsigned_type>(b0 ? -1LL : 0LL),
                                   static_cast<unsigned_type>(b1 ? -1LL : 0LL) };
        }

        /*************
         * from_bool *
         *************/

        template <class A>
        XSIMD_INLINE batch<double, A> from_bool(batch_bool<double, A> const& arg, requires_arch<neon64>) noexcept
        {
            return vreinterpretq_f64_u64(vandq_u64(arg, vreinterpretq_u64_f64(vdupq_n_f64(1.))));
        }

        /********
         * load *
         ********/
#if defined(__clang__) || defined(__GNUC__)
#define xsimd_aligned_load(inst, type, expr) inst((type)__builtin_assume_aligned(expr, 16))
#elif defined(_MSC_VER)
#define xsimd_aligned_load(inst, type, expr) inst##_ex((type)expr, 128)
#else
#define xsimd_aligned_load(inst, type, expr) inst((type)expr)
#endif

        template <class A>
        XSIMD_INLINE batch<double, A> load_aligned(double const* src, convert<double>, requires_arch<neon64>) noexcept
        {
            return xsimd_aligned_load(vld1q_f64, double*, src);
        }

        template <class A>
        XSIMD_INLINE batch<double, A> load_unaligned(double const* src, convert<double>, requires_arch<neon64>) noexcept
        {
            return vld1q_f64(src);
        }
#undef xsimd_aligned_load

        /*********
         * store *
         *********/

        template <class A>
        XSIMD_INLINE void store_aligned(double* dst, batch<double, A> const& src, requires_arch<neon64>) noexcept
        {
            vst1q_f64(dst, src);
        }

        template <class A>
        XSIMD_INLINE void store_unaligned(double* dst, batch<double, A> const& src, requires_arch<neon64>) noexcept
        {
            return store_aligned<A>(dst, src, A {});
        }

        /****************
         * load_complex *
         ****************/

        template <class A>
        XSIMD_INLINE batch<std::complex<double>, A> load_complex_aligned(std::complex<double> const* mem, convert<std::complex<double>>, requires_arch<neon64>) noexcept
        {
            using real_batch = batch<double, A>;
            const double* buf = reinterpret_cast<const double*>(mem);
            float64x2x2_t tmp = vld2q_f64(buf);
            real_batch real = tmp.val[0],
                       imag = tmp.val[1];
            return batch<std::complex<double>, A> { real, imag };
        }

        template <class A>
        XSIMD_INLINE batch<std::complex<double>, A> load_complex_unaligned(std::complex<double> const* mem, convert<std::complex<double>> cvt, requires_arch<neon64>) noexcept
        {
            return load_complex_aligned<A>(mem, cvt, A {});
        }

        /*****************
         * store_complex *
         *****************/

        template <class A>
        XSIMD_INLINE void store_complex_aligned(std::complex<double>* dst, batch<std::complex<double>, A> const& src, requires_arch<neon64>) noexcept
        {
            float64x2x2_t tmp;
            tmp.val[0] = src.real();
            tmp.val[1] = src.imag();
            double* buf = reinterpret_cast<double*>(dst);
            vst2q_f64(buf, tmp);
        }

        template <class A>
        XSIMD_INLINE void store_complex_unaligned(std::complex<double>* dst, batch<std::complex<double>, A> const& src, requires_arch<neon64>) noexcept
        {
            store_complex_aligned(dst, src, A {});
        }

        /*******
         * neg *
         *******/

        template <class A, class T, detail::enable_sized_unsigned_t<T, 8> = 0>
        XSIMD_INLINE batch<T, A> neg(batch<T, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vreinterpretq_u64_s64(vnegq_s64(vreinterpretq_s64_u64(rhs)));
        }

        template <class A, class T, detail::enable_sized_signed_t<T, 8> = 0>
        XSIMD_INLINE batch<T, A> neg(batch<T, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vnegq_s64(rhs);
        }

        template <class A>
        XSIMD_INLINE batch<double, A> neg(batch<double, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vnegq_f64(rhs);
        }

        /*******
         * add *
         *******/

        template <class A>
        XSIMD_INLINE batch<double, A> add(batch<double, A> const& lhs, batch<double, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vaddq_f64(lhs, rhs);
        }

        /********
         * sadd *
         ********/

        template <class A>
        XSIMD_INLINE batch<double, A> sadd(batch<double, A> const& lhs, batch<double, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return add(lhs, rhs, neon64 {});
        }

        /*******
         * sub *
         *******/

        template <class A>
        XSIMD_INLINE batch<double, A> sub(batch<double, A> const& lhs, batch<double, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vsubq_f64(lhs, rhs);
        }

        /********
         * ssub *
         ********/

        template <class A>
        XSIMD_INLINE batch<double, A> ssub(batch<double, A> const& lhs, batch<double, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return sub(lhs, rhs, neon64 {});
        }

        /*******
         * mul *
         *******/

        template <class A>
        XSIMD_INLINE batch<double, A> mul(batch<double, A> const& lhs, batch<double, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vmulq_f64(lhs, rhs);
        }

        /*******
         * div *
         *******/

#if defined(XSIMD_FAST_INTEGER_DIVISION)
        template <class A, class T, detail::enable_sized_unsigned_t<T, 8> = 0>
        XSIMD_INLINE batch<T, A> div(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vcvtq_u64_f64(vcvtq_f64_u64(lhs) / vcvtq_f64_u64(rhs));
        }

        template <class A, class T, detail::enable_sized_signed_t<T, 8> = 0>
        XSIMD_INLINE batch<T, A> div(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vcvtq_s64_f64(vcvtq_f64_s64(lhs) / vcvtq_f64_s64(rhs));
        }
#endif
        template <class A>
        XSIMD_INLINE batch<double, A> div(batch<double, A> const& lhs, batch<double, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vdivq_f64(lhs, rhs);
        }

        /******
         * eq *
         ******/

        template <class A, class T, detail::enable_sized_unsigned_t<T, 8> = 0>
        XSIMD_INLINE batch_bool<T, A> eq(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vceqq_u64(lhs, rhs);
        }

        template <class A, class T, detail::enable_sized_signed_t<T, 8> = 0>
        XSIMD_INLINE batch_bool<T, A> eq(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vceqq_s64(lhs, rhs);
        }

        template <class A>
        XSIMD_INLINE batch_bool<double, A> eq(batch<double, A> const& lhs, batch<double, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vceqq_f64(lhs, rhs);
        }

        template <class A, class T, detail::enable_sized_unsigned_t<T, 8> = 0>
        XSIMD_INLINE batch_bool<T, A> eq(batch_bool<T, A> const& lhs, batch_bool<T, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vceqq_u64(lhs, rhs);
        }

        template <class A, class T, detail::enable_sized_signed_t<T, 8> = 0>
        XSIMD_INLINE batch_bool<T, A> eq(batch_bool<T, A> const& lhs, batch_bool<T, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vceqq_u64(lhs, rhs);
        }

        template <class A>
        XSIMD_INLINE batch_bool<double, A> eq(batch_bool<double, A> const& lhs, batch_bool<double, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vceqq_u64(lhs, rhs);
        }

        /*************
         * fast_cast *
         *************/
        namespace detail
        {
            template <class A>
            XSIMD_INLINE batch<double, A> fast_cast(batch<int64_t, A> const& x, batch<double, A> const&, requires_arch<neon64>) noexcept
            {
                return vcvtq_f64_s64(x);
            }

            template <class A>
            XSIMD_INLINE batch<double, A> fast_cast(batch<uint64_t, A> const& x, batch<double, A> const&, requires_arch<neon64>) noexcept
            {
                return vcvtq_f64_u64(x);
            }

            template <class A>
            XSIMD_INLINE batch<int64_t, A> fast_cast(batch<double, A> const& x, batch<int64_t, A> const&, requires_arch<neon64>) noexcept
            {
                return vcvtq_s64_f64(x);
            }

            template <class A>
            XSIMD_INLINE batch<uint64_t, A> fast_cast(batch<double, A> const& x, batch<uint64_t, A> const&, requires_arch<neon64>) noexcept
            {
                return vcvtq_u64_f64(x);
            }

        }

        /******
         * lt *
         ******/

        template <class A, class T, detail::enable_sized_unsigned_t<T, 8> = 0>
        XSIMD_INLINE batch_bool<T, A> lt(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vcltq_u64(lhs, rhs);
        }

        template <class A, class T, detail::enable_sized_signed_t<T, 8> = 0>
        XSIMD_INLINE batch_bool<T, A> lt(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vcltq_s64(lhs, rhs);
        }

        template <class A>
        XSIMD_INLINE batch_bool<double, A> lt(batch<double, A> const& lhs, batch<double, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vcltq_f64(lhs, rhs);
        }

        /******
         * le *
         ******/

        template <class A, class T, detail::enable_sized_unsigned_t<T, 8> = 0>
        XSIMD_INLINE batch_bool<T, A> le(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vcleq_u64(lhs, rhs);
        }

        template <class A, class T, detail::enable_sized_signed_t<T, 8> = 0>
        XSIMD_INLINE batch_bool<T, A> le(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vcleq_s64(lhs, rhs);
        }

        template <class A>
        XSIMD_INLINE batch_bool<double, A> le(batch<double, A> const& lhs, batch<double, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vcleq_f64(lhs, rhs);
        }

        /******
         * gt *
         ******/

        template <class A, class T, detail::enable_sized_unsigned_t<T, 8> = 0>
        XSIMD_INLINE batch_bool<T, A> gt(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vcgtq_u64(lhs, rhs);
        }

        template <class A, class T, detail::enable_sized_signed_t<T, 8> = 0>
        XSIMD_INLINE batch_bool<T, A> gt(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vcgtq_s64(lhs, rhs);
        }

        template <class A>
        XSIMD_INLINE batch_bool<double, A> gt(batch<double, A> const& lhs, batch<double, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vcgtq_f64(lhs, rhs);
        }

        /******
         * ge *
         ******/

        template <class A, class T, detail::enable_sized_unsigned_t<T, 8> = 0>
        XSIMD_INLINE batch_bool<T, A> ge(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vcgeq_u64(lhs, rhs);
        }

        template <class A, class T, detail::enable_sized_signed_t<T, 8> = 0>
        XSIMD_INLINE batch_bool<T, A> ge(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vcgeq_s64(lhs, rhs);
        }

        template <class A>
        XSIMD_INLINE batch_bool<double, A> ge(batch<double, A> const& lhs, batch<double, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vcgeq_f64(lhs, rhs);
        }

        /*******************
         * batch_bool_cast *
         *******************/

        template <class A, class T_out, class T_in>
        XSIMD_INLINE batch_bool<T_out, A> batch_bool_cast(batch_bool<T_in, A> const& self, batch_bool<T_out, A> const&, requires_arch<neon64>) noexcept
        {
            using register_type = typename batch_bool<T_out, A>::register_type;
            return register_type(self);
        }

        /***************
         * bitwise_and *
         ***************/

        template <class A>
        XSIMD_INLINE batch<double, A> bitwise_and(batch<double, A> const& lhs, batch<double, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vreinterpretq_f64_u64(vandq_u64(vreinterpretq_u64_f64(lhs),
                                                   vreinterpretq_u64_f64(rhs)));
        }

        template <class A>
        XSIMD_INLINE batch_bool<double, A> bitwise_and(batch_bool<double, A> const& lhs, batch_bool<double, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vandq_u64(lhs, rhs);
        }

        /**************
         * bitwise_or *
         **************/

        template <class A>
        XSIMD_INLINE batch<double, A> bitwise_or(batch<double, A> const& lhs, batch<double, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vreinterpretq_f64_u64(vorrq_u64(vreinterpretq_u64_f64(lhs),
                                                   vreinterpretq_u64_f64(rhs)));
        }

        template <class A>
        XSIMD_INLINE batch_bool<double, A> bitwise_or(batch_bool<double, A> const& lhs, batch_bool<double, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vorrq_u64(lhs, rhs);
        }

        /***************
         * bitwise_xor *
         ***************/

        template <class A>
        XSIMD_INLINE batch<double, A> bitwise_xor(batch<double, A> const& lhs, batch<double, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vreinterpretq_f64_u64(veorq_u64(vreinterpretq_u64_f64(lhs),
                                                   vreinterpretq_u64_f64(rhs)));
        }

        template <class A>
        XSIMD_INLINE batch_bool<double, A> bitwise_xor(batch_bool<double, A> const& lhs, batch_bool<double, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return veorq_u64(lhs, rhs);
        }

        /*******
         * neq *
         *******/

        template <class A>
        XSIMD_INLINE batch_bool<double, A> neq(batch_bool<double, A> const& lhs, batch_bool<double, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return bitwise_xor(lhs, rhs, A {});
        }

        /***************
         * bitwise_not *
         ***************/

        template <class A>
        XSIMD_INLINE batch<double, A> bitwise_not(batch<double, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vreinterpretq_f64_u32(vmvnq_u32(vreinterpretq_u32_f64(rhs)));
        }

        template <class A>
        XSIMD_INLINE batch_bool<double, A> bitwise_not(batch_bool<double, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return detail::bitwise_not_u64(rhs);
        }

        /******************
         * bitwise_andnot *
         ******************/

        template <class A>
        XSIMD_INLINE batch<double, A> bitwise_andnot(batch<double, A> const& lhs, batch<double, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vreinterpretq_f64_u64(vbicq_u64(vreinterpretq_u64_f64(lhs),
                                                   vreinterpretq_u64_f64(rhs)));
        }

        template <class A>
        XSIMD_INLINE batch_bool<double, A> bitwise_andnot(batch_bool<double, A> const& lhs, batch_bool<double, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vbicq_u64(lhs, rhs);
        }

        /*******
         * min *
         *******/

        template <class A>
        XSIMD_INLINE batch<double, A> min(batch<double, A> const& lhs, batch<double, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vminq_f64(lhs, rhs);
        }

        /*******
         * max *
         *******/

        template <class A>
        XSIMD_INLINE batch<double, A> max(batch<double, A> const& lhs, batch<double, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vmaxq_f64(lhs, rhs);
        }

        /*******
         * abs *
         *******/

        template <class A, class T, detail::enable_sized_unsigned_t<T, 8> = 0>
        XSIMD_INLINE batch<T, A> abs(batch<T, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return rhs;
        }

        template <class A, class T, detail::enable_sized_signed_t<T, 8> = 0>
        XSIMD_INLINE batch<T, A> abs(batch<T, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vabsq_s64(rhs);
        }

        template <class A>
        XSIMD_INLINE batch<double, A> abs(batch<double, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vabsq_f64(rhs);
        }

        template <class A>
        XSIMD_INLINE batch<int32_t, A> nearbyint_as_int(batch<float, A> const& self,
                                                        requires_arch<neon64>) noexcept
        {
            return vcvtnq_s32_f32(self);
        }

#if !defined(__GNUC__)
        template <class A>
        XSIMD_INLINE batch<int64_t, A> nearbyint_as_int(batch<double, A> const& self,
                                                        requires_arch<neon64>) noexcept
        {
            return vcvtnq_s64_f64(self);
        }
#endif

        /**************
         * reciprocal *
         **************/

        template <class A>
        XSIMD_INLINE batch<double, A>
        reciprocal(const batch<double, A>& x,
                   kernel::requires_arch<neon64>) noexcept
        {
            return vrecpeq_f64(x);
        }

        /********
         * rsqrt *
         ********/

        template <class A>
        XSIMD_INLINE batch<double, A> rsqrt(batch<double, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vrsqrteq_f64(rhs);
        }

        /********
         * sqrt *
         ********/

        template <class A>
        XSIMD_INLINE batch<double, A> sqrt(batch<double, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vsqrtq_f64(rhs);
        }

        /********************
         * Fused operations *
         ********************/

#ifdef __ARM_FEATURE_FMA
        template <class A>
        XSIMD_INLINE batch<double, A> fma(batch<double, A> const& x, batch<double, A> const& y, batch<double, A> const& z, requires_arch<neon64>) noexcept
        {
            return vfmaq_f64(z, x, y);
        }

        template <class A>
        XSIMD_INLINE batch<double, A> fms(batch<double, A> const& x, batch<double, A> const& y, batch<double, A> const& z, requires_arch<neon64>) noexcept
        {
            return vfmaq_f64(-z, x, y);
        }
#endif

        /*********
         * haddp *
         *********/

        template <class A>
        XSIMD_INLINE batch<double, A> haddp(const batch<double, A>* row, requires_arch<neon64>) noexcept
        {
            return vpaddq_f64(row[0], row[1]);
        }

        /**********
         * insert *
         **********/

        template <class A, size_t I>
        XSIMD_INLINE batch<double, A> insert(batch<double, A> const& self, double val, index<I>, requires_arch<neon64>) noexcept
        {
            return vsetq_lane_f64(val, self, I);
        }

        /******************
         * reducer macros *
         ******************/

        // Wrap reducer intrinsics so we can pass them as function pointers
        // - OP: intrinsics name prefix, e.g., vorrq

#define WRAP_REDUCER_INT_EXCLUDING_64(OP)                     \
    namespace wrap                                            \
    {                                                         \
        XSIMD_INLINE uint8_t OP##_u8(uint8x16_t a) noexcept   \
        {                                                     \
            return ::OP##_u8(a);                              \
        }                                                     \
        XSIMD_INLINE int8_t OP##_s8(int8x16_t a) noexcept     \
        {                                                     \
            return ::OP##_s8(a);                              \
        }                                                     \
        XSIMD_INLINE uint16_t OP##_u16(uint16x8_t a) noexcept \
        {                                                     \
            return ::OP##_u16(a);                             \
        }                                                     \
        XSIMD_INLINE int16_t OP##_s16(int16x8_t a) noexcept   \
        {                                                     \
            return ::OP##_s16(a);                             \
        }                                                     \
        XSIMD_INLINE uint32_t OP##_u32(uint32x4_t a) noexcept \
        {                                                     \
            return ::OP##_u32(a);                             \
        }                                                     \
        XSIMD_INLINE int32_t OP##_s32(int32x4_t a) noexcept   \
        {                                                     \
            return ::OP##_s32(a);                             \
        }                                                     \
    }

#define WRAP_REDUCER_INT(OP)                                  \
    WRAP_REDUCER_INT_EXCLUDING_64(OP)                         \
    namespace wrap                                            \
    {                                                         \
        XSIMD_INLINE uint64_t OP##_u64(uint64x2_t a) noexcept \
        {                                                     \
            return ::OP##_u64(a);                             \
        }                                                     \
        XSIMD_INLINE int64_t OP##_s64(int64x2_t a) noexcept   \
        {                                                     \
            return ::OP##_s64(a);                             \
        }                                                     \
    }

#define WRAP_REDUCER_FLOAT(OP)                               \
    namespace wrap                                           \
    {                                                        \
        XSIMD_INLINE float OP##_f32(float32x4_t a) noexcept  \
        {                                                    \
            return ::OP##_f32(a);                            \
        }                                                    \
        XSIMD_INLINE double OP##_f64(float64x2_t a) noexcept \
        {                                                    \
            return ::OP##_f64(a);                            \
        }                                                    \
    }

        namespace detail
        {
            template <class R>
            struct reducer_return_type_impl;

            template <>
            struct reducer_return_type_impl<uint8x16_t>
            {
                using type = uint8_t;
            };

            template <>
            struct reducer_return_type_impl<int8x16_t>
            {
                using type = int8_t;
            };

            template <>
            struct reducer_return_type_impl<uint16x8_t>
            {
                using type = uint16_t;
            };

            template <>
            struct reducer_return_type_impl<int16x8_t>
            {
                using type = int16_t;
            };

            template <>
            struct reducer_return_type_impl<uint32x4_t>
            {
                using type = uint32_t;
            };

            template <>
            struct reducer_return_type_impl<int32x4_t>
            {
                using type = int32_t;
            };

            template <>
            struct reducer_return_type_impl<uint64x2_t>
            {
                using type = uint64_t;
            };

            template <>
            struct reducer_return_type_impl<int64x2_t>
            {
                using type = int64_t;
            };

            template <>
            struct reducer_return_type_impl<float32x4_t>
            {
                using type = float;
            };

            template <>
            struct reducer_return_type_impl<float64x2_t>
            {
                using type = double;
            };

            template <class R>
            using reducer_return_type = typename reducer_return_type_impl<R>::type;

            template <class... T>
            struct neon_reducer_dispatcher_impl : neon_dispatcher_base<reducer_return_type, T...>
            {
            };

            using neon_reducer_dispatcher = neon_reducer_dispatcher_impl<uint8x16_t, int8x16_t,
                                                                         uint16x8_t, int16x8_t,
                                                                         uint32x4_t, int32x4_t,
                                                                         uint64x2_t, int64x2_t,
                                                                         float32x4_t, float64x2_t>;
            template <class T>
            using enable_neon64_type_t = typename std::enable_if<std::is_integral<T>::value || std::is_same<T, float>::value || std::is_same<T, double>::value,
                                                                 int>::type;
        }

        /**************
         * reduce_add *
         **************/

        WRAP_REDUCER_INT(vaddvq)
        WRAP_REDUCER_FLOAT(vaddvq)

        template <class A, class T, detail::enable_neon64_type_t<T> = 0>
        XSIMD_INLINE typename batch<T, A>::value_type reduce_add(batch<T, A> const& arg, requires_arch<neon64>) noexcept
        {
            using register_type = typename batch<T, A>::register_type;
            const detail::neon_reducer_dispatcher::unary dispatcher = {
                std::make_tuple(wrap::vaddvq_u8, wrap::vaddvq_s8, wrap::vaddvq_u16, wrap::vaddvq_s16,
                                wrap::vaddvq_u32, wrap::vaddvq_s32, wrap::vaddvq_u64, wrap::vaddvq_s64,
                                wrap::vaddvq_f32, wrap::vaddvq_f64)
            };
            return dispatcher.apply(register_type(arg));
        }

        /**************
         * reduce_max *
         **************/

        WRAP_REDUCER_INT_EXCLUDING_64(vmaxvq)
        WRAP_REDUCER_FLOAT(vmaxvq)

        namespace wrap
        {
            XSIMD_INLINE uint64_t vmaxvq_u64(uint64x2_t a) noexcept
            {
                return std::max(vdupd_laneq_u64(a, 0), vdupd_laneq_u64(a, 1));
            }

            XSIMD_INLINE int64_t vmaxvq_s64(int64x2_t a) noexcept
            {
                return std::max(vdupd_laneq_s64(a, 0), vdupd_laneq_s64(a, 1));
            }
        }

        template <class A, class T, detail::enable_neon64_type_t<T> = 0>
        XSIMD_INLINE typename batch<T, A>::value_type reduce_max(batch<T, A> const& arg, requires_arch<neon64>) noexcept
        {
            using register_type = typename batch<T, A>::register_type;
            const detail::neon_reducer_dispatcher::unary dispatcher = {
                std::make_tuple(wrap::vmaxvq_u8, wrap::vmaxvq_s8, wrap::vmaxvq_u16, wrap::vmaxvq_s16,
                                wrap::vmaxvq_u32, wrap::vmaxvq_s32, wrap::vmaxvq_u64, wrap::vmaxvq_s64,
                                wrap::vmaxvq_f32, wrap::vmaxvq_f64)
            };
            return dispatcher.apply(register_type(arg));
        }

        /**************
         * reduce_min *
         **************/

        WRAP_REDUCER_INT_EXCLUDING_64(vminvq)
        WRAP_REDUCER_FLOAT(vminvq)

        namespace wrap
        {
            XSIMD_INLINE uint64_t vminvq_u64(uint64x2_t a) noexcept
            {
                return std::min(vdupd_laneq_u64(a, 0), vdupd_laneq_u64(a, 1));
            }

            XSIMD_INLINE int64_t vminvq_s64(int64x2_t a) noexcept
            {
                return std::min(vdupd_laneq_s64(a, 0), vdupd_laneq_s64(a, 1));
            }
        }

        template <class A, class T, detail::enable_neon64_type_t<T> = 0>
        XSIMD_INLINE typename batch<T, A>::value_type reduce_min(batch<T, A> const& arg, requires_arch<neon64>) noexcept
        {
            using register_type = typename batch<T, A>::register_type;
            const detail::neon_reducer_dispatcher::unary dispatcher = {
                std::make_tuple(wrap::vminvq_u8, wrap::vminvq_s8, wrap::vminvq_u16, wrap::vminvq_s16,
                                wrap::vminvq_u32, wrap::vminvq_s32, wrap::vminvq_u64, wrap::vminvq_s64,
                                wrap::vminvq_f32, wrap::vminvq_f64)
            };
            return dispatcher.apply(register_type(arg));
        }

#undef WRAP_REDUCER_INT_EXCLUDING_64
#undef WRAP_REDUCER_INT
#undef WRAP_REDUCER_FLOAT

        /**********
         * select *
         **********/

        template <class A>
        XSIMD_INLINE batch<double, A> select(batch_bool<double, A> const& cond, batch<double, A> const& a, batch<double, A> const& b, requires_arch<neon64>) noexcept
        {
            return vbslq_f64(cond, a, b);
        }

        template <class A, bool... b>
        XSIMD_INLINE batch<double, A> select(batch_bool_constant<double, A, b...> const&,
                                             batch<double, A> const& true_br,
                                             batch<double, A> const& false_br,
                                             requires_arch<neon64>) noexcept
        {
            return select(batch_bool<double, A> { b... }, true_br, false_br, neon64 {});
        }
        /**********
         * zip_lo *
         **********/
        template <class A, class T, detail::enable_sized_unsigned_t<T, 1> = 0>
        XSIMD_INLINE batch<T, A> zip_lo(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vzip1q_u8(lhs, rhs);
        }

        template <class A, class T, detail::enable_sized_signed_t<T, 1> = 0>
        XSIMD_INLINE batch<T, A> zip_lo(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vzip1q_s8(lhs, rhs);
        }

        template <class A, class T, detail::enable_sized_unsigned_t<T, 2> = 0>
        XSIMD_INLINE batch<T, A> zip_lo(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vzip1q_u16(lhs, rhs);
        }

        template <class A, class T, detail::enable_sized_signed_t<T, 2> = 0>
        XSIMD_INLINE batch<T, A> zip_lo(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vzip1q_s16(lhs, rhs);
        }

        template <class A, class T, detail::enable_sized_unsigned_t<T, 4> = 0>
        XSIMD_INLINE batch<T, A> zip_lo(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vzip1q_u32(lhs, rhs);
        }

        template <class A, class T, detail::enable_sized_signed_t<T, 4> = 0>
        XSIMD_INLINE batch<T, A> zip_lo(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vzip1q_s32(lhs, rhs);
        }

        template <class A, class T, detail::enable_sized_unsigned_t<T, 8> = 0>
        XSIMD_INLINE batch<T, A> zip_lo(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vzip1q_u64(lhs, rhs);
        }

        template <class A, class T, detail::enable_sized_signed_t<T, 8> = 0>
        XSIMD_INLINE batch<T, A> zip_lo(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vzip1q_s64(lhs, rhs);
        }

        template <class A>
        XSIMD_INLINE batch<float, A> zip_lo(batch<float, A> const& lhs, batch<float, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vzip1q_f32(lhs, rhs);
        }

        template <class A>
        XSIMD_INLINE batch<double, A> zip_lo(batch<double, A> const& lhs, batch<double, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vzip1q_f64(lhs, rhs);
        }

        /**********
         * zip_hi *
         **********/

        template <class A, class T, detail::enable_sized_unsigned_t<T, 1> = 0>
        XSIMD_INLINE batch<T, A> zip_hi(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vzip2q_u8(lhs, rhs);
        }

        template <class A, class T, detail::enable_sized_signed_t<T, 1> = 0>
        XSIMD_INLINE batch<T, A> zip_hi(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vzip2q_s8(lhs, rhs);
        }

        template <class A, class T, detail::enable_sized_unsigned_t<T, 2> = 0>
        XSIMD_INLINE batch<T, A> zip_hi(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vzip2q_u16(lhs, rhs);
        }

        template <class A, class T, detail::enable_sized_signed_t<T, 2> = 0>
        XSIMD_INLINE batch<T, A> zip_hi(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vzip2q_s16(lhs, rhs);
        }

        template <class A, class T, detail::enable_sized_unsigned_t<T, 4> = 0>
        XSIMD_INLINE batch<T, A> zip_hi(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vzip2q_u32(lhs, rhs);
        }

        template <class A, class T, detail::enable_sized_signed_t<T, 4> = 0>
        XSIMD_INLINE batch<T, A> zip_hi(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vzip2q_s32(lhs, rhs);
        }

        template <class A, class T, detail::enable_sized_unsigned_t<T, 8> = 0>
        XSIMD_INLINE batch<T, A> zip_hi(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vzip2q_u64(lhs, rhs);
        }

        template <class A, class T, detail::enable_sized_signed_t<T, 8> = 0>
        XSIMD_INLINE batch<T, A> zip_hi(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vzip2q_s64(lhs, rhs);
        }

        template <class A>
        XSIMD_INLINE batch<float, A> zip_hi(batch<float, A> const& lhs, batch<float, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vzip2q_f32(lhs, rhs);
        }

        template <class A>
        XSIMD_INLINE batch<double, A> zip_hi(batch<double, A> const& lhs, batch<double, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vzip2q_f64(lhs, rhs);
        }

        /****************
         * extract_pair *
         ****************/

        namespace detail
        {
            template <class A, size_t I, size_t... Is>
            XSIMD_INLINE batch<double, A> extract_pair(batch<double, A> const& lhs, batch<double, A> const& rhs, std::size_t n,
                                                       ::xsimd::detail::index_sequence<I, Is...>) noexcept
            {
                if (n == I)
                {
                    return vextq_f64(rhs, lhs, I);
                }
                else
                {
                    return extract_pair(lhs, rhs, n, ::xsimd::detail::index_sequence<Is...>());
                }
            }
        }

        template <class A>
        XSIMD_INLINE batch<double, A> extract_pair(batch<double, A> const& lhs, batch<double, A> const& rhs, std::size_t n, requires_arch<neon64>) noexcept
        {
            constexpr std::size_t size = batch<double, A>::size;
            assert(n < size && "index in bounds");
            return detail::extract_pair(lhs, rhs, n, ::xsimd::detail::make_index_sequence<size>());
        }

        /******************
         * bitwise_rshift *
         ******************/

        template <class A, class T, detail::enable_sized_unsigned_t<T, 8> = 0>
        XSIMD_INLINE batch<T, A> bitwise_rshift(batch<T, A> const& lhs, int n, requires_arch<neon64>) noexcept
        {
            return bitwise_rshift<A>(lhs, n, neon {});
        }

        template <class A, class T, detail::enable_sized_unsigned_t<T, 8> = 0>
        XSIMD_INLINE batch<T, A> bitwise_rshift(batch<T, A> const& lhs, batch<as_signed_integer_t<T>, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vshlq_u64(lhs, vnegq_s64(rhs));
        }

        template <class A, class T, detail::enable_sized_signed_t<T, 8> = 0>
        XSIMD_INLINE batch<T, A> bitwise_rshift(batch<T, A> const& lhs, int n, requires_arch<neon64>) noexcept
        {
            return bitwise_rshift<A>(lhs, n, neon {});
        }

        template <class A, class T, detail::enable_sized_signed_t<T, 8> = 0>
        XSIMD_INLINE batch<T, A> bitwise_rshift(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon64>) noexcept
        {
            return vshlq_s64(lhs, vnegq_s64(rhs));
        }

        /****************
         * bitwise_cast *
         ****************/

#define WRAP_CAST(SUFFIX, TYPE)                                                \
    namespace wrap                                                             \
    {                                                                          \
        XSIMD_INLINE float64x2_t vreinterpretq_f64_##SUFFIX(TYPE a) noexcept   \
        {                                                                      \
            return ::vreinterpretq_f64_##SUFFIX(a);                            \
        }                                                                      \
        XSIMD_INLINE TYPE vreinterpretq_##SUFFIX##_f64(float64x2_t a) noexcept \
        {                                                                      \
            return ::vreinterpretq_##SUFFIX##_f64(a);                          \
        }                                                                      \
    }

        WRAP_CAST(u8, uint8x16_t)
        WRAP_CAST(s8, int8x16_t)
        WRAP_CAST(u16, uint16x8_t)
        WRAP_CAST(s16, int16x8_t)
        WRAP_CAST(u32, uint32x4_t)
        WRAP_CAST(s32, int32x4_t)
        WRAP_CAST(u64, uint64x2_t)
        WRAP_CAST(s64, int64x2_t)
        WRAP_CAST(f32, float32x4_t)

#undef WRAP_CAST

        template <class A, class T>
        XSIMD_INLINE batch<double, A> bitwise_cast(batch<T, A> const& arg, batch<double, A> const&, requires_arch<neon64>) noexcept
        {
            using caster_type = detail::bitwise_caster_impl<float64x2_t,
                                                            uint8x16_t, int8x16_t,
                                                            uint16x8_t, int16x8_t,
                                                            uint32x4_t, int32x4_t,
                                                            uint64x2_t, int64x2_t,
                                                            float32x4_t>;
            const caster_type caster = {
                std::make_tuple(wrap::vreinterpretq_f64_u8, wrap::vreinterpretq_f64_s8, wrap::vreinterpretq_f64_u16, wrap::vreinterpretq_f64_s16,
                                wrap::vreinterpretq_f64_u32, wrap::vreinterpretq_f64_s32, wrap::vreinterpretq_f64_u64, wrap::vreinterpretq_f64_s64,
                                wrap::vreinterpretq_f64_f32)
            };
            using register_type = typename batch<T, A>::register_type;
            return caster.apply(register_type(arg));
        }

        namespace detail
        {
            template <class S, class... R>
            struct bitwise_caster_neon64
            {
                using container_type = std::tuple<R (*)(S)...>;
                container_type m_func;

                template <class V>
                V apply(float64x2_t rhs) const
                {
                    using func_type = V (*)(float64x2_t);
                    auto func = xsimd::detail::get<func_type>(m_func);
                    return func(rhs);
                }
            };
        }

        template <class A, class R>
        XSIMD_INLINE batch<R, A> bitwise_cast(batch<double, A> const& arg, batch<R, A> const&, requires_arch<neon64>) noexcept
        {
            using caster_type = detail::bitwise_caster_neon64<float64x2_t,
                                                              uint8x16_t, int8x16_t,
                                                              uint16x8_t, int16x8_t,
                                                              uint32x4_t, int32x4_t,
                                                              uint64x2_t, int64x2_t,
                                                              float32x4_t>;
            const caster_type caster = {
                std::make_tuple(wrap::vreinterpretq_u8_f64, wrap::vreinterpretq_s8_f64, wrap::vreinterpretq_u16_f64, wrap::vreinterpretq_s16_f64,
                                wrap::vreinterpretq_u32_f64, wrap::vreinterpretq_s32_f64, wrap::vreinterpretq_u64_f64, wrap::vreinterpretq_s64_f64,
                                wrap::vreinterpretq_f32_f64)
            };
            using src_register_type = typename batch<double, A>::register_type;
            using dst_register_type = typename batch<R, A>::register_type;
            return caster.apply<dst_register_type>(src_register_type(arg));
        }

        template <class A>
        XSIMD_INLINE batch<double, A> bitwise_cast(batch<double, A> const& arg, batch<double, A> const&, requires_arch<neon64>) noexcept
        {
            return arg;
        }

        /*********
         * isnan *
         *********/

        template <class A>
        XSIMD_INLINE batch_bool<double, A> isnan(batch<double, A> const& arg, requires_arch<neon64>) noexcept
        {
            return !(arg == arg);
        }

        /****************
         * rotate_right *
         ****************/
        template <size_t N, class A>
        XSIMD_INLINE batch<double, A> rotate_right(batch<double, A> const& a, requires_arch<neon64>) noexcept
        {
            return vextq_f64(a, a, N);
        }
    }

    template <typename T, class A, T... Values>
    struct batch_constant;

    namespace kernel
    {
        /*********************
         * swizzle (dynamic) *
         *********************/
        template <class A>
        XSIMD_INLINE batch<uint8_t, A> swizzle(batch<uint8_t, A> const& self, batch<uint8_t, A> idx,
                                               requires_arch<neon64>) noexcept
        {
            return vqtbl1q_u8(self, idx);
        }

        template <class A>
        XSIMD_INLINE batch<int8_t, A> swizzle(batch<int8_t, A> const& self, batch<uint8_t, A> idx,
                                              requires_arch<neon64>) noexcept
        {
            return vqtbl1q_s8(self, idx);
        }

        template <class A>
        XSIMD_INLINE batch<uint16_t, A> swizzle(batch<uint16_t, A> const& self,
                                                batch<uint16_t, A> idx,
                                                requires_arch<neon64>) noexcept
        {
            using batch_type = batch<uint8_t, A>;
            using index_type = batch<uint8_t, A>;
            return vreinterpretq_u16_u8(swizzle(batch_type(vreinterpretq_u8_u16(self)),
                                                index_type(vreinterpretq_u8_u16(idx * 0x0202 + 0x0100)),
                                                neon64 {}));
        }

        template <class A>
        XSIMD_INLINE batch<int16_t, A> swizzle(batch<int16_t, A> const& self,
                                               batch<uint16_t, A> idx,
                                               requires_arch<neon64>) noexcept
        {
            return bitwise_cast<int16_t>(swizzle(bitwise_cast<uint16_t>(self), idx, neon64 {}));
        }

        template <class A>
        XSIMD_INLINE batch<uint32_t, A> swizzle(batch<uint32_t, A> const& self,
                                                batch<uint32_t, A> idx,
                                                requires_arch<neon64>) noexcept
        {
            using batch_type = batch<uint8_t, A>;
            using index_type = batch<uint8_t, A>;
            return vreinterpretq_u32_u8(swizzle(batch_type(vreinterpretq_u8_u32(self)),
                                                index_type(vreinterpretq_u8_u32(idx * 0x04040404 + 0x03020100)),
                                                neon64 {}));
        }

        template <class A>
        XSIMD_INLINE batch<int32_t, A> swizzle(batch<int32_t, A> const& self,
                                               batch<uint32_t, A> idx,
                                               requires_arch<neon64>) noexcept
        {
            return bitwise_cast<int32_t>(swizzle(bitwise_cast<uint32_t>(self), idx, neon64 {}));
        }

        template <class A>
        XSIMD_INLINE batch<uint64_t, A> swizzle(batch<uint64_t, A> const& self,
                                                batch<uint64_t, A> idx,
                                                requires_arch<neon64>) noexcept
        {
            using batch_type = batch<uint8_t, A>;
            using index_type = batch<uint8_t, A>;
            return vreinterpretq_u64_u8(swizzle(batch_type(vreinterpretq_u8_u64(self)),
                                                index_type(vreinterpretq_u8_u64(idx * 0x0808080808080808ull + 0x0706050403020100ull)),
                                                neon64 {}));
        }

        template <class A>
        XSIMD_INLINE batch<int64_t, A> swizzle(batch<int64_t, A> const& self,
                                               batch<uint64_t, A> idx,
                                               requires_arch<neon64>) noexcept
        {
            return bitwise_cast<int64_t>(swizzle(bitwise_cast<uint64_t>(self), idx, neon64 {}));
        }

        template <class A>
        XSIMD_INLINE batch<float, A> swizzle(batch<float, A> const& self,
                                             batch<uint32_t, A> idx,
                                             requires_arch<neon64>) noexcept
        {
            return bitwise_cast<float>(swizzle(bitwise_cast<uint32_t>(self), idx, neon64 {}));
        }

        template <class A>
        XSIMD_INLINE batch<double, A> swizzle(batch<double, A> const& self,
                                              batch<uint64_t, A> idx,
                                              requires_arch<neon64>) noexcept
        {
            return bitwise_cast<double>(swizzle(bitwise_cast<uint64_t>(self), idx, neon64 {}));
        }

        /********************
         * swizzle (static) *
         ********************/

        namespace detail
        {
            using ::xsimd::batch_constant;
            using ::xsimd::detail::integer_sequence;
            using ::xsimd::detail::make_integer_sequence;

            template <class CB1, class CB2, class IS>
            struct index_burst_impl;

            template <typename T1, class A, typename T2, T2... V,
                      T2... incr>
            struct index_burst_impl<batch_constant<T1, A>, batch_constant<T2, A, V...>,
                                    integer_sequence<T2, incr...>>
            {
                using type = batch_constant<T2, A, V...>;
            };

            template <typename T1, class A, T1 V0, T1... V1,
                      typename T2, T2... V2, T2... incr>
            struct index_burst_impl<batch_constant<T1, A, V0, V1...>, batch_constant<T2, A, V2...>,
                                    integer_sequence<T2, incr...>>
            {
                using next_input = batch_constant<T1, A, V1...>;
                using next_output = batch_constant<T2, A, V2..., (V0 + incr)...>;
                using type = typename index_burst_impl<next_input, next_output, integer_sequence<T2, incr...>>::type;
            };

            template <class B, class T>
            struct index_burst;

            template <typename Tp, class A, Tp... V, typename T>
            struct index_burst<batch_constant<Tp, A, V...>, T>
            {
                static constexpr size_t mul = sizeof(Tp) / sizeof(T);
                using input = batch_constant<Tp, A, (mul * V)...>;
                using output = batch_constant<T, A>;
                using type = typename index_burst_impl<input, output, make_integer_sequence<T, mul>>::type;
            };

            template <class B, typename T>
            using index_burst_t = typename index_burst<B, T>::type;

            template <typename T, class B>
            XSIMD_INLINE index_burst_t<B, T> burst_index(B)
            {
                return index_burst_t<B, T>();
            }
        }

        template <class A, uint8_t V0, uint8_t V1, uint8_t V2, uint8_t V3, uint8_t V4, uint8_t V5, uint8_t V6, uint8_t V7,
                  uint8_t V8, uint8_t V9, uint8_t V10, uint8_t V11, uint8_t V12, uint8_t V13, uint8_t V14, uint8_t V15>
        XSIMD_INLINE batch<uint8_t, A> swizzle(batch<uint8_t, A> const& self,
                                               batch_constant<uint8_t, A, V0, V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13, V14, V15> idx,
                                               requires_arch<neon64>) noexcept
        {
            return vqtbl1q_u8(self, batch<uint8_t, A>(idx));
        }

        template <class A, uint8_t V0, uint8_t V1, uint8_t V2, uint8_t V3, uint8_t V4, uint8_t V5, uint8_t V6, uint8_t V7,
                  uint8_t V8, uint8_t V9, uint8_t V10, uint8_t V11, uint8_t V12, uint8_t V13, uint8_t V14, uint8_t V15>
        XSIMD_INLINE batch<int8_t, A> swizzle(batch<int8_t, A> const& self,
                                              batch_constant<uint8_t, A, V0, V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13, V14, V15> idx,
                                              requires_arch<neon64>) noexcept
        {
            return vqtbl1q_s8(self, batch<uint8_t, A>(idx));
        }

        template <class A, uint16_t V0, uint16_t V1, uint16_t V2, uint16_t V3, uint16_t V4, uint16_t V5, uint16_t V6, uint16_t V7>
        XSIMD_INLINE batch<uint16_t, A> swizzle(batch<uint16_t, A> const& self,
                                                batch_constant<uint16_t, A, V0, V1, V2, V3, V4, V5, V6, V7> idx,
                                                requires_arch<neon64>) noexcept
        {
            using batch_type = batch<uint8_t, A>;
            return vreinterpretq_u16_u8(swizzle<A>(batch_type(vreinterpretq_u8_u16(self)), detail::burst_index<uint8_t>(idx), A()));
        }

        template <class A, uint16_t V0, uint16_t V1, uint16_t V2, uint16_t V3, uint16_t V4, uint16_t V5, uint16_t V6, uint16_t V7>
        XSIMD_INLINE batch<int16_t, A> swizzle(batch<int16_t, A> const& self,
                                               batch_constant<uint16_t, A, V0, V1, V2, V3, V4, V5, V6, V7> idx,
                                               requires_arch<neon64>) noexcept
        {
            using batch_type = batch<int8_t, A>;
            return vreinterpretq_s16_s8(swizzle<A>(batch_type(vreinterpretq_s8_s16(self)), detail::burst_index<uint8_t>(idx), A()));
        }

        template <class A, uint32_t V0, uint32_t V1, uint32_t V2, uint32_t V3>
        XSIMD_INLINE batch<uint32_t, A> swizzle(batch<uint32_t, A> const& self,
                                                batch_constant<uint32_t, A, V0, V1, V2, V3> idx,
                                                requires_arch<neon64>) noexcept
        {
            using batch_type = batch<uint8_t, A>;
            return vreinterpretq_u32_u8(swizzle<A>(batch_type(vreinterpretq_u8_u32(self)), detail::burst_index<uint8_t>(idx), A()));
        }

        template <class A, uint32_t V0, uint32_t V1, uint32_t V2, uint32_t V3>
        XSIMD_INLINE batch<int32_t, A> swizzle(batch<int32_t, A> const& self,
                                               batch_constant<uint32_t, A, V0, V1, V2, V3> idx,
                                               requires_arch<neon64>) noexcept
        {
            using batch_type = batch<int8_t, A>;
            return vreinterpretq_s32_s8(swizzle<A>(batch_type(vreinterpretq_s8_s32(self)), detail::burst_index<uint8_t>(idx), A()));
        }

        template <class A, uint64_t V0, uint64_t V1>
        XSIMD_INLINE batch<uint64_t, A> swizzle(batch<uint64_t, A> const& self,
                                                batch_constant<uint64_t, A, V0, V1> idx,
                                                requires_arch<neon64>) noexcept
        {
            using batch_type = batch<uint8_t, A>;
            return vreinterpretq_u64_u8(swizzle<A>(batch_type(vreinterpretq_u8_u64(self)), detail::burst_index<uint8_t>(idx), A()));
        }

        template <class A, uint64_t V0, uint64_t V1>
        XSIMD_INLINE batch<int64_t, A> swizzle(batch<int64_t, A> const& self,
                                               batch_constant<uint64_t, A, V0, V1> idx,
                                               requires_arch<neon64>) noexcept
        {
            using batch_type = batch<int8_t, A>;
            return vreinterpretq_s64_s8(swizzle<A>(batch_type(vreinterpretq_s8_s64(self)), detail::burst_index<uint8_t>(idx), A()));
        }

        template <class A, uint32_t V0, uint32_t V1, uint32_t V2, uint32_t V3>
        XSIMD_INLINE batch<float, A> swizzle(batch<float, A> const& self,
                                             batch_constant<uint32_t, A, V0, V1, V2, V3> idx,
                                             requires_arch<neon64>) noexcept
        {
            using batch_type = batch<uint8_t, A>;
            return vreinterpretq_f32_u8(swizzle<A>(batch_type(vreinterpretq_u8_f32(self)), detail::burst_index<uint8_t>(idx), A()));
        }

        template <class A, uint64_t V0, uint64_t V1>
        XSIMD_INLINE batch<double, A> swizzle(batch<double, A> const& self,
                                              batch_constant<uint64_t, A, V0, V1> idx,
                                              requires_arch<neon64>) noexcept
        {
            using batch_type = batch<uint8_t, A>;
            return vreinterpretq_f64_u8(swizzle<A>(batch_type(vreinterpretq_u8_f64(self)), detail::burst_index<uint8_t>(idx), A()));
        }

        template <class A, uint32_t V0, uint32_t V1, uint32_t V2, uint32_t V3>
        XSIMD_INLINE batch<std::complex<float>, A> swizzle(batch<std::complex<float>, A> const& self,
                                                           batch_constant<uint32_t, A, V0, V1, V2, V3> idx,
                                                           requires_arch<neon64>) noexcept
        {
            return batch<std::complex<float>>(swizzle(self.real(), idx, A()), swizzle(self.imag(), idx, A()));
        }

        template <class A, uint64_t V0, uint64_t V1>
        XSIMD_INLINE batch<std::complex<double>, A> swizzle(batch<std::complex<double>, A> const& self,
                                                            batch_constant<uint64_t, A, V0, V1> idx,
                                                            requires_arch<neon64>) noexcept
        {
            return batch<std::complex<double>>(swizzle(self.real(), idx, A()), swizzle(self.imag(), idx, A()));
        }
    }
}

#endif
