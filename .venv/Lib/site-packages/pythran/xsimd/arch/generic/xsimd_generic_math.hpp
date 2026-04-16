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

#ifndef XSIMD_GENERIC_MATH_HPP
#define XSIMD_GENERIC_MATH_HPP

#include "../xsimd_scalar.hpp"
#include "./xsimd_generic_details.hpp"
#include "./xsimd_generic_trigo.hpp"

#include <type_traits>

namespace xsimd
{

    namespace kernel
    {

        using namespace types;
        // abs
        template <class A, class T, class>
        XSIMD_INLINE batch<T, A> abs(batch<T, A> const& self, requires_arch<generic>) noexcept
        {
            if (std::is_unsigned<T>::value)
                return self;
            else
            {
                auto sign = bitofsign(self);
                auto inv = self ^ sign;
                return inv - sign;
            }
        }

        template <class A, class T>
        XSIMD_INLINE batch<T, A> abs(batch<std::complex<T>, A> const& z, requires_arch<generic>) noexcept
        {
            return hypot(z.real(), z.imag());
        }

        // avg
        namespace detail
        {
            template <class A, class T>
            XSIMD_INLINE batch<T, A> avg(batch<T, A> const& x, batch<T, A> const& y, std::true_type, std::false_type) noexcept
            {
                return (x & y) + ((x ^ y) >> 1);
            }

            template <class A, class T>
            XSIMD_INLINE batch<T, A> avg(batch<T, A> const& x, batch<T, A> const& y, std::true_type, std::true_type) noexcept
            {
                // Inspired by
                // https://stackoverflow.com/questions/5697500/take-the-average-of-two-signed-numbers-in-c
                auto t = (x & y) + ((x ^ y) >> 1);
                auto t_u = bitwise_cast<typename std::make_unsigned<T>::type>(t);
                auto avg = t + (bitwise_cast<T>(t_u >> (8 * sizeof(T) - 1)) & (x ^ y));
                return avg;
            }

            template <class A, class T>
            XSIMD_INLINE batch<T, A> avg(batch<T, A> const& x, batch<T, A> const& y, std::false_type, std::true_type) noexcept
            {
                return (x + y) / 2;
            }
        }

        template <class A, class T>
        XSIMD_INLINE batch<T, A> avg(batch<T, A> const& x, batch<T, A> const& y, requires_arch<generic>) noexcept
        {
            return detail::avg(x, y, typename std::is_integral<T>::type {}, typename std::is_signed<T>::type {});
        }

        // avgr
        namespace detail
        {
            template <class A, class T>
            XSIMD_INLINE batch<T, A> avgr(batch<T, A> const& x, batch<T, A> const& y, std::true_type) noexcept
            {
                constexpr unsigned shift = 8 * sizeof(T) - 1;
                auto adj = std::is_signed<T>::value ? ((x ^ y) & 0x1) : (((x ^ y) << shift) >> shift);
                return ::xsimd::kernel::avg(x, y, A {}) + adj;
            }

            template <class A, class T>
            XSIMD_INLINE batch<T, A> avgr(batch<T, A> const& x, batch<T, A> const& y, std::false_type) noexcept
            {
                return ::xsimd::kernel::avg(x, y, A {});
            }
        }

        template <class A, class T>
        XSIMD_INLINE batch<T, A> avgr(batch<T, A> const& x, batch<T, A> const& y, requires_arch<generic>) noexcept
        {
            return detail::avgr(x, y, typename std::is_integral<T>::type {});
        }

        // batch_cast
        template <class A, class T>
        XSIMD_INLINE batch<T, A> batch_cast(batch<T, A> const& self, batch<T, A> const&, requires_arch<generic>) noexcept
        {
            return self;
        }

        namespace detail
        {
            template <class A, class T_out, class T_in>
            XSIMD_INLINE batch<T_out, A> batch_cast(batch<T_in, A> const& self, batch<T_out, A> const& out, requires_arch<generic>, with_fast_conversion) noexcept
            {
                return fast_cast(self, out, A {});
            }
            template <class A, class T_out, class T_in>
            XSIMD_INLINE batch<T_out, A> batch_cast(batch<T_in, A> const& self, batch<T_out, A> const&, requires_arch<generic>, with_slow_conversion) noexcept
            {
                static_assert(!std::is_same<T_in, T_out>::value, "there should be no conversion for this type combination");
                using batch_type_in = batch<T_in, A>;
                using batch_type_out = batch<T_out, A>;
                static_assert(batch_type_in::size == batch_type_out::size, "compatible sizes");
                alignas(A::alignment()) T_in buffer_in[batch_type_in::size];
                alignas(A::alignment()) T_out buffer_out[batch_type_out::size];
                self.store_aligned(&buffer_in[0]);
                std::copy(std::begin(buffer_in), std::end(buffer_in), std::begin(buffer_out));
                return batch_type_out::load_aligned(buffer_out);
            }

        }

        template <class A, class T_out, class T_in>
        XSIMD_INLINE batch<T_out, A> batch_cast(batch<T_in, A> const& self, batch<T_out, A> const& out, requires_arch<generic>) noexcept
        {
            return detail::batch_cast(self, out, A {}, detail::conversion_type<A, T_in, T_out> {});
        }

        // bitofsign
        template <class A, class T>
        XSIMD_INLINE batch<T, A> bitofsign(batch<T, A> const& self, requires_arch<generic>) noexcept
        {
            static_assert(std::is_integral<T>::value, "int type implementation");
            if (std::is_unsigned<T>::value)
                return batch<T, A>(0);
            else
                return self >> (T)(8 * sizeof(T) - 1);
        }

        template <class A>
        XSIMD_INLINE batch<float, A> bitofsign(batch<float, A> const& self, requires_arch<generic>) noexcept
        {
            return self & constants::signmask<batch<float, A>>();
        }
        template <class A>
        XSIMD_INLINE batch<double, A> bitofsign(batch<double, A> const& self, requires_arch<generic>) noexcept
        {
            return self & constants::signmask<batch<double, A>>();
        }

        // bitwise_cast
        template <class A, class T>
        XSIMD_INLINE batch<T, A> bitwise_cast(batch<T, A> const& self, batch<T, A> const&, requires_arch<generic>) noexcept
        {
            return self;
        }

        // cbrt
        /* origin: boost/simd/arch/common/simd/function/cbrt.hpp */
        /*
         * ====================================================
         * copyright 2016 NumScale SAS
         *
         * Distributed under the Boost Software License, Version 1.0.
         * (See copy at http://boost.org/LICENSE_1_0.txt)
         * ====================================================
         */
        template <class A>
        XSIMD_INLINE batch<float, A> cbrt(batch<float, A> const& self, requires_arch<generic>) noexcept
        {
            using batch_type = batch<float, A>;
            batch_type z = abs(self);
#ifndef XSIMD_NO_DENORMALS
            auto denormal = z < constants::smallestposval<batch_type>();
            z = select(denormal, z * constants::twotonmb<batch_type>(), z);
            batch_type f = select(denormal, constants::twotonmbo3<batch_type>(), batch_type(1.));
#endif
            const batch_type CBRT2(bit_cast<float>(0x3fa14518));
            const batch_type CBRT4(bit_cast<float>(0x3fcb2ff5));
            const batch_type CBRT2I(bit_cast<float>(0x3f4b2ff5));
            const batch_type CBRT4I(bit_cast<float>(0x3f214518));
            using i_type = as_integer_t<batch_type>;
            i_type e;
            batch_type x = frexp(z, e);
            x = detail::horner<batch_type,
                               0x3ece0609,
                               0x3f91eb77,
                               0xbf745265,
                               0x3f0bf0fe,
                               0xbe09e49a>(x);
            auto flag = e >= i_type(0);
            i_type e1 = abs(e);
            i_type rem = e1;
            e1 /= i_type(3);
            rem -= e1 * i_type(3);
            e = e1 * sign(e);
            const batch_type cbrt2 = select(batch_bool_cast<float>(flag), CBRT2, CBRT2I);
            const batch_type cbrt4 = select(batch_bool_cast<float>(flag), CBRT4, CBRT4I);
            batch_type fact = select(batch_bool_cast<float>(rem == i_type(1)), cbrt2, batch_type(1.));
            fact = select(batch_bool_cast<float>(rem == i_type(2)), cbrt4, fact);
            x = ldexp(x * fact, e);
            x -= (x - z / (x * x)) * batch_type(1.f / 3.f);
#ifndef XSIMD_NO_DENORMALS
            x = (x | bitofsign(self)) * f;
#else
            x = x | bitofsign(self);
#endif
#ifndef XSIMD_NO_INFINITIES
            return select(self == batch_type(0.) || isinf(self), self, x);
#else
            return select(self == batch_type(0.), self, x);
#endif
        }

        template <class A>
        XSIMD_INLINE batch<double, A> cbrt(batch<double, A> const& self, requires_arch<generic>) noexcept
        {
            using batch_type = batch<double, A>;
            batch_type z = abs(self);
#ifndef XSIMD_NO_DENORMALS
            auto denormal = z < constants::smallestposval<batch_type>();
            z = select(denormal, z * constants::twotonmb<batch_type>(), z);
            batch_type f = select(denormal, constants::twotonmbo3<batch_type>(), batch_type(1.));
#endif
            const batch_type CBRT2(bit_cast<double>(int64_t(0x3ff428a2f98d728b)));
            const batch_type CBRT4(bit_cast<double>(int64_t(0x3ff965fea53d6e3d)));
            const batch_type CBRT2I(bit_cast<double>(int64_t(0x3fe965fea53d6e3d)));
            const batch_type CBRT4I(bit_cast<double>(int64_t(0x3fe428a2f98d728b)));
            using i_type = as_integer_t<batch_type>;
            i_type e;
            batch_type x = frexp(z, e);
            x = detail::horner<batch_type,
                               0x3fd9c0c12122a4feull,
                               0x3ff23d6ee505873aull,
                               0xbfee8a4ca3ba37b8ull,
                               0x3fe17e1fc7e59d58ull,
                               0xbfc13c93386fdff6ull>(x);
            auto flag = e >= typename i_type::value_type(0);
            i_type e1 = abs(e);
            i_type rem = e1;
            e1 /= i_type(3);
            rem -= e1 * i_type(3);
            e = e1 * sign(e);
            const batch_type cbrt2 = select(batch_bool_cast<double>(flag), CBRT2, CBRT2I);
            const batch_type cbrt4 = select(batch_bool_cast<double>(flag), CBRT4, CBRT4I);
            batch_type fact = select(batch_bool_cast<double>(rem == i_type(1)), cbrt2, batch_type(1.));
            fact = select(batch_bool_cast<double>(rem == i_type(2)), cbrt4, fact);
            x = ldexp(x * fact, e);
            x -= (x - z / (x * x)) * batch_type(1. / 3.);
            x -= (x - z / (x * x)) * batch_type(1. / 3.);
#ifndef XSIMD_NO_DENORMALS
            x = (x | bitofsign(self)) * f;
#else
            x = x | bitofsign(self);
#endif
#ifndef XSIMD_NO_INFINITIES
            return select(self == batch_type(0.) || isinf(self), self, x);
#else
            return select(self == batch_type(0.), self, x);
#endif
        }

        // clip
        template <class A, class T>
        XSIMD_INLINE batch<T, A> clip(batch<T, A> const& self, batch<T, A> const& lo, batch<T, A> const& hi, requires_arch<generic>) noexcept
        {
            return min(hi, max(self, lo));
        }

        // copysign
        template <class A, class T, class _ = typename std::enable_if<std::is_floating_point<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> copysign(batch<T, A> const& self, batch<T, A> const& other, requires_arch<generic>) noexcept
        {
            return abs(self) | bitofsign(other);
        }

        // erf

        namespace detail
        {
            /* origin: boost/simd/arch/common/detail/generic/erf_kernel.hpp */
            /*
             * ====================================================
             * copyright 2016 NumScale SAS
             *
             * Distributed under the Boost Software License, Version 1.0.
             * (See copy at http://boost.org/LICENSE_1_0.txt)
             * ====================================================
             */
            template <class B>
            struct erf_kernel;

            template <class A>
            struct erf_kernel<batch<float, A>>
            {
                using batch_type = batch<float, A>;
                // computes erf(a0)/a0
                // x is sqr(a0) and 0 <= abs(a0) <= 2/3
                static XSIMD_INLINE batch_type erf1(const batch_type& x) noexcept
                {
                    return detail::horner<batch_type,
                                          0x3f906eba, //   1.128379154774254e+00
                                          0xbec0937e, //  -3.761252839094832e-01
                                          0x3de70f22, //   1.128218315189123e-01
                                          0xbcdb61f4, //  -2.678010670585737e-02
                                          0x3ba4468d, //   5.013293006147870e-03
                                          0xba1fc83b //  -6.095205117313012e-04
                                          >(x);
                }

                // computes erfc(x)*exp(sqr(x))
                // x >=  2/3
                static XSIMD_INLINE batch_type erfc2(const batch_type& x) noexcept
                {
                    return detail::horner<batch_type,
                                          0x3f0a0e8b, //   5.392844046572836e-01
                                          0xbf918a62, //  -1.137035586823118e+00
                                          0x3e243828, //   1.603704761054187e-01
                                          0x3ec4ca6e, //   3.843569094305250e-01
                                          0x3e1175c7, //   1.420508523645926e-01
                                          0x3e2006f0, //   1.562764709849380e-01
                                          0xbfaea865, //  -1.364514006347145e+00
                                          0x4050b063, //   3.260765682222576e+00
                                          0xc0cd1a85, //  -6.409487379234005e+00
                                          0x40d67e3b, //   6.702908785399893e+00
                                          0xc0283611 //  -2.628299919293280e+00
                                          >(x);
                }

                static XSIMD_INLINE batch_type erfc3(const batch_type& x) noexcept
                {
                    return (batch_type(1.) - x) * detail::horner<batch_type,
                                                                 0x3f7ffffe, //   9.9999988e-01
                                                                 0xbe036d7e, //  -1.2834737e-01
                                                                 0xbfa11698, //  -1.2585020e+00
                                                                 0xbffc9284, //  -1.9732213e+00
                                                                 0xc016c985, //  -2.3560498e+00
                                                                 0x3f2cff3b, //   6.7576951e-01
                                                                 0xc010d956, //  -2.2632651e+00
                                                                 0x401b5680, //   2.4271545e+00
                                                                 0x41aa8e55 //   2.1319498e+01
                                                                 >(x);
                }
            };

            template <class A>
            struct erf_kernel<batch<double, A>>
            {
                using batch_type = batch<double, A>;
                // computes erf(a0)/a0
                // x is sqr(a0) and 0 <= abs(a0) <= 0.65
                static XSIMD_INLINE batch_type erf1(const batch_type& x) noexcept
                {
                    return detail::horner<batch_type,
                                          0x3ff20dd750429b61ull, // 1.12837916709551
                                          0x3fc16500f106c0a5ull, // 0.135894887627278
                                          0x3fa4a59a4f02579cull, // 4.03259488531795E-02
                                          0x3f53b7664358865aull, // 1.20339380863079E-03
                                          0x3f110512d5b20332ull // 6.49254556481904E-05
                                          >(x)
                        / detail::horner<batch_type,
                                         0x3ff0000000000000ull, // 1
                                         0x3fdd0a84eb1ca867ull, // 0.453767041780003
                                         0x3fb64536ca92ea2full, // 8.69936222615386E-02
                                         0x3f8166f75999dbd1ull, // 8.49717371168693E-03
                                         0x3f37ea4332348252ull // 3.64915280629351E-04
                                         >(x);
                }

                // computes erfc(x)*exp(x*x)
                // 0.65 <= abs(x) <= 2.2
                static XSIMD_INLINE batch_type erfc2(const batch_type& x) noexcept
                {
                    return detail::horner<batch_type,
                                          0x3feffffffbbb552bull, // 0.999999992049799
                                          0x3ff54dfe9b258a60ull, // 1.33154163936765
                                          0x3fec1986509e687bull, // 0.878115804155882
                                          0x3fd53dd7a67c7e9full, // 0.331899559578213
                                          0x3fb2488a6b5cb5e5ull, // 7.14193832506776E-02
                                          0x3f7cf4cfe0aacbb4ull, // 7.06940843763253E-03
                                          0x0ull // 0
                                          >(x)
                        / detail::horner<batch_type,
                                         0x3ff0000000000000ull, // 1
                                         0x4003adeae79b9708ull, // 2.45992070144246
                                         0x40053b1052dca8bdull, // 2.65383972869776
                                         0x3ff9e677c2777c3cull, // 1.61876655543871
                                         0x3fe307622fcff772ull, // 0.594651311286482
                                         0x3fc033c113a7deeeull, // 0.126579413030178
                                         0x3f89a996639b0d00ull // 1.25304936549413E-02
                                         >(x);
                }

                // computes erfc(x)*exp(x*x)
                // 2.2 <= abs(x) <= 6
                static XSIMD_INLINE batch_type erfc3(const batch_type& x) noexcept
                {
                    return detail::horner<batch_type,
                                          0x3fefff5a9e697ae2ull, // 0.99992114009714
                                          0x3ff9fa202deb88e5ull, // 1.62356584489367
                                          0x3ff44744306832aeull, // 1.26739901455873
                                          0x3fe29be1cff90d94ull, // 0.581528574177741
                                          0x3fc42210f88b9d43ull, // 0.157289620742839
                                          0x3f971d0907ea7a92ull, // 2.25716982919218E-02
                                          0x0ll // 0
                                          >(x)
                        / detail::horner<batch_type,
                                         0x3ff0000000000000ull, // 1
                                         0x400602f24bf3fdb6ull, // 2.75143870676376
                                         0x400afd487397568full, // 3.37367334657285
                                         0x400315ffdfd5ce91ull, // 2.38574194785344
                                         0x3ff0cfd4cb6cde9full, // 1.05074004614827
                                         0x3fd1d7ab774bb837ull, // 0.278788439273629
                                         0x3fa47bd61bbb3843ull // 4.00072964526861E-02
                                         >(x);
                }

                // computes erfc(rx)*exp(rx*rx)
                // x >=  6 rx = 1/x
                static XSIMD_INLINE batch_type erfc4(const batch_type& x) noexcept
                {
                    return detail::horner<batch_type,
                                          0xbc7e4ad1ec7d0000ll, // -2.627435221016534e-17
                                          0x3fe20dd750429a16ll, // 5.641895835477182e-01
                                          0x3db60000e984b501ll, // 2.000889609806154e-11
                                          0xbfd20dd753ae5dfdll, // -2.820947949598745e-01
                                          0x3e907e71e046a820ll, // 2.457786367990903e-07
                                          0x3fdb1494cac06d39ll, // 4.231311779019112e-01
                                          0x3f34a451701654f1ll, // 3.149699042180451e-04
                                          0xbff105e6b8ef1a63ll, // -1.063940737150596e+00
                                          0x3fb505a857e9ccc8ll, // 8.211757799454056e-02
                                          0x40074fbabc514212ll, // 2.913930388669777e+00
                                          0x4015ac7631f7ac4fll, // 5.418419628850713e+00
                                          0xc0457e03041e9d8bll, // -4.298446704382794e+01
                                          0x4055803d26c4ec4fll, // 8.600373238783617e+01
                                          0xc0505fce04ec4ec5ll // -6.549694941594051e+01
                                          >(x);
                }
            };
        }
        /* origin: boost/simd/arch/common/simd/function/erf.hpp */
        /*
         * ====================================================
         * copyright 2016 NumScale SAS
         *
         * Distributed under the Boost Software License, Version 1.0.
         * (See copy at http://boost.org/LICENSE_1_0.txt)
         * ====================================================
         */

        template <class A>
        XSIMD_INLINE batch<float, A> erf(batch<float, A> const& self, requires_arch<generic>) noexcept
        {
            using batch_type = batch<float, A>;
            batch_type x = abs(self);
            batch_type r1(0.);
            auto test1 = x < batch_type(2.f / 3.f);
            if (any(test1))
            {
                r1 = self * detail::erf_kernel<batch_type>::erf1(x * x);
                if (all(test1))
                    return r1;
            }
            batch_type z = x / (batch_type(1.) + x);
            z -= batch_type(0.4f);
            batch_type r2 = batch_type(1.) - exp(-x * x) * detail::erf_kernel<batch_type>::erfc2(z);
            r2 = select(self < batch_type(0.), -r2, r2);
            r1 = select(test1, r1, r2);
#ifndef XSIMD_NO_INFINITIES
            r1 = select(xsimd::isinf(self), sign(self), r1);
#endif
            return r1;
        }

        template <class A>
        XSIMD_INLINE batch<double, A> erf(batch<double, A> const& self, requires_arch<generic>) noexcept
        {
            using batch_type = batch<double, A>;
            batch_type x = abs(self);
            batch_type xx = x * x;
            batch_type lim1(0.65);
            batch_type lim2(2.2);
            auto test1 = x < lim1;
            batch_type r1(0.);
            if (any(test1))
            {
                r1 = self * detail::erf_kernel<batch_type>::erf1(xx);
                if (all(test1))
                    return r1;
            }
            auto test2 = x < lim2;
            auto test3 = test2 && !test1;
            batch_type ex = exp(-xx);
            if (any(test3))
            {
                batch_type z = batch_type(1.) - ex * detail::erf_kernel<batch_type>::erfc2(x);
                batch_type r2 = select(self < batch_type(0.), -z, z);
                r1 = select(test1, r1, r2);
                if (all(test1 || test3))
                    return r1;
            }
            batch_type z = batch_type(1.) - ex * detail::erf_kernel<batch_type>::erfc3(x);
            z = select(self < batch_type(0.), -z, z);
#ifndef XSIMD_NO_INFINITIES
            z = select(xsimd::isinf(self), sign(self), z);
#endif
            return select(test2, r1, z);
        }

        // erfc
        template <class A>
        XSIMD_INLINE batch<float, A> erfc(batch<float, A> const& self, requires_arch<generic>) noexcept
        {
            using batch_type = batch<float, A>;
            batch_type x = abs(self);
            auto test0 = self < batch_type(0.);
            batch_type r1(0.);
            auto test1 = 3.f * x < 2.f;
            batch_type z = x / (batch_type(1.) + x);
            if (any(test1))
            {
                r1 = detail::erf_kernel<batch_type>::erfc3(z);
                if (all(test1))
                    return select(test0, batch_type(2.) - r1, r1);
            }

            z -= batch_type(0.4f);
            batch_type r2 = exp(-x * x) * detail::erf_kernel<batch_type>::erfc2(z);
            r1 = select(test1, r1, r2);
#ifndef XSIMD_NO_INFINITIES
            r1 = select(x == constants::infinity<batch_type>(), batch_type(0.), r1);
#endif
            return select(test0, batch_type(2.) - r1, r1);
        }

        template <class A>
        XSIMD_INLINE batch<double, A> erfc(batch<double, A> const& self, requires_arch<generic>) noexcept
        {
            using batch_type = batch<double, A>;
            batch_type x = abs(self);
            batch_type xx = x * x;
            batch_type lim1(0.65);
            batch_type lim2(2.2);
            auto test0 = self < batch_type(0.);
            auto test1 = x < lim1;
            batch_type r1(0.);
            if (any(test1))
            {
                r1 = batch_type(1.) - x * detail::erf_kernel<batch_type>::erf1(xx);
                if (all(test1))
                    return select(test0, batch_type(2.) - r1, r1);
            }
            auto test2 = x < lim2;
            auto test3 = test2 && !test1;
            batch_type ex = exp(-xx);
            if (any(test3))
            {
                batch_type z = ex * detail::erf_kernel<batch_type>::erfc2(x);
                r1 = select(test1, r1, z);
                if (all(test1 || test3))
                    return select(test0, batch_type(2.) - r1, r1);
            }
            batch_type z = ex * detail::erf_kernel<batch_type>::erfc3(x);
            r1 = select(test2, r1, z);
#ifndef XSIMD_NO_INFINITIES
            r1 = select(x == constants::infinity<batch_type>(), batch_type(0.), r1);
#endif
            return select(test0, batch_type(2.) - r1, r1);
        }

        // estrin
        namespace detail
        {

            template <class B>
            struct estrin
            {
                B x;

                template <typename... Ts>
                XSIMD_INLINE B operator()(const Ts&... coefs) noexcept
                {
                    return eval(coefs...);
                }

            private:
                XSIMD_INLINE B eval(const B& c0) noexcept
                {
                    return c0;
                }

                XSIMD_INLINE B eval(const B& c0, const B& c1) noexcept
                {
                    return fma(x, c1, c0);
                }

                template <size_t... Is, class Tuple>
                XSIMD_INLINE B eval(::xsimd::detail::index_sequence<Is...>, const Tuple& tuple)
                {
                    return estrin { x * x }(std::get<Is>(tuple)...);
                }

                template <class... Args>
                XSIMD_INLINE B eval(const std::tuple<Args...>& tuple) noexcept
                {
                    return eval(::xsimd::detail::make_index_sequence<sizeof...(Args)>(), tuple);
                }

                template <class... Args>
                XSIMD_INLINE B eval(const std::tuple<Args...>& tuple, const B& c0) noexcept
                {
                    return eval(std::tuple_cat(tuple, std::make_tuple(eval(c0))));
                }

                template <class... Args>
                XSIMD_INLINE B eval(const std::tuple<Args...>& tuple, const B& c0, const B& c1) noexcept
                {
                    return eval(std::tuple_cat(tuple, std::make_tuple(eval(c0, c1))));
                }

                template <class... Args, class... Ts>
                XSIMD_INLINE B eval(const std::tuple<Args...>& tuple, const B& c0, const B& c1, const Ts&... coefs) noexcept
                {
                    return eval(std::tuple_cat(tuple, std::make_tuple(eval(c0, c1))), coefs...);
                }

                template <class... Ts>
                XSIMD_INLINE B eval(const B& c0, const B& c1, const Ts&... coefs) noexcept
                {
                    return eval(std::make_tuple(eval(c0, c1)), coefs...);
                }
            };
        }

        template <class T, class A, uint64_t... Coefs>
        XSIMD_INLINE batch<T, A> estrin(const batch<T, A>& self) noexcept
        {
            using batch_type = batch<T, A>;
            return detail::estrin<batch_type> { self }(detail::coef<batch_type, Coefs>()...);
        }

        // exp
        /* origin: boost/simd/arch/common/detail/simd/expo_base.hpp */
        /*
         * ====================================================
         * copyright 2016 NumScale SAS
         *
         * Distributed under the Boost Software License, Version 1.0.
         * (See copy at http://boost.org/LICENSE_1_0.txt)
         * ====================================================
         */
        namespace detail
        {
            enum exp_reduction_tag
            {
                exp_tag,
                exp2_tag,
                exp10_tag
            };

            template <class B, exp_reduction_tag Tag>
            struct exp_reduction_base;

            template <class B>
            struct exp_reduction_base<B, exp_tag>
            {
                static constexpr B maxlog() noexcept
                {
                    return constants::maxlog<B>();
                }

                static constexpr B minlog() noexcept
                {
                    return constants::minlog<B>();
                }
            };

            template <class B>
            struct exp_reduction_base<B, exp10_tag>
            {
                static constexpr B maxlog() noexcept
                {
                    return constants::maxlog10<B>();
                }

                static constexpr B minlog() noexcept
                {
                    return constants::minlog10<B>();
                }
            };

            template <class B>
            struct exp_reduction_base<B, exp2_tag>
            {
                static constexpr B maxlog() noexcept
                {
                    return constants::maxlog2<B>();
                }

                static constexpr B minlog() noexcept
                {
                    return constants::minlog2<B>();
                }
            };

            template <class T, class A, exp_reduction_tag Tag>
            struct exp_reduction;

            template <class A>
            struct exp_reduction<float, A, exp_tag> : exp_reduction_base<batch<float, A>, exp_tag>
            {
                using batch_type = batch<float, A>;
                static XSIMD_INLINE batch_type approx(const batch_type& x) noexcept
                {
                    batch_type y = detail::horner<batch_type,
                                                  0x3f000000, //  5.0000000e-01
                                                  0x3e2aa9a5, //  1.6666277e-01
                                                  0x3d2aa957, //  4.1665401e-02
                                                  0x3c098d8b, //  8.3955629e-03
                                                  0x3ab778cf //  1.3997796e-03
                                                  >(x);
                    return ++fma(y, x * x, x);
                }

                static XSIMD_INLINE batch_type reduce(const batch_type& a, batch_type& x) noexcept
                {
                    batch_type k = nearbyint(constants::invlog_2<batch_type>() * a);
                    x = fnma(k, constants::log_2hi<batch_type>(), a);
                    x = fnma(k, constants::log_2lo<batch_type>(), x);
                    return k;
                }
            };

            template <class A>
            struct exp_reduction<float, A, exp10_tag> : exp_reduction_base<batch<float, A>, exp10_tag>
            {
                using batch_type = batch<float, A>;
                static XSIMD_INLINE batch_type approx(const batch_type& x) noexcept
                {
                    return ++(detail::horner<batch_type,
                                             0x40135d8e, //    2.3025851e+00
                                             0x4029a926, //    2.6509490e+00
                                             0x400237da, //    2.0346589e+00
                                             0x3f95eb4c, //    1.1712432e+00
                                             0x3f0aacef, //    5.4170126e-01
                                             0x3e54dff1 //    2.0788552e-01
                                             >(x)
                              * x);
                }

                static XSIMD_INLINE batch_type reduce(const batch_type& a, batch_type& x) noexcept
                {
                    batch_type k = nearbyint(constants::invlog10_2<batch_type>() * a);
                    x = fnma(k, constants::log10_2hi<batch_type>(), a);
                    x -= k * constants::log10_2lo<batch_type>();
                    return k;
                }
            };

            template <class A>
            struct exp_reduction<float, A, exp2_tag> : exp_reduction_base<batch<float, A>, exp2_tag>
            {
                using batch_type = batch<float, A>;
                static XSIMD_INLINE batch_type approx(const batch_type& x) noexcept
                {
                    batch_type y = detail::horner<batch_type,
                                                  0x3e75fdf1, //    2.4022652e-01
                                                  0x3d6356eb, //    5.5502813e-02
                                                  0x3c1d9422, //    9.6178371e-03
                                                  0x3ab01218, //    1.3433127e-03
                                                  0x3922c8c4 //    1.5524315e-04
                                                  >(x);
                    return ++fma(y, x * x, x * constants::log_2<batch_type>());
                }

                static XSIMD_INLINE batch_type reduce(const batch_type& a, batch_type& x) noexcept
                {
                    batch_type k = nearbyint(a);
                    x = (a - k);
                    return k;
                }
            };

            template <class A>
            struct exp_reduction<double, A, exp_tag> : exp_reduction_base<batch<double, A>, exp_tag>
            {
                using batch_type = batch<double, A>;
                static XSIMD_INLINE batch_type approx(const batch_type& x) noexcept
                {
                    batch_type t = x * x;
                    return fnma(t,
                                detail::horner<batch_type,
                                               0x3fc555555555553eull,
                                               0xbf66c16c16bebd93ull,
                                               0x3f11566aaf25de2cull,
                                               0xbebbbd41c5d26bf1ull,
                                               0x3e66376972bea4d0ull>(t),
                                x);
                }

                static XSIMD_INLINE batch_type reduce(const batch_type& a, batch_type& hi, batch_type& lo, batch_type& x) noexcept
                {
                    batch_type k = nearbyint(constants::invlog_2<batch_type>() * a);
                    hi = fnma(k, constants::log_2hi<batch_type>(), a);
                    lo = k * constants::log_2lo<batch_type>();
                    x = hi - lo;
                    return k;
                }

                static XSIMD_INLINE batch_type finalize(const batch_type& x, const batch_type& c, const batch_type& hi, const batch_type& lo) noexcept
                {
                    return batch_type(1.) - (((lo - (x * c) / (batch_type(2.) - c)) - hi));
                }
            };

            template <class A>
            struct exp_reduction<double, A, exp10_tag> : exp_reduction_base<batch<double, A>, exp10_tag>
            {
                using batch_type = batch<double, A>;
                static XSIMD_INLINE batch_type approx(const batch_type& x) noexcept
                {
                    batch_type xx = x * x;
                    batch_type px = x * detail::horner<batch_type, 0x40a2b4798e134a01ull, 0x40796b7a050349e4ull, 0x40277d9474c55934ull, 0x3fa4fd75f3062dd4ull>(xx);
                    batch_type x2 = px / (detail::horner1<batch_type, 0x40a03f37650df6e2ull, 0x4093e05eefd67782ull, 0x405545fdce51ca08ull>(xx) - px);
                    return ++(x2 + x2);
                }

                static XSIMD_INLINE batch_type reduce(const batch_type& a, batch_type&, batch_type&, batch_type& x) noexcept
                {
                    batch_type k = nearbyint(constants::invlog10_2<batch_type>() * a);
                    x = fnma(k, constants::log10_2hi<batch_type>(), a);
                    x = fnma(k, constants::log10_2lo<batch_type>(), x);
                    return k;
                }

                static XSIMD_INLINE batch_type finalize(const batch_type&, const batch_type& c, const batch_type&, const batch_type&) noexcept
                {
                    return c;
                }
            };

            template <class A>
            struct exp_reduction<double, A, exp2_tag> : exp_reduction_base<batch<double, A>, exp2_tag>
            {
                using batch_type = batch<double, A>;
                static XSIMD_INLINE batch_type approx(const batch_type& x) noexcept
                {
                    batch_type t = x * x;
                    return fnma(t,
                                detail::horner<batch_type,
                                               0x3fc555555555553eull,
                                               0xbf66c16c16bebd93ull,
                                               0x3f11566aaf25de2cull,
                                               0xbebbbd41c5d26bf1ull,
                                               0x3e66376972bea4d0ull>(t),
                                x);
                }

                static XSIMD_INLINE batch_type reduce(const batch_type& a, batch_type&, batch_type&, batch_type& x) noexcept
                {
                    batch_type k = nearbyint(a);
                    x = (a - k) * constants::log_2<batch_type>();
                    return k;
                }

                static XSIMD_INLINE batch_type finalize(const batch_type& x, const batch_type& c, const batch_type&, const batch_type&) noexcept
                {
                    return batch_type(1.) + x + x * c / (batch_type(2.) - c);
                }
            };

            template <exp_reduction_tag Tag, class A>
            XSIMD_INLINE batch<float, A> exp(batch<float, A> const& self) noexcept
            {
                using batch_type = batch<float, A>;
                using reducer_t = exp_reduction<float, A, Tag>;
                batch_type x;
                batch_type k = reducer_t::reduce(self, x);
                x = reducer_t::approx(x);
                x = select(self <= reducer_t::minlog(), batch_type(0.), ldexp(x, to_int(k)));
                x = select(self >= reducer_t::maxlog(), constants::infinity<batch_type>(), x);
                return x;
            }

            template <exp_reduction_tag Tag, class A>
            XSIMD_INLINE batch<double, A> exp(batch<double, A> const& self) noexcept
            {
                using batch_type = batch<double, A>;
                using reducer_t = exp_reduction<double, A, Tag>;
                batch_type hi, lo, x;
                batch_type k = reducer_t::reduce(self, hi, lo, x);
                batch_type c = reducer_t::approx(x);
                c = reducer_t::finalize(x, c, hi, lo);
                c = select(self <= reducer_t::minlog(), batch_type(0.), ldexp(c, to_int(k)));
                c = select(self >= reducer_t::maxlog(), constants::infinity<batch_type>(), c);
                return c;
            }
        }

        template <class A, class T>
        XSIMD_INLINE batch<T, A> exp(batch<T, A> const& self, requires_arch<generic>) noexcept
        {
            return detail::exp<detail::exp_tag>(self);
        }

        template <class A, class T>
        XSIMD_INLINE batch<std::complex<T>, A> exp(batch<std::complex<T>, A> const& self, requires_arch<generic>) noexcept
        {
            using batch_type = batch<std::complex<T>, A>;
            auto isincos = sincos(self.imag());
            return exp(self.real()) * batch_type(std::get<1>(isincos), std::get<0>(isincos));
        }

        // exp10
        template <class A, class T>
        XSIMD_INLINE batch<T, A> exp10(batch<T, A> const& self, requires_arch<generic>) noexcept
        {
            return detail::exp<detail::exp10_tag>(self);
        }

        // exp2
        template <class A, class T>
        XSIMD_INLINE batch<T, A> exp2(batch<T, A> const& self, requires_arch<generic>) noexcept
        {
            return detail::exp<detail::exp2_tag>(self);
        }

        // expm1
        namespace detail
        {
            /* origin: boost/simd/arch/common/detail/generic/expm1_kernel.hpp */
            /*
             * ====================================================
             * copyright 2016 NumScale SAS
             *
             * Distributed under the Boost Software License, Version 1.0.
             * (See copy at http://boost.org/LICENSE_1_0.txt)
             * ====================================================
             */
            template <class A>
            static XSIMD_INLINE batch<float, A> expm1(const batch<float, A>& a) noexcept
            {
                using batch_type = batch<float, A>;
                batch_type k = nearbyint(constants::invlog_2<batch_type>() * a);
                batch_type x = fnma(k, constants::log_2hi<batch_type>(), a);
                x = fnma(k, constants::log_2lo<batch_type>(), x);
                batch_type hx = x * batch_type(0.5);
                batch_type hxs = x * hx;
                batch_type r = detail::horner<batch_type,
                                              0X3F800000UL, // 1
                                              0XBD08887FUL, // -3.3333298E-02
                                              0X3ACF6DB4UL // 1.582554
                                              >(hxs);
                batch_type t = fnma(r, hx, batch_type(3.));
                batch_type e = hxs * ((r - t) / (batch_type(6.) - x * t));
                e = fms(x, e, hxs);
                using i_type = as_integer_t<batch_type>;
                i_type ik = to_int(k);
                batch_type two2mk = ::xsimd::bitwise_cast<float>((constants::maxexponent<batch_type>() - ik) << constants::nmb<batch_type>());
                batch_type y = batch_type(1.) - two2mk - (e - x);
                return ldexp(y, ik);
            }

            template <class A>
            static XSIMD_INLINE batch<double, A> expm1(const batch<double, A>& a) noexcept
            {
                using batch_type = batch<double, A>;
                batch_type k = nearbyint(constants::invlog_2<batch_type>() * a);
                batch_type hi = fnma(k, constants::log_2hi<batch_type>(), a);
                batch_type lo = k * constants::log_2lo<batch_type>();
                batch_type x = hi - lo;
                batch_type hxs = x * x * batch_type(0.5);
                batch_type r = detail::horner<batch_type,
                                              0X3FF0000000000000ULL,
                                              0XBFA11111111110F4ULL,
                                              0X3F5A01A019FE5585ULL,
                                              0XBF14CE199EAADBB7ULL,
                                              0X3ED0CFCA86E65239ULL,
                                              0XBE8AFDB76E09C32DULL>(hxs);
                batch_type t = batch_type(3.) - r * batch_type(0.5) * x;
                batch_type e = hxs * ((r - t) / (batch_type(6) - x * t));
                batch_type c = (hi - x) - lo;
                e = (x * (e - c) - c) - hxs;
                using i_type = as_integer_t<batch_type>;
                i_type ik = to_int(k);
                batch_type two2mk = ::xsimd::bitwise_cast<double>((constants::maxexponent<batch_type>() - ik) << constants::nmb<batch_type>());
                batch_type ct1 = batch_type(1.) - two2mk - (e - x);
                batch_type ct2 = ++(x - (e + two2mk));
                batch_type y = select(k < batch_type(20.), ct1, ct2);
                return ldexp(y, ik);
            }

        }

        template <class A, class T>
        XSIMD_INLINE batch<T, A> expm1(batch<T, A> const& self, requires_arch<generic>) noexcept
        {
            using batch_type = batch<T, A>;
            return select(self < constants::logeps<batch_type>(),
                          batch_type(-1.),
                          select(self > constants::maxlog<batch_type>(),
                                 constants::infinity<batch_type>(),
                                 detail::expm1(self)));
        }

        template <class A, class T>
        XSIMD_INLINE batch<std::complex<T>, A> expm1(const batch<std::complex<T>, A>& z, requires_arch<generic>) noexcept
        {
            using batch_type = batch<std::complex<T>, A>;
            using real_batch = typename batch_type::real_batch;
            real_batch isin = sin(z.imag());
            real_batch rem1 = expm1(z.real());
            real_batch re = rem1 + 1.;
            real_batch si = sin(z.imag() * 0.5);
            return { rem1 - 2. * re * si * si, re * isin };
        }

        // polar
        template <class A, class T>
        XSIMD_INLINE batch<std::complex<T>, A> polar(const batch<T, A>& r, const batch<T, A>& theta, requires_arch<generic>) noexcept
        {
            auto sincosTheta = sincos(theta);
            return { r * sincosTheta.second, r * sincosTheta.first };
        }

        // fdim
        template <class A, class T>
        XSIMD_INLINE batch<T, A> fdim(batch<T, A> const& self, batch<T, A> const& other, requires_arch<generic>) noexcept
        {
            return fmax(batch<T, A>(0), self - other);
        }

        // fmod
        template <class A, class T>
        XSIMD_INLINE batch<T, A> fmod(batch<T, A> const& self, batch<T, A> const& other, requires_arch<generic>) noexcept
        {
            return fnma(trunc(self / other), other, self);
        }

        // frexp
        /* origin: boost/simd/arch/common/simd/function/ifrexp.hpp */
        /*
         * ====================================================
         * copyright 2016 NumScale SAS
         *
         * Distributed under the Boost Software License, Version 1.0.
         * (See copy at http://boost.org/LICENSE_1_0.txt)
         * ====================================================
         */
        template <class A, class T>
        XSIMD_INLINE batch<T, A> frexp(const batch<T, A>& self, batch<as_integer_t<T>, A>& exp, requires_arch<generic>) noexcept
        {
            using batch_type = batch<T, A>;
            using int_type = as_integer_t<T>;
            using i_type = batch<int_type, A>;
            i_type m1f = constants::mask1frexp<batch_type>();
            i_type r1 = m1f & ::xsimd::bitwise_cast<int_type>(self);
            batch_type x = self & ::xsimd::bitwise_cast<T>(~m1f);
            exp = (r1 >> constants::nmb<batch_type>()) - constants::maxexponentm1<batch_type>();
            exp = select(batch_bool_cast<typename i_type::value_type>(self != batch_type(0.)), exp, i_type(typename i_type::value_type(0)));
            return select((self != batch_type(0.)), x | ::xsimd::bitwise_cast<T>(constants::mask2frexp<batch_type>()), batch_type(0.));
        }

        // from bool
        template <class A, class T>
        XSIMD_INLINE batch<T, A> from_bool(batch_bool<T, A> const& self, requires_arch<generic>) noexcept
        {
            return batch<T, A>(self.data) & batch<T, A>(1);
        }

        // horner
        template <class T, class A, uint64_t... Coefs>
        XSIMD_INLINE batch<T, A> horner(const batch<T, A>& self) noexcept
        {
            return detail::horner<batch<T, A>, Coefs...>(self);
        }

        // hypot
        template <class A, class T>
        XSIMD_INLINE batch<T, A> hypot(batch<T, A> const& self, batch<T, A> const& other, requires_arch<generic>) noexcept
        {
            return sqrt(fma(self, self, other * other));
        }

        // ipow
        template <class A, class T, class ITy>
        XSIMD_INLINE batch<T, A> ipow(batch<T, A> const& self, ITy other, requires_arch<generic>) noexcept
        {
            return ::xsimd::detail::ipow(self, other);
        }

        // ldexp
        /* origin: boost/simd/arch/common/simd/function/ldexp.hpp */
        /*
         * ====================================================
         * copyright 2016 NumScale SAS
         *
         * Distributed under the Boost Software License, Version 1.0.
         * (See copy at http://boost.org/LICENSE_1_0.txt)
         * ====================================================
         */
        template <class A, class T>
        XSIMD_INLINE batch<T, A> ldexp(const batch<T, A>& self, const batch<as_integer_t<T>, A>& other, requires_arch<generic>) noexcept
        {
            using batch_type = batch<T, A>;
            using itype = as_integer_t<batch_type>;
            itype ik = other + constants::maxexponent<T>();
            ik = ik << constants::nmb<T>();
            return self * ::xsimd::bitwise_cast<T>(ik);
        }

        // lgamma
        template <class A, class T>
        XSIMD_INLINE batch<T, A> lgamma(batch<T, A> const& self, requires_arch<generic>) noexcept;

        namespace detail
        {
            /* origin: boost/simd/arch/common/detail/generic/gammaln_kernel.hpp */
            /*
             * ====================================================
             * copyright 2016 NumScale SAS
             *
             * Distributed under the Boost Software License, Version 1.0.
             * (See copy at http://boost.org/LICENSE_1_0.txt)
             * ====================================================
             */
            template <class A>
            static XSIMD_INLINE batch<float, A> gammalnB(const batch<float, A>& x) noexcept
            {
                return horner<batch<float, A>,
                              0x3ed87730, //    4.227843421859038E-001
                              0x3ea51a64, //    3.224669577325661E-001,
                              0xbd89f07e, //   -6.735323259371034E-002,
                              0x3ca89ed8, //    2.058355474821512E-002,
                              0xbbf164fd, //   -7.366775108654962E-003,
                              0x3b3ba883, //    2.863437556468661E-003,
                              0xbaabeab1, //   -1.311620815545743E-003,
                              0x3a1ebb94 //    6.055172732649237E-004
                              >(x);
            }

            template <class A>
            static XSIMD_INLINE batch<float, A> gammalnC(const batch<float, A>& x) noexcept
            {
                return horner<batch<float, A>,
                              0xbf13c468, //   -5.772156501719101E-001
                              0x3f528d34, //    8.224670749082976E-001,
                              0xbecd27a8, //   -4.006931650563372E-001,
                              0x3e8a898b, //    2.705806208275915E-001,
                              0xbe53c04f, //   -2.067882815621965E-001,
                              0x3e2d4dab, //    1.692415923504637E-001,
                              0xbe22d329, //   -1.590086327657347E-001,
                              0x3e0c3c4f //    1.369488127325832E-001
                              >(x);
            }

            template <class A>
            static XSIMD_INLINE batch<float, A> gammaln2(const batch<float, A>& x) noexcept
            {
                return horner<batch<float, A>,
                              0x3daaaa94, //   8.333316229807355E-002f
                              0xbb358701, //  -2.769887652139868E-003f,
                              0x3a31fd69 //   6.789774945028216E-004f
                              >(x);
            }

            template <class A>
            static XSIMD_INLINE batch<double, A> gammaln1(const batch<double, A>& x) noexcept
            {
                return horner<batch<double, A>,
                              0xc12a0c675418055eull, //  -8.53555664245765465627E5
                              0xc13a45890219f20bull, //  -1.72173700820839662146E6,
                              0xc131bc82f994db51ull, //  -1.16237097492762307383E6,
                              0xc1143d73f89089e5ull, //  -3.31612992738871184744E5,
                              0xc0e2f234355bb93eull, //  -3.88016315134637840924E4,
                              0xc09589018ff36761ull //  -1.37825152569120859100E3
                              >(x)
                    / horner<batch<double, A>,
                             0xc13ece4b6a11e14aull, //  -2.01889141433532773231E6
                             0xc1435255892ff34cull, //  -2.53252307177582951285E6,
                             0xc131628671950043ull, //  -1.13933444367982507207E6,
                             0xc10aeb84b9744c9bull, //  -2.20528590553854454839E5,
                             0xc0d0aa0d7b89d757ull, //  -1.70642106651881159223E4,
                             0xc075fd0d1cf312b2ull, //  -3.51815701436523470549E2,
                             0x3ff0000000000000ull //   1.00000000000000000000E0
                             >(x);
            }

            template <class A>
            static XSIMD_INLINE batch<double, A> gammalnA(const batch<double, A>& x) noexcept
            {
                return horner<batch<double, A>,
                              0x3fb555555555554bull, //    8.33333333333331927722E-2
                              0xbf66c16c16b0a5a1ull, //   -2.77777777730099687205E-3,
                              0x3f4a019f20dc5ebbull, //    7.93650340457716943945E-4,
                              0xbf437fbdb580e943ull, //   -5.95061904284301438324E-4,
                              0x3f4a985027336661ull //    8.11614167470508450300E-4
                              >(x);
            }

            /* origin: boost/simd/arch/common/simd/function/gammaln.hpp */
            /*
             * ====================================================
             * copyright 2016 NumScale SAS
             *
             * Distributed under the Boost Software License, Version 1.0.
             * (See copy at http://boost.org/LICENSE_1_0.txt)
             * ====================================================
             */
            template <class B>
            struct lgamma_impl;

            template <class A>
            struct lgamma_impl<batch<float, A>>
            {
                using batch_type = batch<float, A>;
                static XSIMD_INLINE batch_type compute(const batch_type& a) noexcept
                {
                    auto inf_result = (a <= batch_type(0.)) && is_flint(a);
                    batch_type x = select(inf_result, constants::nan<batch_type>(), a);
                    batch_type q = abs(x);
#ifndef XSIMD_NO_INFINITIES
                    inf_result = (x == constants::infinity<batch_type>()) || inf_result;
#endif
                    auto ltza = a < batch_type(0.);
                    batch_type r;
                    batch_type r1 = other(q);
                    if (any(ltza))
                    {
                        r = select(inf_result, constants::infinity<batch_type>(), negative(q, r1));
                        if (all(ltza))
                            return r;
                    }
                    batch_type r2 = select(ltza, r, r1);
                    return select(a == constants::minusinfinity<batch_type>(), constants::nan<batch_type>(), select(inf_result, constants::infinity<batch_type>(), r2));
                }

            private:
                static XSIMD_INLINE batch_type negative(const batch_type& q, const batch_type& w) noexcept
                {
                    batch_type p = floor(q);
                    batch_type z = q - p;
                    auto test2 = z < batch_type(0.5);
                    z = select(test2, z - batch_type(1.), z);
                    z = q * sin(z, trigo_pi_tag());
                    return -log(constants::invpi<batch_type>() * abs(z)) - w;
                }

                static XSIMD_INLINE batch_type other(const batch_type& x) noexcept
                {
                    auto xlt650 = (x < batch_type(6.5));
                    batch_type r0x = x;
                    batch_type r0z = x;
                    batch_type r0s = batch_type(1.);
                    batch_type r1 = batch_type(0.);
                    batch_type p = constants::nan<batch_type>();
                    if (any(xlt650))
                    {
                        batch_type z = batch_type(1.);
                        batch_type tx = select(xlt650, x, batch_type(0.));
                        batch_type nx = batch_type(0.);
                        const batch_type _075 = batch_type(0.75);
                        const batch_type _150 = batch_type(1.50);
                        const batch_type _125 = batch_type(1.25);
                        const batch_type _250 = batch_type(2.50);
                        auto xge150 = (x >= _150);
                        auto txgt250 = (tx > _250);

                        // x >= 1.5
                        while (any(xge150 && txgt250))
                        {
                            nx = select(txgt250, nx - batch_type(1.), nx);
                            tx = select(txgt250, x + nx, tx);
                            z = select(txgt250, z * tx, z);
                            txgt250 = (tx > _250);
                        }
                        r0x = select(xge150, x + nx - batch_type(2.), x);
                        r0z = select(xge150, z, r0z);
                        r0s = select(xge150, batch_type(1.), r0s);

                        // x >= 1.25 && x < 1.5
                        auto xge125 = (x >= _125);
                        auto xge125t = xge125 && !xge150;
                        if (any(xge125))
                        {
                            r0x = select(xge125t, x - batch_type(1.), r0x);
                            r0z = select(xge125t, z * x, r0z);
                            r0s = select(xge125t, batch_type(-1.), r0s);
                        }

                        // x >= 0.75 && x < 1.5
                        batch_bool<float, A> kernelC(false);
                        auto xge075 = (x >= _075);
                        auto xge075t = xge075 && !xge125;
                        if (any(xge075t))
                        {
                            kernelC = xge075t;
                            r0x = select(xge075t, x - batch_type(1.), x);
                            r0z = select(xge075t, batch_type(1.), r0z);
                            r0s = select(xge075t, batch_type(-1.), r0s);
                            p = gammalnC(r0x);
                        }

                        // tx < 1.5 && x < 0.75
                        auto txlt150 = (tx < _150) && !xge075;
                        if (any(txlt150))
                        {
                            auto orig = txlt150;
                            while (any(txlt150))
                            {
                                z = select(txlt150, z * tx, z);
                                nx = select(txlt150, nx + batch_type(1.), nx);
                                tx = select(txlt150, x + nx, tx);
                                txlt150 = (tx < _150) && !xge075;
                            }
                            r0x = select(orig, r0x + nx - batch_type(2.), r0x);
                            r0z = select(orig, z, r0z);
                            r0s = select(orig, batch_type(-1.), r0s);
                        }
                        p = select(kernelC, p, gammalnB(r0x));
                        if (all(xlt650))
                            return fma(r0x, p, r0s * log(abs(r0z)));
                    }
                    r0z = select(xlt650, abs(r0z), x);
                    batch_type m = log(r0z);
                    r1 = fma(r0x, p, r0s * m);
                    batch_type r2 = fma(x - batch_type(0.5), m, constants::logsqrt2pi<batch_type>() - x);
                    r2 += gammaln2(batch_type(1.) / (x * x)) / x;
                    return select(xlt650, r1, r2);
                }
            };

            template <class A>
            struct lgamma_impl<batch<double, A>>
            {
                using batch_type = batch<double, A>;

                static XSIMD_INLINE batch_type compute(const batch_type& a) noexcept
                {
                    auto inf_result = (a <= batch_type(0.)) && is_flint(a);
                    batch_type x = select(inf_result, constants::nan<batch_type>(), a);
                    batch_type q = abs(x);
#ifndef XSIMD_NO_INFINITIES
                    inf_result = (q == constants::infinity<batch_type>());
#endif
                    auto test = (a < batch_type(-34.));
                    batch_type r = constants::nan<batch_type>();
                    if (any(test))
                    {
                        r = large_negative(q);
                        if (all(test))
                            return select(inf_result, constants::nan<batch_type>(), r);
                    }
                    batch_type r1 = other(a);
                    batch_type r2 = select(test, r, r1);
                    return select(a == constants::minusinfinity<batch_type>(), constants::nan<batch_type>(), select(inf_result, constants::infinity<batch_type>(), r2));
                }

            private:
                // FIXME: cannot mark this one as XSIMD_INLINE because there's a
                // recursive loop on `lgamma'.
                static inline batch_type large_negative(const batch_type& q) noexcept
                {
                    batch_type w = lgamma(q);
                    batch_type p = floor(q);
                    batch_type z = q - p;
                    auto test2 = (z < batch_type(0.5));
                    z = select(test2, z - batch_type(1.), z);
                    z = q * sin(z, trigo_pi_tag());
                    z = abs(z);
                    return constants::logpi<batch_type>() - log(z) - w;
                }

                static XSIMD_INLINE batch_type other(const batch_type& xx) noexcept
                {
                    batch_type x = xx;
                    auto test = (x < batch_type(13.));
                    batch_type r1 = batch_type(0.);
                    if (any(test))
                    {
                        batch_type z = batch_type(1.);
                        batch_type p = batch_type(0.);
                        batch_type u = select(test, x, batch_type(0.));
                        auto test1 = (u >= batch_type(3.));
                        while (any(test1))
                        {
                            p = select(test1, p - batch_type(1.), p);
                            u = select(test1, x + p, u);
                            z = select(test1, z * u, z);
                            test1 = (u >= batch_type(3.));
                        }

                        auto test2 = (u < batch_type(2.));
                        while (any(test2))
                        {
                            z = select(test2, z / u, z);
                            p = select(test2, p + batch_type(1.), p);
                            u = select(test2, x + p, u);
                            test2 = (u < batch_type(2.));
                        }

                        z = abs(z);
                        x += p - batch_type(2.);
                        r1 = x * gammaln1(x) + log(z);
                        if (all(test))
                            return r1;
                    }
                    batch_type r2 = fma(xx - batch_type(0.5), log(xx), constants::logsqrt2pi<batch_type>() - xx);
                    batch_type p = batch_type(1.) / (xx * xx);
                    r2 += gammalnA(p) / xx;
                    return select(test, r1, r2);
                }
            };
        }

        template <class A, class T>
        XSIMD_INLINE batch<T, A> lgamma(batch<T, A> const& self, requires_arch<generic>) noexcept
        {
            return detail::lgamma_impl<batch<T, A>>::compute(self);
        }

        // log
        /* origin: boost/simd/arch/common/simd/function/log.hpp */
        /*
         * ====================================================
         * copyright 2016 NumScale SAS
         *
         * Distributed under the Boost Software License, Version 1.0.
         * (See copy at http://boost.org/LICENSE_1_0.txt)
         * ====================================================
         */
        template <class A>
        XSIMD_INLINE batch<float, A> log(batch<float, A> const& self, requires_arch<generic>) noexcept
        {
            using batch_type = batch<float, A>;
            using int_type = as_integer_t<float>;
            using i_type = batch<int_type, A>;
            batch_type x = self;
            i_type k(0);
            auto isnez = (self != batch_type(0.));
#ifndef XSIMD_NO_DENORMALS
            auto test = (self < constants::smallestposval<batch_type>()) && isnez;
            if (any(test))
            {
                k = select(batch_bool_cast<int_type>(test), k - i_type(23), k);
                x = select(test, x * batch_type(8388608ul), x);
            }
#endif
            i_type ix = ::xsimd::bitwise_cast<int_type>(x);
            ix += 0x3f800000 - 0x3f3504f3;
            k += (ix >> 23) - 0x7f;
            ix = (ix & i_type(0x007fffff)) + 0x3f3504f3;
            x = ::xsimd::bitwise_cast<float>(ix);
            batch_type f = --x;
            batch_type s = f / (batch_type(2.) + f);
            batch_type z = s * s;
            batch_type w = z * z;
            batch_type t1 = w * detail::horner<batch_type, 0x3eccce13, 0x3e789e26>(w);
            batch_type t2 = z * detail::horner<batch_type, 0x3f2aaaaa, 0x3e91e9ee>(w);
            batch_type R = t2 + t1;
            batch_type hfsq = batch_type(0.5) * f * f;
            batch_type dk = to_float(k);
            batch_type r = fma(dk, constants::log_2hi<batch_type>(), fma(s, (hfsq + R), dk * constants::log_2lo<batch_type>()) - hfsq + f);
#ifndef XSIMD_NO_INFINITIES
            batch_type zz = select(isnez, select(self == constants::infinity<batch_type>(), constants::infinity<batch_type>(), r), constants::minusinfinity<batch_type>());
#else
            batch_type zz = select(isnez, r, constants::minusinfinity<batch_type>());
#endif
            return select(!(self >= batch_type(0.)), constants::nan<batch_type>(), zz);
        }

        template <class A>
        XSIMD_INLINE batch<double, A> log(batch<double, A> const& self, requires_arch<generic>) noexcept
        {
            using batch_type = batch<double, A>;
            using int_type = as_integer_t<double>;
            using i_type = batch<int_type, A>;

            batch_type x = self;
            i_type hx = ::xsimd::bitwise_cast<int_type>(x) >> 32;
            i_type k(0);
            auto isnez = (self != batch_type(0.));
#ifndef XSIMD_NO_DENORMALS
            auto test = (self < constants::smallestposval<batch_type>()) && isnez;
            if (any(test))
            {
                k = select(batch_bool_cast<int_type>(test), k - i_type(54), k);
                x = select(test, x * batch_type(18014398509481984ull), x);
            }
#endif
            hx += 0x3ff00000 - 0x3fe6a09e;
            k += (hx >> 20) - 0x3ff;
            batch_type dk = to_float(k);
            hx = (hx & i_type(0x000fffff)) + 0x3fe6a09e;
            x = ::xsimd::bitwise_cast<double>(hx << 32 | (i_type(0xffffffff) & ::xsimd::bitwise_cast<int_type>(x)));

            batch_type f = --x;
            batch_type hfsq = batch_type(0.5) * f * f;
            batch_type s = f / (batch_type(2.) + f);
            batch_type z = s * s;
            batch_type w = z * z;

            batch_type t1 = w * detail::horner<batch_type, 0x3fd999999997fa04ll, 0x3fcc71c51d8e78afll, 0x3fc39a09d078c69fll>(w);
            batch_type t2 = z * detail::horner<batch_type, 0x3fe5555555555593ll, 0x3fd2492494229359ll, 0x3fc7466496cb03dell, 0x3fc2f112df3e5244ll>(w);
            batch_type R = t2 + t1;
            batch_type r = fma(dk, constants::log_2hi<batch_type>(), fma(s, (hfsq + R), dk * constants::log_2lo<batch_type>()) - hfsq + f);
#ifndef XSIMD_NO_INFINITIES
            batch_type zz = select(isnez, select(self == constants::infinity<batch_type>(), constants::infinity<batch_type>(), r), constants::minusinfinity<batch_type>());
#else
            batch_type zz = select(isnez, r, constants::minusinfinity<batch_type>());
#endif
            return select(!(self >= batch_type(0.)), constants::nan<batch_type>(), zz);
        }

        template <class A, class T>
        XSIMD_INLINE batch<std::complex<T>, A> log(const batch<std::complex<T>, A>& z, requires_arch<generic>) noexcept
        {
            return batch<std::complex<T>, A>(log(abs(z)), atan2(z.imag(), z.real()));
        }

        // log2
        template <class A>
        XSIMD_INLINE batch<float, A> log2(batch<float, A> const& self, requires_arch<generic>) noexcept
        {
            using batch_type = batch<float, A>;
            using int_type = as_integer_t<float>;
            using i_type = batch<int_type, A>;
            batch_type x = self;
            i_type k(0);
            auto isnez = (self != batch_type(0.));
#ifndef XSIMD_NO_DENORMALS
            auto test = (self < constants::smallestposval<batch_type>()) && isnez;
            if (any(test))
            {
                k = select(batch_bool_cast<int_type>(test), k - i_type(25), k);
                x = select(test, x * batch_type(33554432ul), x);
            }
#endif
            i_type ix = ::xsimd::bitwise_cast<int_type>(x);
            ix += 0x3f800000 - 0x3f3504f3;
            k += (ix >> 23) - 0x7f;
            ix = (ix & i_type(0x007fffff)) + 0x3f3504f3;
            x = ::xsimd::bitwise_cast<float>(ix);
            batch_type f = --x;
            batch_type s = f / (batch_type(2.) + f);
            batch_type z = s * s;
            batch_type w = z * z;
            batch_type t1 = w * detail::horner<batch_type, 0x3eccce13, 0x3e789e26>(w);
            batch_type t2 = z * detail::horner<batch_type, 0x3f2aaaaa, 0x3e91e9ee>(w);
            batch_type R = t1 + t2;
            batch_type hfsq = batch_type(0.5) * f * f;
            batch_type dk = to_float(k);
            batch_type r = fma(fms(s, hfsq + R, hfsq) + f, constants::invlog_2<batch_type>(), dk);
#ifndef XSIMD_NO_INFINITIES
            batch_type zz = select(isnez, select(self == constants::infinity<batch_type>(), constants::infinity<batch_type>(), r), constants::minusinfinity<batch_type>());
#else
            batch_type zz = select(isnez, r, constants::minusinfinity<batch_type>());
#endif
            return select(!(self >= batch_type(0.)), constants::nan<batch_type>(), zz);
        }

        template <class A>
        XSIMD_INLINE batch<double, A> log2(batch<double, A> const& self, requires_arch<generic>) noexcept
        {
            using batch_type = batch<double, A>;
            using int_type = as_integer_t<double>;
            using i_type = batch<int_type, A>;
            batch_type x = self;
            i_type hx = ::xsimd::bitwise_cast<int_type>(x) >> 32;
            i_type k(0);
            auto isnez = (self != batch_type(0.));
#ifndef XSIMD_NO_DENORMALS
            auto test = (self < constants::smallestposval<batch_type>()) && isnez;
            if (any(test))
            {
                k = select(batch_bool_cast<typename i_type::value_type>(test), k - i_type(54), k);
                x = select(test, x * batch_type(18014398509481984ull), x);
            }
#endif
            hx += 0x3ff00000 - 0x3fe6a09e;
            k += (hx >> 20) - 0x3ff;
            hx = (hx & i_type(0x000fffff)) + 0x3fe6a09e;
            x = ::xsimd::bitwise_cast<double>(hx << 32 | (i_type(0xffffffff) & ::xsimd::bitwise_cast<int_type>(x)));
            batch_type f = --x;
            batch_type s = f / (batch_type(2.) + f);
            batch_type z = s * s;
            batch_type w = z * z;
            batch_type t1 = w * detail::horner<batch_type, 0x3fd999999997fa04ll, 0x3fcc71c51d8e78afll, 0x3fc39a09d078c69fll>(w);
            batch_type t2 = z * detail::horner<batch_type, 0x3fe5555555555593ll, 0x3fd2492494229359ll, 0x3fc7466496cb03dell, 0x3fc2f112df3e5244ll>(w);
            batch_type R = t2 + t1;
            batch_type hfsq = batch_type(0.5) * f * f;
            batch_type hi = f - hfsq;
            hi = hi & ::xsimd::bitwise_cast<double>((constants::allbits<i_type>() << 32));
            batch_type lo = fma(s, hfsq + R, f - hi - hfsq);
            batch_type val_hi = hi * constants::invlog_2hi<batch_type>();
            batch_type val_lo = fma(lo + hi, constants::invlog_2lo<batch_type>(), lo * constants::invlog_2hi<batch_type>());
            batch_type dk = to_float(k);
            batch_type w1 = dk + val_hi;
            val_lo += (dk - w1) + val_hi;
            val_hi = w1;
            batch_type r = val_lo + val_hi;
#ifndef XSIMD_NO_INFINITIES
            batch_type zz = select(isnez, select(self == constants::infinity<batch_type>(), constants::infinity<batch_type>(), r), constants::minusinfinity<batch_type>());
#else
            batch_type zz = select(isnez, r, constants::minusinfinity<batch_type>());
#endif
            return select(!(self >= batch_type(0.)), constants::nan<batch_type>(), zz);
        }

        namespace detail
        {
            template <class T, class A>
            XSIMD_INLINE batch<T, A> logN_complex_impl(const batch<T, A>& z, typename batch<T, A>::value_type base) noexcept
            {
                using batch_type = batch<T, A>;
                using rv_type = typename batch_type::value_type;
                return log(z) / batch_type(rv_type(base));
            }
        }

        template <class A, class T>
        XSIMD_INLINE batch<std::complex<T>, A> log2(batch<std::complex<T>, A> const& self, requires_arch<generic>) noexcept
        {
            return detail::logN_complex_impl(self, std::log(2));
        }

        // log10
        /* origin: FreeBSD /usr/src/lib/msun/src/e_log10f.c */
        /*
         * ====================================================
         * Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
         *
         * Developed at SunPro, a Sun Microsystems, Inc. business.
         * Permission to use, copy, modify, and distribute this
         * software is freely granted, provided that this notice
         * is preserved.
         * ====================================================
         */
        template <class A>
        XSIMD_INLINE batch<float, A> log10(batch<float, A> const& self, requires_arch<generic>) noexcept
        {
            using batch_type = batch<float, A>;
            const batch_type
                ivln10hi(4.3432617188e-01f),
                ivln10lo(-3.1689971365e-05f),
                log10_2hi(3.0102920532e-01f),
                log10_2lo(7.9034151668e-07f);
            using int_type = as_integer_t<float>;
            using i_type = batch<int_type, A>;
            batch_type x = self;
            i_type k(0);
            auto isnez = (self != batch_type(0.));
#ifndef XSIMD_NO_DENORMALS
            auto test = (self < constants::smallestposval<batch_type>()) && isnez;
            if (any(test))
            {
                k = select(batch_bool_cast<int_type>(test), k - i_type(25), k);
                x = select(test, x * batch_type(33554432ul), x);
            }
#endif
            i_type ix = ::xsimd::bitwise_cast<int_type>(x);
            ix += 0x3f800000 - 0x3f3504f3;
            k += (ix >> 23) - 0x7f;
            ix = (ix & i_type(0x007fffff)) + 0x3f3504f3;
            x = ::xsimd::bitwise_cast<float>(ix);
            batch_type f = --x;
            batch_type s = f / (batch_type(2.) + f);
            batch_type z = s * s;
            batch_type w = z * z;
            batch_type t1 = w * detail::horner<batch_type, 0x3eccce13, 0x3e789e26>(w);
            batch_type t2 = z * detail::horner<batch_type, 0x3f2aaaaa, 0x3e91e9ee>(w);
            batch_type R = t2 + t1;
            batch_type dk = to_float(k);
            batch_type hfsq = batch_type(0.5) * f * f;
            batch_type hibits = f - hfsq;
            hibits &= ::xsimd::bitwise_cast<float>(i_type(0xfffff000));
            batch_type lobits = fma(s, hfsq + R, f - hibits - hfsq);
            batch_type r = fma(dk, log10_2hi,
                               fma(hibits, ivln10hi,
                                   fma(lobits, ivln10hi,
                                       fma(lobits + hibits, ivln10lo, dk * log10_2lo))));
#ifndef XSIMD_NO_INFINITIES
            batch_type zz = select(isnez, select(self == constants::infinity<batch_type>(), constants::infinity<batch_type>(), r), constants::minusinfinity<batch_type>());
#else
            batch_type zz = select(isnez, r, constants::minusinfinity<batch_type>());
#endif
            return select(!(self >= batch_type(0.)), constants::nan<batch_type>(), zz);
        }

        template <class A>
        XSIMD_INLINE batch<double, A> log10(batch<double, A> const& self, requires_arch<generic>) noexcept
        {
            using batch_type = batch<double, A>;
            const batch_type
                ivln10hi(4.34294481878168880939e-01),
                ivln10lo(2.50829467116452752298e-11),
                log10_2hi(3.01029995663611771306e-01),
                log10_2lo(3.69423907715893078616e-13);
            using int_type = as_integer_t<double>;
            using i_type = batch<int_type, A>;
            batch_type x = self;
            i_type hx = ::xsimd::bitwise_cast<int_type>(x) >> 32;
            i_type k(0);
            auto isnez = (self != batch_type(0.));
#ifndef XSIMD_NO_DENORMALS
            auto test = (self < constants::smallestposval<batch_type>()) && isnez;
            if (any(test))
            {
                k = select(batch_bool_cast<int_type>(test), k - i_type(54), k);
                x = select(test, x * batch_type(18014398509481984ull), x);
            }
#endif
            hx += 0x3ff00000 - 0x3fe6a09e;
            k += (hx >> 20) - 0x3ff;
            hx = (hx & i_type(0x000fffff)) + 0x3fe6a09e;
            x = ::xsimd::bitwise_cast<double>(hx << 32 | (i_type(0xffffffff) & ::xsimd::bitwise_cast<int_type>(x)));
            batch_type f = --x;
            batch_type dk = to_float(k);
            batch_type s = f / (batch_type(2.) + f);
            batch_type z = s * s;
            batch_type w = z * z;
            batch_type t1 = w * detail::horner<batch_type, 0x3fd999999997fa04ll, 0x3fcc71c51d8e78afll, 0x3fc39a09d078c69fll>(w);
            batch_type t2 = z * detail::horner<batch_type, 0x3fe5555555555593ll, 0x3fd2492494229359ll, 0x3fc7466496cb03dell, 0x3fc2f112df3e5244ll>(w);
            batch_type R = t2 + t1;
            batch_type hfsq = batch_type(0.5) * f * f;
            batch_type hi = f - hfsq;
            hi = hi & ::xsimd::bitwise_cast<double>(constants::allbits<i_type>() << 32);
            batch_type lo = f - hi - hfsq + s * (hfsq + R);
            batch_type val_hi = hi * ivln10hi;
            batch_type y = dk * log10_2hi;
            batch_type val_lo = dk * log10_2lo + (lo + hi) * ivln10lo + lo * ivln10hi;
            batch_type w1 = y + val_hi;
            val_lo += (y - w1) + val_hi;
            val_hi = w1;
            batch_type r = val_lo + val_hi;
#ifndef XSIMD_NO_INFINITIES
            batch_type zz = select(isnez, select(self == constants::infinity<batch_type>(), constants::infinity<batch_type>(), r), constants::minusinfinity<batch_type>());
#else
            batch_type zz = select(isnez, r, constants::minusinfinity<batch_type>());
#endif
            return select(!(self >= batch_type(0.)), constants::nan<batch_type>(), zz);
        }

        template <class A, class T>
        XSIMD_INLINE batch<std::complex<T>, A> log10(const batch<std::complex<T>, A>& z, requires_arch<generic>) noexcept
        {
            return detail::logN_complex_impl(z, std::log(10));
        }

        // log1p
        /* origin: boost/simd/arch/common/simd/function/log1p.hpp */
        /*
         * ====================================================
         * copyright 2016 NumScale SAS
         *
         * Distributed under the Boost Software License, Version 1.0.
         * (See copy at http://boost.org/LICENSE_1_0.txt)
         * ====================================================
         */
        template <class A>
        XSIMD_INLINE batch<float, A> log1p(batch<float, A> const& self, requires_arch<generic>) noexcept
        {
            using batch_type = batch<float, A>;
            using int_type = as_integer_t<float>;
            using i_type = batch<int_type, A>;
            const batch_type uf = self + batch_type(1.);
            auto isnez = (uf != batch_type(0.));
            i_type iu = ::xsimd::bitwise_cast<int_type>(uf);
            iu += 0x3f800000 - 0x3f3504f3;
            i_type k = (iu >> 23) - 0x7f;
            iu = (iu & i_type(0x007fffff)) + 0x3f3504f3;
            batch_type f = --(::xsimd::bitwise_cast<float>(iu));
            batch_type s = f / (batch_type(2.) + f);
            batch_type z = s * s;
            batch_type w = z * z;
            batch_type t1 = w * detail::horner<batch_type, 0x3eccce13, 0x3e789e26>(w);
            batch_type t2 = z * detail::horner<batch_type, 0x3f2aaaaa, 0x3e91e9ee>(w);
            batch_type R = t2 + t1;
            batch_type hfsq = batch_type(0.5) * f * f;
            batch_type dk = to_float(k);
            /* correction term ~ log(1+x)-log(u), avoid underflow in c/u */
            batch_type c = select(batch_bool_cast<float>(k >= i_type(2)), batch_type(1.) - (uf - self), self - (uf - batch_type(1.))) / uf;
            batch_type r = fma(dk, constants::log_2hi<batch_type>(), fma(s, (hfsq + R), dk * constants::log_2lo<batch_type>() + c) - hfsq + f);
#ifndef XSIMD_NO_INFINITIES
            batch_type zz = select(isnez, select(self == constants::infinity<batch_type>(), constants::infinity<batch_type>(), r), constants::minusinfinity<batch_type>());
#else
            batch_type zz = select(isnez, r, constants::minusinfinity<batch_type>());
#endif
            return select(!(uf >= batch_type(0.)), constants::nan<batch_type>(), zz);
        }

        template <class A>
        XSIMD_INLINE batch<double, A> log1p(batch<double, A> const& self, requires_arch<generic>) noexcept
        {
            using batch_type = batch<double, A>;
            using int_type = as_integer_t<double>;
            using i_type = batch<int_type, A>;
            const batch_type uf = self + batch_type(1.);
            auto isnez = (uf != batch_type(0.));
            i_type hu = ::xsimd::bitwise_cast<int_type>(uf) >> 32;
            hu += 0x3ff00000 - 0x3fe6a09e;
            i_type k = (hu >> 20) - 0x3ff;
            /* correction term ~ log(1+x)-log(u), avoid underflow in c/u */
            batch_type c = select(batch_bool_cast<double>(k >= i_type(2)), batch_type(1.) - (uf - self), self - (uf - batch_type(1.))) / uf;
            hu = (hu & i_type(0x000fffff)) + 0x3fe6a09e;
            batch_type f = ::xsimd::bitwise_cast<double>((hu << 32) | (i_type(0xffffffff) & ::xsimd::bitwise_cast<int_type>(uf)));
            f = --f;
            batch_type hfsq = batch_type(0.5) * f * f;
            batch_type s = f / (batch_type(2.) + f);
            batch_type z = s * s;
            batch_type w = z * z;
            batch_type t1 = w * detail::horner<batch_type, 0x3fd999999997fa04ll, 0x3fcc71c51d8e78afll, 0x3fc39a09d078c69fll>(w);
            batch_type t2 = z * detail::horner<batch_type, 0x3fe5555555555593ll, 0x3fd2492494229359ll, 0x3fc7466496cb03dell, 0x3fc2f112df3e5244ll>(w);
            batch_type R = t2 + t1;
            batch_type dk = to_float(k);
            batch_type r = fma(dk, constants::log_2hi<batch_type>(), fma(s, hfsq + R, dk * constants::log_2lo<batch_type>() + c) - hfsq + f);
#ifndef XSIMD_NO_INFINITIES
            batch_type zz = select(isnez, select(self == constants::infinity<batch_type>(), constants::infinity<batch_type>(), r), constants::minusinfinity<batch_type>());
#else
            batch_type zz = select(isnez, r, constants::minusinfinity<batch_type>());
#endif
            return select(!(uf >= batch_type(0.)), constants::nan<batch_type>(), zz);
        }

        template <class A, class T>
        XSIMD_INLINE batch<std::complex<T>, A> log1p(batch<std::complex<T>, A> const& self, requires_arch<generic>) noexcept
        {
            using batch_type = batch<std::complex<T>, A>;
            using real_batch = typename batch_type::real_batch;
            batch_type u = 1 + self;
            batch_type logu = log(u);
            return select(u == batch_type(1.),
                          self,
                          select(u.real() <= real_batch(0.),
                                 logu,
                                 logu * self / (u - batch_type(1.))));
        }

        // mod
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> mod(batch<T, A> const& self, batch<T, A> const& other, requires_arch<generic>) noexcept
        {
            return detail::apply([](T x, T y) noexcept -> T
                                 { return x % y; },
                                 self, other);
        }

        // nearbyint
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> nearbyint(batch<T, A> const& self, requires_arch<generic>) noexcept
        {
            return self;
        }
        namespace detail
        {
            template <class A, class T>
            XSIMD_INLINE batch<T, A> nearbyintf(batch<T, A> const& self) noexcept
            {
                using batch_type = batch<T, A>;
                batch_type s = bitofsign(self);
                batch_type v = self ^ s;
                batch_type t2n = constants::twotonmb<batch_type>();
                // Under fast-math, reordering is possible and the compiler optimizes d
                // to v. That's not what we want, so prevent compiler optimization here.
                // FIXME: it may be better to emit a memory barrier here (?).
#ifdef __FAST_MATH__
                volatile batch_type d0 = v + t2n;
                batch_type d = *(batch_type*)(void*)(&d0) - t2n;
#else
                batch_type d0 = v + t2n;
                batch_type d = d0 - t2n;
#endif
                return s ^ select(v < t2n, d, v);
            }
        }
        template <class A>
        XSIMD_INLINE batch<float, A> nearbyint(batch<float, A> const& self, requires_arch<generic>) noexcept
        {
            return detail::nearbyintf(self);
        }
        template <class A>
        XSIMD_INLINE batch<double, A> nearbyint(batch<double, A> const& self, requires_arch<generic>) noexcept
        {
            return detail::nearbyintf(self);
        }

        // nearbyint_as_int
        template <class T, class A, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> nearbyint_as_int(batch<T, A> const& self, requires_arch<generic>) noexcept
        {
            return self;
        }

        // nearbyint_as_int
        template <class A>
        XSIMD_INLINE batch<as_integer_t<float>, A>
        nearbyint_as_int(batch<float, A> const& self, requires_arch<generic>) noexcept
        {
            using U = as_integer_t<float>;
            return kernel::detail::apply_transform<U>([](float x) noexcept -> U
                                                      { return std::nearbyintf(x); },
                                                      self);
        }

        template <class A>
        XSIMD_INLINE batch<as_integer_t<double>, A>
        nearbyint_as_int(batch<double, A> const& self, requires_arch<generic>) noexcept
        {
            using U = as_integer_t<double>;
            return kernel::detail::apply_transform<U>([](double x) noexcept -> U
                                                      { return std::nearbyint(x); },
                                                      self);
        }

        // nextafter
        namespace detail
        {
            template <class T, class A, bool is_int = std::is_integral<T>::value>
            struct nextafter_kernel
            {
                using batch_type = batch<T, A>;

                static XSIMD_INLINE batch_type next(batch_type const& b) noexcept
                {
                    return b;
                }

                static XSIMD_INLINE batch_type prev(batch_type const& b) noexcept
                {
                    return b;
                }
            };

            template <class T, class A>
            struct bitwise_cast_batch;

            template <class A>
            struct bitwise_cast_batch<float, A>
            {
                using type = batch<int32_t, A>;
            };

            template <class A>
            struct bitwise_cast_batch<double, A>
            {
                using type = batch<int64_t, A>;
            };

            template <class T, class A>
            struct nextafter_kernel<T, A, false>
            {
                using batch_type = batch<T, A>;
                using int_batch = typename bitwise_cast_batch<T, A>::type;
                using int_type = typename int_batch::value_type;

                static XSIMD_INLINE batch_type next(const batch_type& b) noexcept
                {
                    batch_type n = ::xsimd::bitwise_cast<T>(::xsimd::bitwise_cast<int_type>(b) + int_type(1));
                    return select(b == constants::infinity<batch_type>(), b, n);
                }

                static XSIMD_INLINE batch_type prev(const batch_type& b) noexcept
                {
                    batch_type p = ::xsimd::bitwise_cast<T>(::xsimd::bitwise_cast<int_type>(b) - int_type(1));
                    return select(b == constants::minusinfinity<batch_type>(), b, p);
                }
            };
        }
        template <class A, class T>
        XSIMD_INLINE batch<T, A> nextafter(batch<T, A> const& from, batch<T, A> const& to, requires_arch<generic>) noexcept
        {
            using kernel = detail::nextafter_kernel<T, A>;
            return select(from == to, from,
                          select(to > from, kernel::next(from), kernel::prev(from)));
        }

        // pow
        /* origin: boost/simd/arch/common/simd/function/pow.hpp*/
        /*
         * ====================================================
         * copyright 2016 NumScale SAS
         *
         * Distributed under the Boost Software License, Version 1.0.
         * (See copy at http://boost.org/LICENSE_1_0.txt)
         * ====================================================
         */
        template <class A, class T>
        XSIMD_INLINE batch<T, A> pow(batch<T, A> const& self, batch<T, A> const& other, requires_arch<generic>) noexcept
        {
            using batch_type = batch<T, A>;
            const auto zero = batch_type(0.);
            auto negself = self < zero;
            auto iszeropowpos = self == zero && other >= zero;
            auto adj_self = select(iszeropowpos, batch_type(1), abs(self));
            batch_type z = exp(other * log(adj_self));
            z = select(iszeropowpos, zero, z);
            z = select(is_odd(other) && negself, -z, z);
            auto invalid = negself && !(is_flint(other) || isinf(other));
            return select(invalid, constants::nan<batch_type>(), z);
        }

        template <class A, class T>
        XSIMD_INLINE batch<std::complex<T>, A> pow(const batch<std::complex<T>, A>& a, const batch<std::complex<T>, A>& z, requires_arch<generic>) noexcept
        {
            using cplx_batch = batch<std::complex<T>, A>;
            using real_batch = typename cplx_batch::real_batch;
            real_batch absa = abs(a);
            real_batch arga = arg(a);
            real_batch x = z.real();
            real_batch y = z.imag();
            real_batch r = pow(absa, x);
            real_batch theta = x * arga;
            real_batch ze(0);
            auto cond = (y == ze);
            r = select(cond, r, r * exp(-y * arga));
            theta = select(cond, theta, theta + y * log(absa));
            return select(absa == ze, cplx_batch(ze), cplx_batch(r * cos(theta), r * sin(theta)));
        }

        // reciprocal
        template <class T, class A, class = typename std::enable_if<std::is_floating_point<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> reciprocal(batch<T, A> const& self,
                                            requires_arch<generic>) noexcept
        {
            using batch_type = batch<T, A>;
            return div(batch_type(1), self);
        }

        // reduce_add
        template <class A, class T>
        XSIMD_INLINE std::complex<T> reduce_add(batch<std::complex<T>, A> const& self, requires_arch<generic>) noexcept
        {
            return { reduce_add(self.real()), reduce_add(self.imag()) };
        }

        namespace detail
        {
            template <class T, T N>
            struct split_high
            {
                static constexpr T get(T i, T)
                {
                    return i >= N ? (i % 2) : i + N;
                }
            };

            template <class Op, class A, class T>
            XSIMD_INLINE T reduce(Op, batch<T, A> const& self, std::integral_constant<unsigned, 1>) noexcept
            {
                return self.get(0);
            }

            template <class Op, class A, class T, unsigned Lvl>
            XSIMD_INLINE T reduce(Op op, batch<T, A> const& self, std::integral_constant<unsigned, Lvl>) noexcept
            {
                using index_type = as_unsigned_integer_t<T>;
                batch<T, A> split = swizzle(self, make_batch_constant<index_type, A, split_high<index_type, Lvl / 2>>());
                return reduce(op, op(split, self), std::integral_constant<unsigned, Lvl / 2>());
            }
        }

        // reduce_max
        template <class A, class T>
        XSIMD_INLINE T reduce_max(batch<T, A> const& self, requires_arch<generic>) noexcept
        {
            return detail::reduce([](batch<T, A> const& x, batch<T, A> const& y)
                                  { return max(x, y); },
                                  self, std::integral_constant<unsigned, batch<T, A>::size>());
        }

        // reduce_min
        template <class A, class T>
        XSIMD_INLINE T reduce_min(batch<T, A> const& self, requires_arch<generic>) noexcept
        {
            return detail::reduce([](batch<T, A> const& x, batch<T, A> const& y)
                                  { return min(x, y); },
                                  self, std::integral_constant<unsigned, batch<T, A>::size>());
        }

        // remainder
        template <class A>
        XSIMD_INLINE batch<float, A> remainder(batch<float, A> const& self, batch<float, A> const& other, requires_arch<generic>) noexcept
        {
            return fnma(nearbyint(self / other), other, self);
        }
        template <class A>
        XSIMD_INLINE batch<double, A> remainder(batch<double, A> const& self, batch<double, A> const& other, requires_arch<generic>) noexcept
        {
            return fnma(nearbyint(self / other), other, self);
        }
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> remainder(batch<T, A> const& self, batch<T, A> const& other, requires_arch<generic>) noexcept
        {
            auto mod = self % other;
            return select(mod <= other / 2, mod, mod - other);
        }

        // select
        template <class A, class T>
        XSIMD_INLINE batch<std::complex<T>, A> select(batch_bool<T, A> const& cond, batch<std::complex<T>, A> const& true_br, batch<std::complex<T>, A> const& false_br, requires_arch<generic>) noexcept
        {
            return { select(cond, true_br.real(), false_br.real()), select(cond, true_br.imag(), false_br.imag()) };
        }

        // sign
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> sign(batch<T, A> const& self, requires_arch<generic>) noexcept
        {
            using batch_type = batch<T, A>;
            batch_type res = select(self > batch_type(0), batch_type(1), batch_type(0)) - select(self < batch_type(0), batch_type(1), batch_type(0));
            return res;
        }

        namespace detail
        {
            template <class T, class A>
            XSIMD_INLINE batch<T, A> signf(batch<T, A> const& self) noexcept
            {
                using batch_type = batch<T, A>;
                batch_type res = select(self > batch_type(0.f), batch_type(1.f), batch_type(0.f)) - select(self < batch_type(0.f), batch_type(1.f), batch_type(0.f));
#ifdef XSIMD_NO_NANS
                return res;
#else
                return select(isnan(self), constants::nan<batch_type>(), res);
#endif
            }
        }

        template <class A>
        XSIMD_INLINE batch<float, A> sign(batch<float, A> const& self, requires_arch<generic>) noexcept
        {
            return detail::signf(self);
        }
        template <class A>
        XSIMD_INLINE batch<double, A> sign(batch<double, A> const& self, requires_arch<generic>) noexcept
        {
            return detail::signf(self);
        }
        template <class A, class T>
        XSIMD_INLINE batch<std::complex<T>, A> sign(const batch<std::complex<T>, A>& z, requires_arch<generic>) noexcept
        {
            using batch_type = batch<std::complex<T>, A>;
            using real_batch = typename batch_type::real_batch;
            auto rz = z.real();
            auto iz = z.imag();
            return select(rz != real_batch(0.),
                          batch_type(sign(rz)),
                          batch_type(sign(iz)));
        }

        // signnz
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> signnz(batch<T, A> const& self, requires_arch<generic>) noexcept
        {
            using batch_type = batch<T, A>;
            return (self >> (sizeof(T) * 8 - 1)) | batch_type(1.);
        }

        namespace detail
        {
            template <class T, class A>
            XSIMD_INLINE batch<T, A> signnzf(batch<T, A> const& self) noexcept
            {
                using batch_type = batch<T, A>;
#ifndef XSIMD_NO_NANS
                return select(isnan(self), constants::nan<batch_type>(), batch_type(1.) | (constants::signmask<batch_type>() & self));
#else
                return batch_type(1.) | (constants::signmask<batch_type>() & self);
#endif
            }
        }

        template <class A>
        XSIMD_INLINE batch<float, A> signnz(batch<float, A> const& self, requires_arch<generic>) noexcept
        {
            return detail::signnzf(self);
        }
        template <class A>
        XSIMD_INLINE batch<double, A> signnz(batch<double, A> const& self, requires_arch<generic>) noexcept
        {
            return detail::signnzf(self);
        }

        // sqrt
        template <class A, class T>
        XSIMD_INLINE batch<std::complex<T>, A> sqrt(batch<std::complex<T>, A> const& z, requires_arch<generic>) noexcept
        {

            constexpr T csqrt_scale_factor = std::is_same<T, float>::value ? 6.7108864e7f : 1.8014398509481984e16;
            constexpr T csqrt_scale = std::is_same<T, float>::value ? 1.220703125e-4f : 7.450580596923828125e-9;
            using batch_type = batch<std::complex<T>, A>;
            using real_batch = batch<T, A>;
            real_batch x = z.real();
            real_batch y = z.imag();
            real_batch sqrt_x = sqrt(fabs(x));
            real_batch sqrt_hy = sqrt(0.5 * fabs(y));
            auto cond = (fabs(x) > real_batch(4.) || fabs(y) > real_batch(4.));
            x = select(cond, x * 0.25, x * csqrt_scale_factor);
            y = select(cond, y * 0.25, y * csqrt_scale_factor);
            real_batch scale = select(cond, real_batch(2.), real_batch(csqrt_scale));
            real_batch r = abs(batch_type(x, y));

            auto condxp = x > real_batch(0.);
            real_batch t0 = select(condxp, xsimd::sqrt(0.5 * (r + x)), xsimd::sqrt(0.5 * (r - x)));
            real_batch r0 = scale * fabs((0.5 * y) / t0);
            t0 *= scale;
            real_batch t = select(condxp, t0, r0);
            r = select(condxp, r0, t0);
            batch_type resg = select(y < real_batch(0.), batch_type(t, -r), batch_type(t, r));
            real_batch ze(0.);

            return select(y == ze,
                          select(x == ze,
                                 batch_type(ze, ze),
                                 select(x < ze, batch_type(ze, sqrt_x), batch_type(sqrt_x, ze))),
                          select(x == ze,
                                 select(y > ze, batch_type(sqrt_hy, sqrt_hy), batch_type(sqrt_hy, -sqrt_hy)),
                                 resg));
        }

        // tgamma

        namespace detail
        {
            /* origin: boost/simd/arch/common/detail/generic/stirling_kernel.hpp */
            /*
             * ====================================================
             * copyright 2016 NumScale SAS
             *
             * Distributed under the Boost Software License, Version 1.0.
             * (See copy at http://boost.org/LICENSE_1_0.txt)
             * ====================================================
             */
            template <class B>
            struct stirling_kernel;

            template <class A>
            struct stirling_kernel<batch<float, A>>
            {
                using batch_type = batch<float, A>;
                static XSIMD_INLINE batch_type compute(const batch_type& x) noexcept
                {
                    return horner<batch_type,
                                  0x3daaaaab,
                                  0x3b638e39,
                                  0xbb2fb930,
                                  0xb970b359>(x);
                }

                static XSIMD_INLINE batch_type split_limit() noexcept
                {
                    return batch_type(bit_cast<float>(uint32_t(0x41d628f6)));
                }

                static XSIMD_INLINE batch_type large_limit() noexcept
                {
                    return batch_type(bit_cast<float>(uint32_t(0x420c28f3)));
                }
            };

            template <class A>
            struct stirling_kernel<batch<double, A>>
            {
                using batch_type = batch<double, A>;
                static XSIMD_INLINE batch_type compute(const batch_type& x) noexcept
                {
                    return horner<batch_type,
                                  0x3fb5555555555986ull, //   8.33333333333482257126E-2
                                  0x3f6c71c71b98c5fdull, //   3.47222221605458667310E-3
                                  0xbf65f72607d44fd7ull, //  -2.68132617805781232825E-3
                                  0xbf2e166b27e61d7cull, //  -2.29549961613378126380E-4
                                  0x3f49cc72592d7293ull //   7.87311395793093628397E-4
                                  >(x);
                }

                static XSIMD_INLINE batch_type split_limit() noexcept
                {
                    return batch_type(bit_cast<double>(uint64_t(0x4061e083ba3443d4)));
                }

                static XSIMD_INLINE batch_type large_limit() noexcept
                {
                    return batch_type(bit_cast<double>(uint64_t(0x4065800000000000)));
                }
            };

            /* origin: boost/simd/arch/common/simd/function/stirling.hpp */
            /*
             * ====================================================
             * copyright 2016 NumScale SAS
             *
             * Distributed under the Boost Software License, Version 1.0.
             * (See copy at http://boost.org/LICENSE_1_0.txt)
             * ====================================================
             */
            template <class T, class A>
            XSIMD_INLINE batch<T, A> stirling(const batch<T, A>& a) noexcept
            {
                using batch_type = batch<T, A>;
                const batch_type stirlingsplitlim = stirling_kernel<batch_type>::split_limit();
                const batch_type stirlinglargelim = stirling_kernel<batch_type>::large_limit();
                batch_type x = select(a >= batch_type(0.), a, constants::nan<batch_type>());
                batch_type w = batch_type(1.) / x;
                w = fma(w, stirling_kernel<batch_type>::compute(w), batch_type(1.));
                batch_type y = exp(-x);
                auto test = (x < stirlingsplitlim);
                batch_type z = x - batch_type(0.5);
                z = select(test, z, batch_type(0.5) * z);
                batch_type v = exp(z * log(abs(x)));
                y *= v;
                y = select(test, y, y * v);
                y *= constants::sqrt_2pi<batch_type>() * w;
#ifndef XSIMD_NO_INFINITIES
                y = select(isinf(x), x, y);
#endif
                return select(x > stirlinglargelim, constants::infinity<batch_type>(), y);
            }

            /* origin: boost/simd/arch/common/detail/generic/gamma_kernel.hpp */
            /*
             * ====================================================
             * copyright 2016 NumScale SAS
             *
             * Distributed under the Boost Software License, Version 1.0.
             * (See copy at http://boost.org/LICENSE_1_0.txt)
             * ====================================================
             */
            template <class B>
            struct tgamma_kernel;

            template <class A>
            struct tgamma_kernel<batch<float, A>>
            {
                using batch_type = batch<float, A>;
                static XSIMD_INLINE batch_type compute(const batch_type& x) noexcept
                {
                    return horner<batch_type,
                                  0x3f800000UL, //  9.999999757445841E-01
                                  0x3ed87799UL, //  4.227874605370421E-01
                                  0x3ed2d411UL, //  4.117741948434743E-01
                                  0x3da82a34UL, //  8.211174403261340E-02
                                  0x3d93ae7cUL, //  7.211014349068177E-02
                                  0x3b91db14UL, //  4.451165155708328E-03
                                  0x3ba90c99UL, //  5.158972571345137E-03
                                  0x3ad28b22UL //  1.606319369134976E-03
                                  >(x);
                }
            };

            template <class A>
            struct tgamma_kernel<batch<double, A>>
            {
                using batch_type = batch<double, A>;
                static XSIMD_INLINE batch_type compute(const batch_type& x) noexcept
                {
                    return horner<batch_type,
                                  0x3ff0000000000000ULL, // 9.99999999999999996796E-1
                                  0x3fdfa1373993e312ULL, // 4.94214826801497100753E-1
                                  0x3fca8da9dcae7d31ULL, // 2.07448227648435975150E-1
                                  0x3fa863d918c423d3ULL, // 4.76367800457137231464E-2
                                  0x3f8557cde9db14b0ULL, // 1.04213797561761569935E-2
                                  0x3f5384e3e686bfabULL, // 1.19135147006586384913E-3
                                  0x3f24fcb839982153ULL // 1.60119522476751861407E-4
                                  >(x)
                        / horner<batch_type,
                                 0x3ff0000000000000ULL, //  1.00000000000000000320E00
                                 0x3fb24944c9cd3c51ULL, //  7.14304917030273074085E-2
                                 0xbfce071a9d4287c2ULL, // -2.34591795718243348568E-1
                                 0x3fa25779e33fde67ULL, //  3.58236398605498653373E-2
                                 0x3f8831ed5b1bb117ULL, //  1.18139785222060435552E-2
                                 0xBf7240e4e750b44aULL, // -4.45641913851797240494E-3
                                 0x3f41ae8a29152573ULL, //  5.39605580493303397842E-4
                                 0xbef8487a8400d3aFULL // -2.31581873324120129819E-5
                                 >(x);
                }
            };

            /* origin: boost/simd/arch/common/simd/function/gamma.hpp */
            /*
             * ====================================================
             * copyright 2016 NumScale SAS
             *
             * Distributed under the Boost Software License, Version 1.0.
             * (See copy at http://boost.org/LICENSE_1_0.txt)
             * ====================================================
             */
            template <class B>
            XSIMD_INLINE B tgamma_large_negative(const B& a) noexcept
            {
                B st = stirling(a);
                B p = floor(a);
                B sgngam = select(is_even(p), -B(1.), B(1.));
                B z = a - p;
                auto test2 = z < B(0.5);
                z = select(test2, z - B(1.), z);
                z = a * sin(z, trigo_pi_tag());
                z = abs(z);
                return sgngam * constants::pi<B>() / (z * st);
            }

            template <class B, class BB>
            XSIMD_INLINE B tgamma_other(const B& a, const BB& test) noexcept
            {
                B x = select(test, B(2.), a);
#ifndef XSIMD_NO_INFINITIES
                auto inf_result = (a == constants::infinity<B>());
                x = select(inf_result, B(2.), x);
#endif
                B z = B(1.);
                auto test1 = (x >= B(3.));
                while (any(test1))
                {
                    x = select(test1, x - B(1.), x);
                    z = select(test1, z * x, z);
                    test1 = (x >= B(3.));
                }
                test1 = (x < B(0.));
                while (any(test1))
                {
                    z = select(test1, z / x, z);
                    x = select(test1, x + B(1.), x);
                    test1 = (x < B(0.));
                }
                auto test2 = (x < B(2.));
                while (any(test2))
                {
                    z = select(test2, z / x, z);
                    x = select(test2, x + B(1.), x);
                    test2 = (x < B(2.));
                }
                x = z * tgamma_kernel<B>::compute(x - B(2.));
#ifndef XSIMD_NO_INFINITIES
                return select(inf_result, a, x);
#else
                return x;
#endif
            }
        }

        template <class A, class T>
        XSIMD_INLINE batch<T, A> tgamma(batch<T, A> const& self, requires_arch<generic>) noexcept
        {
            using batch_type = batch<T, A>;
            auto nan_result = (self < batch_type(0.) && is_flint(self));
#ifndef XSIMD_NO_INVALIDS
            nan_result = isnan(self) || nan_result;
#endif
            batch_type q = abs(self);
            auto test = (self < batch_type(-33.));
            batch_type r = constants::nan<batch_type>();
            if (any(test))
            {
                r = detail::tgamma_large_negative(q);
                if (all(test))
                    return select(nan_result, constants::nan<batch_type>(), r);
            }
            batch_type r1 = detail::tgamma_other(self, test);
            batch_type r2 = select(test, r, r1);
            return select(self == batch_type(0.), copysign(constants::infinity<batch_type>(), self), select(nan_result, constants::nan<batch_type>(), r2));
        }

    }

}

#endif
