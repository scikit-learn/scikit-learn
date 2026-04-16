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

#ifndef XSIMD_GENERIC_TRIGO_HPP
#define XSIMD_GENERIC_TRIGO_HPP

#include "./xsimd_generic_details.hpp"

#include <array>

namespace xsimd
{

    namespace kernel
    {
        /* origin: boost/simd/arch/common/detail/simd/trig_base.hpp */
        /*
         * ====================================================
         * copyright 2016 NumScale SAS
         *
         * Distributed under the Boost Software License, Version 1.0.
         * (See copy at http://boost.org/LICENSE_1_0.txt)
         * ====================================================
         */

        using namespace types;

        // acos
        template <class A, class T>
        XSIMD_INLINE batch<T, A> acos(batch<T, A> const& self, requires_arch<generic>) noexcept
        {
            using batch_type = batch<T, A>;
            batch_type x = abs(self);
            auto x_larger_05 = x > batch_type(0.5);
            x = select(x_larger_05, sqrt(fma(batch_type(-0.5), x, batch_type(0.5))), self);
            x = asin(x);
            x = select(x_larger_05, x + x, x);
            x = select(self < batch_type(-0.5), constants::pi<batch_type>() - x, x);
            return select(x_larger_05, x, constants::pio2<batch_type>() - x);
        }
        template <class A, class T>
        XSIMD_INLINE batch<std::complex<T>, A> acos(const batch<std::complex<T>, A>& z, requires_arch<generic>) noexcept
        {
            using batch_type = batch<std::complex<T>, A>;
            using real_batch = typename batch_type::real_batch;
            batch_type tmp = asin(z);
            return { constants::pio2<real_batch>() - tmp.real(), -tmp.imag() };
        }

        // acosh
        /* origin: boost/simd/arch/common/simd/function/acosh.hpp */
        /*
         * ====================================================
         * copyright 2016 NumScale SAS
         *
         * Distributed under the Boost Software License, Version 1.0.
         * (See copy at http://boost.org/LICENSE_1_0.txt)
         * ====================================================
         */
        template <class A, class T>
        XSIMD_INLINE batch<T, A> acosh(batch<T, A> const& self, requires_arch<generic>) noexcept
        {
            using batch_type = batch<T, A>;
            batch_type x = self - batch_type(1.);
            auto test = x > constants::oneotwoeps<batch_type>();
            batch_type z = select(test, self, x + sqrt(x + x + x * x));
            batch_type l1pz = log1p(z);
            return select(test, l1pz + constants::log_2<batch_type>(), l1pz);
        }
        template <class A, class T>
        XSIMD_INLINE batch<std::complex<T>, A> acosh(const batch<std::complex<T>, A>& z, requires_arch<generic>) noexcept
        {
            using batch_type = batch<std::complex<T>, A>;
            batch_type w = acos(z);
            w = batch_type(-w.imag(), w.real());
            return w;
        }

        // asin
        template <class A>
        XSIMD_INLINE batch<float, A> asin(batch<float, A> const& self, requires_arch<generic>) noexcept
        {
            using batch_type = batch<float, A>;
            batch_type x = abs(self);
            batch_type sign = bitofsign(self);
            auto x_larger_05 = x > batch_type(0.5);
            batch_type z = select(x_larger_05, batch_type(0.5) * (batch_type(1.) - x), x * x);
            x = select(x_larger_05, sqrt(z), x);
            batch_type z1 = detail::horner<batch_type,
                                           0x3e2aaae4,
                                           0x3d9980f6,
                                           0x3d3a3ec7,
                                           0x3cc617e3,
                                           0x3d2cb352>(z);
            z1 = fma(z1, z * x, x);
            z = select(x_larger_05, constants::pio2<batch_type>() - (z1 + z1), z1);
            return z ^ sign;
        }
        template <class A>
        XSIMD_INLINE batch<double, A> asin(batch<double, A> const& self, requires_arch<generic>) noexcept
        {
            using batch_type = batch<double, A>;
            batch_type x = abs(self);
            auto small_cond = x < constants::sqrteps<batch_type>();
            batch_type ct1 = batch_type(bit_cast<double>(int64_t(0x3fe4000000000000)));
            batch_type zz1 = batch_type(1.) - x;
            batch_type vp = zz1 * detail::horner<batch_type, 0x403c896240f3081dull, 0xc03991aaac01ab68ull, 0x401bdff5baf33e6aull, 0xbfe2079259f9290full, 0x3f684fc3988e9f08ull>(zz1) / detail::horner1<batch_type, 0x40756709b0b644beull, 0xc077fe08959063eeull, 0x40626219af6a7f42ull, 0xc035f2a2b6bf5d8cull>(zz1);
            zz1 = sqrt(zz1 + zz1);
            batch_type z = constants::pio4<batch_type>() - zz1;
            zz1 = fms(zz1, vp, constants::pio_2lo<batch_type>());
            z = z - zz1;
            zz1 = z + constants::pio4<batch_type>();
            batch_type zz2 = self * self;
            z = zz2 * detail::horner<batch_type, 0xc020656c06ceafd5ull, 0x40339007da779259ull, 0xc0304331de27907bull, 0x4015c74b178a2dd9ull, 0xbfe34341333e5c16ull, 0x3f716b9b0bd48ad3ull>(zz2) / detail::horner1<batch_type, 0xc04898220a3607acull, 0x4061705684ffbf9dull, 0xc06265bb6d3576d7ull, 0x40519fc025fe9054ull, 0xc02d7b590b5e0eabull>(zz2);
            zz2 = fma(x, z, x);
            return select(x > batch_type(1.), constants::nan<batch_type>(),
                          select(small_cond, x,
                                 select(x > ct1, zz1, zz2))
                              ^ bitofsign(self));
        }
        template <class A, class T>
        XSIMD_INLINE batch<std::complex<T>, A> asin(const batch<std::complex<T>, A>& z, requires_arch<generic>) noexcept
        {
            using batch_type = batch<std::complex<T>, A>;
            using real_batch = typename batch_type::real_batch;
            real_batch x = z.real();
            real_batch y = z.imag();

            batch_type ct(-y, x);
            batch_type zz(real_batch(1.) - (x - y) * (x + y), -2 * x * y);
            zz = log(ct + sqrt(zz));
            batch_type resg(zz.imag(), -zz.real());

            return select(y == real_batch(0.),
                          select(fabs(x) > real_batch(1.),
                                 batch_type(constants::pio2<real_batch>(), real_batch(0.)),
                                 batch_type(asin(x), real_batch(0.))),
                          resg);
        }

        // asinh
        /* origin: boost/simd/arch/common/simd/function/asinh.hpp */
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
            template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
            XSIMD_INLINE batch<T, A>
            average(const batch<T, A>& x1, const batch<T, A>& x2) noexcept
            {
                return (x1 & x2) + ((x1 ^ x2) >> 1);
            }

            template <class A, class T>
            XSIMD_INLINE batch<T, A>
            averagef(const batch<T, A>& x1, const batch<T, A>& x2) noexcept
            {
                using batch_type = batch<T, A>;
                return fma(x1, batch_type(0.5), x2 * batch_type(0.5));
            }
            template <class A>
            XSIMD_INLINE batch<float, A> average(batch<float, A> const& x1, batch<float, A> const& x2) noexcept
            {
                return averagef(x1, x2);
            }
            template <class A>
            XSIMD_INLINE batch<double, A> average(batch<double, A> const& x1, batch<double, A> const& x2) noexcept
            {
                return averagef(x1, x2);
            }
        }
        template <class A>
        XSIMD_INLINE batch<float, A> asinh(batch<float, A> const& self, requires_arch<generic>) noexcept
        {
            using batch_type = batch<float, A>;
            batch_type x = abs(self);
            auto lthalf = x < batch_type(0.5);
            batch_type x2 = x * x;
            batch_type bts = bitofsign(self);
            batch_type z(0.);
            if (any(lthalf))
            {
                z = detail::horner<batch_type,
                                   0x3f800000,
                                   0xbe2aa9ad,
                                   0x3d9949b1,
                                   0xbd2ee581,
                                   0x3ca4d6e6>(x2)
                    * x;
                if (all(lthalf))
                    return z ^ bts;
            }
            batch_type tmp = select(x > constants::oneosqrteps<batch_type>(), x, detail::average(x, hypot(batch_type(1.), x)));
#ifndef XSIMD_NO_NANS
            return select(isnan(self), constants::nan<batch_type>(), select(lthalf, z, log(tmp) + constants::log_2<batch_type>()) ^ bts);
#else
            return select(lthalf, z, log(tmp) + constants::log_2<batch_type>()) ^ bts;
#endif
        }
        template <class A>
        XSIMD_INLINE batch<double, A> asinh(batch<double, A> const& self, requires_arch<generic>) noexcept
        {
            using batch_type = batch<double, A>;
            batch_type x = abs(self);
            auto test = x > constants::oneosqrteps<batch_type>();
            batch_type z = select(test, x - batch_type(1.), x + x * x / (batch_type(1.) + hypot(batch_type(1.), x)));
#ifndef XSIMD_NO_INFINITIES
            z = select(x == constants::infinity<batch_type>(), x, z);
#endif
            batch_type l1pz = log1p(z);
            z = select(test, l1pz + constants::log_2<batch_type>(), l1pz);
            return bitofsign(self) ^ z;
        }
        template <class A, class T>
        XSIMD_INLINE batch<std::complex<T>, A> asinh(const batch<std::complex<T>, A>& z, requires_arch<generic>) noexcept
        {
            using batch_type = batch<std::complex<T>, A>;
            batch_type w = asin(batch_type(-z.imag(), z.real()));
            w = batch_type(w.imag(), -w.real());
            return w;
        }

        // atan
        namespace detail
        {
            template <class A>
            static XSIMD_INLINE batch<float, A> kernel_atan(const batch<float, A>& x, const batch<float, A>& recx) noexcept
            {
                using batch_type = batch<float, A>;
                const auto flag1 = x < constants::tan3pio8<batch_type>();
                const auto flag2 = (x >= batch_type(bit_cast<float>((uint32_t)0x3ed413cd))) && flag1;
                batch_type yy = select(flag1, batch_type(0.), constants::pio2<batch_type>());
                yy = select(flag2, constants::pio4<batch_type>(), yy);
                batch_type xx = select(flag1, x, -recx);
                xx = select(flag2, (x - batch_type(1.)) / (x + batch_type(1.)), xx);
                const batch_type z = xx * xx;
                batch_type z1 = detail::horner<batch_type,
                                               0xbeaaaa2aul,
                                               0x3e4c925ful,
                                               0xbe0e1b85ul,
                                               0x3da4f0d1ul>(z);
                z1 = fma(xx, z1 * z, xx);
                z1 = select(flag2, z1 + constants::pio_4lo<batch_type>(), z1);
                z1 = select(!flag1, z1 + constants::pio_2lo<batch_type>(), z1);
                return yy + z1;
            }
            template <class A>
            static XSIMD_INLINE batch<double, A> kernel_atan(const batch<double, A>& x, const batch<double, A>& recx) noexcept
            {
                using batch_type = batch<double, A>;
                const auto flag1 = x < constants::tan3pio8<batch_type>();
                const auto flag2 = (x >= constants::tanpio8<batch_type>()) && flag1;
                batch_type yy = select(flag1, batch_type(0.), constants::pio2<batch_type>());
                yy = select(flag2, constants::pio4<batch_type>(), yy);
                batch_type xx = select(flag1, x, -recx);
                xx = select(flag2, (x - batch_type(1.)) / (x + batch_type(1.)), xx);
                batch_type z = xx * xx;
                z *= detail::horner<batch_type,
                                    0xc0503669fd28ec8eull,
                                    0xc05eb8bf2d05ba25ull,
                                    0xc052c08c36880273ull,
                                    0xc03028545b6b807aull,
                                    0xbfec007fa1f72594ull>(z)
                    / detail::horner1<batch_type,
                                      0x4068519efbbd62ecull,
                                      0x407e563f13b049eaull,
                                      0x407b0e18d2e2be3bull,
                                      0x4064a0dd43b8fa25ull,
                                      0x4038dbc45b14603cull>(z);
                z = fma(xx, z, xx);
                z = select(flag2, z + constants::pio_4lo<batch_type>(), z);
                z = z + select(flag1, batch_type(0.), constants::pio_2lo<batch_type>());
                return yy + z;
            }
        }
        template <class A, class T>
        XSIMD_INLINE batch<T, A> atan(batch<T, A> const& self, requires_arch<generic>) noexcept
        {
            using batch_type = batch<T, A>;
            const batch_type absa = abs(self);
            const batch_type x = detail::kernel_atan(absa, batch_type(1.) / absa);
            return x ^ bitofsign(self);
        }
        template <class A, class T>
        XSIMD_INLINE batch<std::complex<T>, A> atan(const batch<std::complex<T>, A>& z, requires_arch<generic>) noexcept
        {
            using batch_type = batch<std::complex<T>, A>;
            using real_batch = typename batch_type::real_batch;
            real_batch x = z.real();
            real_batch y = z.imag();
            real_batch x2 = x * x;
            real_batch one(1.);
            real_batch a = one - x2 - (y * y);
            real_batch w = 0.5 * atan2(2. * x, a);
            real_batch num = y + one;
            num = x2 + num * num;
            real_batch den = y - one;
            den = x2 + den * den;
            batch_type res = select((x == real_batch(0.)) && (y == real_batch(1.)),
                                    batch_type(real_batch(0.), constants::infinity<real_batch>()),
                                    batch_type(w, 0.25 * log(num / den)));
            return res;
        }

        // atanh
        /* origin: boost/simd/arch/common/simd/function/acosh.hpp */
        /*
         * ====================================================
         * copyright 2016 NumScale SAS
         *
         * Distributed under the Boost Software License, Version 1.0.
         * (See copy at http://boost.org/LICENSE_1_0.txt)
         * ====================================================
         */
        template <class A, class T>
        XSIMD_INLINE batch<T, A> atanh(batch<T, A> const& self, requires_arch<generic>) noexcept
        {
            using batch_type = batch<T, A>;
            batch_type x = abs(self);
            batch_type t = x + x;
            batch_type z = batch_type(1.) - x;
            auto test = x < batch_type(0.5);
            batch_type tmp = select(test, x, t) / z;
            return bitofsign(self) ^ (batch_type(0.5) * log1p(select(test, fma(t, tmp, t), tmp)));
        }
        template <class A, class T>
        XSIMD_INLINE batch<std::complex<T>, A> atanh(const batch<std::complex<T>, A>& z, requires_arch<generic>) noexcept
        {
            using batch_type = batch<std::complex<T>, A>;
            batch_type w = atan(batch_type(-z.imag(), z.real()));
            w = batch_type(w.imag(), -w.real());
            return w;
        }

        // atan2
        template <class A, class T>
        XSIMD_INLINE batch<T, A> atan2(batch<T, A> const& self, batch<T, A> const& other, requires_arch<generic>) noexcept
        {
            using batch_type = batch<T, A>;
            const batch_type q = abs(self / other);
            const batch_type z = detail::kernel_atan(q, batch_type(1.) / q);
            return select(other > batch_type(0.), z, constants::pi<batch_type>() - z) * signnz(self);
        }

        // cos
        namespace detail
        {
            template <class T, class A>
            XSIMD_INLINE batch<T, A> quadrant(const batch<T, A>& x) noexcept
            {
                return x & batch<T, A>(3);
            }

            template <class A>
            XSIMD_INLINE batch<float, A> quadrant(const batch<float, A>& x) noexcept
            {
                return to_float(quadrant(to_int(x)));
            }

            template <class A>
            XSIMD_INLINE batch<double, A> quadrant(const batch<double, A>& x) noexcept
            {
                using batch_type = batch<double, A>;
                batch_type a = x * batch_type(0.25);
                return (a - floor(a)) * batch_type(4.);
            }
            /* origin: boost/simd/arch/common/detail/simd/f_trig_evaluation.hpp */
            /*
             * ====================================================
             * copyright 2016 NumScale SAS
             *
             * Distributed under the Boost Software License, Version 1.0.
             * (See copy at http://boost.org/LICENSE_1_0.txt)
             * ====================================================
             */

            template <class A>
            XSIMD_INLINE batch<float, A> cos_eval(const batch<float, A>& z) noexcept
            {
                using batch_type = batch<float, A>;
                batch_type y = detail::horner<batch_type,
                                              0x3d2aaaa5,
                                              0xbab60619,
                                              0x37ccf5ce>(z);
                return batch_type(1.) + fma(z, batch_type(-0.5), y * z * z);
            }

            template <class A>
            XSIMD_INLINE batch<float, A> sin_eval(const batch<float, A>& z, const batch<float, A>& x) noexcept
            {
                using batch_type = batch<float, A>;
                batch_type y = detail::horner<batch_type,
                                              0xbe2aaaa2,
                                              0x3c08839d,
                                              0xb94ca1f9>(z);
                return fma(y * z, x, x);
            }

            template <class A>
            static XSIMD_INLINE batch<float, A> base_tancot_eval(const batch<float, A>& z) noexcept
            {
                using batch_type = batch<float, A>;
                batch_type zz = z * z;
                batch_type y = detail::horner<batch_type,
                                              0x3eaaaa6f,
                                              0x3e0896dd,
                                              0x3d5ac5c9,
                                              0x3cc821b5,
                                              0x3b4c779c,
                                              0x3c19c53b>(zz);
                return fma(y, zz * z, z);
            }

            template <class A, class BB>
            static XSIMD_INLINE batch<float, A> tan_eval(const batch<float, A>& z, const BB& test) noexcept
            {
                using batch_type = batch<float, A>;
                batch_type y = base_tancot_eval(z);
                return select(test, y, -batch_type(1.) / y);
            }

            template <class A, class BB>
            static XSIMD_INLINE batch<float, A> cot_eval(const batch<float, A>& z, const BB& test) noexcept
            {
                using batch_type = batch<float, A>;
                batch_type y = base_tancot_eval(z);
                return select(test, batch_type(1.) / y, -y);
            }

            /* origin: boost/simd/arch/common/detail/simd/d_trig_evaluation.hpp */
            /*
             * ====================================================
             * copyright 2016 NumScale SAS
             *
             * Distributed under the Boost Software License, Version 1.0.
             * (See copy at http://boost.org/LICENSE_1_0.txt)
             * ====================================================
             */
            template <class A>
            static XSIMD_INLINE batch<double, A> cos_eval(const batch<double, A>& z) noexcept
            {
                using batch_type = batch<double, A>;
                batch_type y = detail::horner<batch_type,
                                              0x3fe0000000000000ull,
                                              0xbfa5555555555551ull,
                                              0x3f56c16c16c15d47ull,
                                              0xbefa01a019ddbcd9ull,
                                              0x3e927e4f8e06d9a5ull,
                                              0xbe21eea7c1e514d4ull,
                                              0x3da8ff831ad9b219ull>(z);
                return batch_type(1.) - y * z;
            }

            template <class A>
            static XSIMD_INLINE batch<double, A> sin_eval(const batch<double, A>& z, const batch<double, A>& x) noexcept
            {
                using batch_type = batch<double, A>;
                batch_type y = detail::horner<batch_type,
                                              0xbfc5555555555548ull,
                                              0x3f8111111110f7d0ull,
                                              0xbf2a01a019bfdf03ull,
                                              0x3ec71de3567d4896ull,
                                              0xbe5ae5e5a9291691ull,
                                              0x3de5d8fd1fcf0ec1ull>(z);
                return fma(y * z, x, x);
            }

            template <class A>
            static XSIMD_INLINE batch<double, A> base_tancot_eval(const batch<double, A>& z) noexcept
            {
                using batch_type = batch<double, A>;
                batch_type zz = z * z;
                batch_type num = detail::horner<batch_type,
                                                0xc1711fead3299176ull,
                                                0x413199eca5fc9dddull,
                                                0xc0c992d8d24f3f38ull>(zz);
                batch_type den = detail::horner1<batch_type,
                                                 0xc189afe03cbe5a31ull,
                                                 0x4177d98fc2ead8efull,
                                                 0xc13427bc582abc96ull,
                                                 0x40cab8a5eeb36572ull>(zz);
                return fma(z, (zz * (num / den)), z);
            }

            template <class A, class BB>
            static XSIMD_INLINE batch<double, A> tan_eval(const batch<double, A>& z, const BB& test) noexcept
            {
                using batch_type = batch<double, A>;
                batch_type y = base_tancot_eval(z);
                return select(test, y, -batch_type(1.) / y);
            }

            template <class A, class BB>
            static XSIMD_INLINE batch<double, A> cot_eval(const batch<double, A>& z, const BB& test) noexcept
            {
                using batch_type = batch<double, A>;
                batch_type y = base_tancot_eval(z);
                return select(test, batch_type(1.) / y, -y);
            }
            /* origin: boost/simd/arch/common/detail/simd/trig_reduction.hpp */
            /*
             * ====================================================
             * copyright 2016 NumScale SAS
             *
             * Distributed under the Boost Software License, Version 1.0.
             * (See copy at http://boost.org/LICENSE_1_0.txt)
             * ====================================================
             */

            struct trigo_radian_tag
            {
            };
            struct trigo_pi_tag
            {
            };

            template <class B, class Tag = trigo_radian_tag>
            struct trigo_reducer
            {
                static XSIMD_INLINE B reduce(const B& x, B& xr) noexcept
                {
                    if (all(x <= constants::pio4<B>()))
                    {
                        xr = x;
                        return B(0.);
                    }
                    else if (all(x <= constants::pio2<B>()))
                    {
                        auto test = x > constants::pio4<B>();
                        xr = x - constants::pio2_1<B>();
                        xr -= constants::pio2_2<B>();
                        xr -= constants::pio2_3<B>();
                        xr = select(test, xr, x);
                        return select(test, B(1.), B(0.));
                    }
                    else if (all(x <= constants::twentypi<B>()))
                    {
                        B xi = nearbyint(x * constants::twoopi<B>());
                        xr = fnma(xi, constants::pio2_1<B>(), x);
                        xr -= xi * constants::pio2_2<B>();
                        xr -= xi * constants::pio2_3<B>();
                        return quadrant(xi);
                    }
                    else if (all(x <= constants::mediumpi<B>()))
                    {
                        B fn = nearbyint(x * constants::twoopi<B>());
                        B r = x - fn * constants::pio2_1<B>();
                        B w = fn * constants::pio2_1t<B>();
                        B t = r;
                        w = fn * constants::pio2_2<B>();
                        r = t - w;
                        w = fn * constants::pio2_2t<B>() - ((t - r) - w);
                        t = r;
                        w = fn * constants::pio2_3<B>();
                        r = t - w;
                        w = fn * constants::pio2_3t<B>() - ((t - r) - w);
                        xr = r - w;
                        return quadrant(fn);
                    }
                    else
                    {
                        static constexpr std::size_t size = B::size;
                        using value_type = typename B::value_type;
                        alignas(B) std::array<value_type, size> tmp;
                        alignas(B) std::array<value_type, size> txr;
                        alignas(B) std::array<value_type, size> args;
                        x.store_aligned(args.data());

                        for (std::size_t i = 0; i < size; ++i)
                        {
                            double arg = args[i];
                            if (arg == std::numeric_limits<value_type>::infinity())
                            {
                                tmp[i] = 0.;
                                txr[i] = std::numeric_limits<value_type>::quiet_NaN();
                            }
                            else
                            {
                                double y[2];
                                std::int32_t n = ::xsimd::detail::__ieee754_rem_pio2(arg, y);
                                tmp[i] = value_type(n & 3);
                                txr[i] = value_type(y[0]);
                            }
                        }
                        xr = B::load_aligned(&txr[0]);
                        B res = B::load_aligned(&tmp[0]);
                        return res;
                    }
                }
            };

            template <class B>
            struct trigo_reducer<B, trigo_pi_tag>
            {
                static XSIMD_INLINE B reduce(const B& x, B& xr) noexcept
                {
                    B xi = nearbyint(x * B(2.));
                    B x2 = x - xi * B(0.5);
                    xr = x2 * constants::pi<B>();
                    return quadrant(xi);
                }
            };

        }
        template <class A, class T>
        XSIMD_INLINE batch<T, A> cos(batch<T, A> const& self, requires_arch<generic>) noexcept
        {
            using batch_type = batch<T, A>;
            const batch_type x = abs(self);
            batch_type xr = constants::nan<batch_type>();
            const batch_type n = detail::trigo_reducer<batch_type>::reduce(x, xr);
            auto tmp = select(n >= batch_type(2.), batch_type(1.), batch_type(0.));
            auto swap_bit = fma(batch_type(-2.), tmp, n);
            auto sign_bit = select((swap_bit ^ tmp) != batch_type(0.), constants::signmask<batch_type>(), batch_type(0.));
            const batch_type z = xr * xr;
            const batch_type se = detail::sin_eval(z, xr);
            const batch_type ce = detail::cos_eval(z);
            const batch_type z1 = select(swap_bit != batch_type(0.), se, ce);
            return z1 ^ sign_bit;
        }

        template <class A, class T>
        XSIMD_INLINE batch<std::complex<T>, A> cos(batch<std::complex<T>, A> const& z, requires_arch<generic>) noexcept
        {
            return { cos(z.real()) * cosh(z.imag()), -sin(z.real()) * sinh(z.imag()) };
        }

        // cosh

        /* origin: boost/simd/arch/common/simd/function/cosh.hpp */
        /*
         * ====================================================
         * copyright 2016 NumScale SAS
         *
         * Distributed under the Boost Software License, Version 1.0.
         * (See copy at http://boost.org/LICENSE_1_0.txt)
         * ====================================================
         */

        template <class A, class T>
        XSIMD_INLINE batch<T, A> cosh(batch<T, A> const& self, requires_arch<generic>) noexcept
        {
            using batch_type = batch<T, A>;
            batch_type x = abs(self);
            auto test1 = x > (constants::maxlog<batch_type>() - constants::log_2<batch_type>());
            batch_type fac = select(test1, batch_type(0.5), batch_type(1.));
            batch_type tmp = exp(x * fac);
            batch_type tmp1 = batch_type(0.5) * tmp;
            return select(test1, tmp1 * tmp, detail::average(tmp, batch_type(1.) / tmp));
        }
        template <class A, class T>
        XSIMD_INLINE batch<std::complex<T>, A> cosh(const batch<std::complex<T>, A>& z, requires_arch<generic>) noexcept
        {
            auto x = z.real();
            auto y = z.imag();
            return { cosh(x) * cos(y), sinh(x) * sin(y) };
        }

        // sin
        namespace detail
        {
            template <class A, class T, class Tag = trigo_radian_tag>
            XSIMD_INLINE batch<T, A> sin(batch<T, A> const& self, Tag = Tag()) noexcept
            {
                using batch_type = batch<T, A>;
                const batch_type x = abs(self);
                batch_type xr = constants::nan<batch_type>();
                const batch_type n = detail::trigo_reducer<batch_type, Tag>::reduce(x, xr);
                auto tmp = select(n >= batch_type(2.), batch_type(1.), batch_type(0.));
                auto swap_bit = fma(batch_type(-2.), tmp, n);
                auto sign_bit = bitofsign(self) ^ select(tmp != batch_type(0.), constants::signmask<batch_type>(), batch_type(0.));
                const batch_type z = xr * xr;
                const batch_type se = detail::sin_eval(z, xr);
                const batch_type ce = detail::cos_eval(z);
                const batch_type z1 = select(swap_bit == batch_type(0.), se, ce);
                return z1 ^ sign_bit;
            }
        }

        template <class A, class T>
        XSIMD_INLINE batch<T, A> sin(batch<T, A> const& self, requires_arch<generic>) noexcept
        {
            return detail::sin(self);
        }

        template <class A, class T>
        XSIMD_INLINE batch<std::complex<T>, A> sin(batch<std::complex<T>, A> const& z, requires_arch<generic>) noexcept
        {
            return { sin(z.real()) * cosh(z.imag()), cos(z.real()) * sinh(z.imag()) };
        }

        // sincos
        template <class A, class T>
        XSIMD_INLINE std::pair<batch<T, A>, batch<T, A>> sincos(batch<T, A> const& self, requires_arch<generic>) noexcept
        {
            using batch_type = batch<T, A>;
            const batch_type x = abs(self);
            batch_type xr = constants::nan<batch_type>();
            const batch_type n = detail::trigo_reducer<batch_type>::reduce(x, xr);
            auto tmp = select(n >= batch_type(2.), batch_type(1.), batch_type(0.));
            auto swap_bit = fma(batch_type(-2.), tmp, n);
            const batch_type z = xr * xr;
            const batch_type se = detail::sin_eval(z, xr);
            const batch_type ce = detail::cos_eval(z);
            auto sin_sign_bit = bitofsign(self) ^ select(tmp != batch_type(0.), constants::signmask<batch_type>(), batch_type(0.));
            const batch_type sin_z1 = select(swap_bit == batch_type(0.), se, ce);
            auto cos_sign_bit = select((swap_bit ^ tmp) != batch_type(0.), constants::signmask<batch_type>(), batch_type(0.));
            const batch_type cos_z1 = select(swap_bit != batch_type(0.), se, ce);
            return std::make_pair(sin_z1 ^ sin_sign_bit, cos_z1 ^ cos_sign_bit);
        }

        template <class A, class T>
        XSIMD_INLINE std::pair<batch<std::complex<T>, A>, batch<std::complex<T>, A>>
        sincos(batch<std::complex<T>, A> const& z, requires_arch<generic>) noexcept
        {
            using batch_type = batch<std::complex<T>, A>;
            using real_batch = typename batch_type::real_batch;
            real_batch rcos = cos(z.real());
            real_batch rsin = sin(z.real());
            real_batch icosh = cosh(z.imag());
            real_batch isinh = sinh(z.imag());
            return std::make_pair(batch_type(rsin * icosh, rcos * isinh), batch_type(rcos * icosh, -rsin * isinh));
        }

        // sinh
        namespace detail
        {
            /* origin: boost/simd/arch/common/detail/generic/sinh_kernel.hpp */
            /*
             * ====================================================
             * copyright 2016 NumScale SAS
             *
             * Distributed under the Boost Software License, Version 1.0.
             * (See copy at http://boost.org/LICENSE_1_0.txt)
             * ====================================================
             */
            template <class A>
            XSIMD_INLINE batch<float, A> sinh_kernel(batch<float, A> const& self) noexcept
            {
                using batch_type = batch<float, A>;
                batch_type sqr_self = self * self;
                return detail::horner<batch_type,
                                      0x3f800000, // 1.0f
                                      0x3e2aaacc, // 1.66667160211E-1f
                                      0x3c087bbe, // 8.33028376239E-3f
                                      0x39559e2f // 2.03721912945E-4f
                                      >(sqr_self)
                    * self;
            }

            template <class A>
            XSIMD_INLINE batch<double, A> sinh_kernel(batch<double, A> const& self) noexcept
            {
                using batch_type = batch<double, A>;
                batch_type sqrself = self * self;
                return fma(self, (detail::horner<batch_type,
                                                 0xc115782bdbf6ab05ull, //  -3.51754964808151394800E5
                                                 0xc0c694b8c71d6182ull, //  -1.15614435765005216044E4,
                                                 0xc064773a398ff4feull, //  -1.63725857525983828727E2,
                                                 0xbfe9435fe8bb3cd6ull //  -7.89474443963537015605E-1
                                                 >(sqrself)
                                  / detail::horner1<batch_type,
                                                    0xc1401a20e4f90044ull, //  -2.11052978884890840399E6
                                                    0x40e1a7ba7ed72245ull, //   3.61578279834431989373E4,
                                                    0xc0715b6096e96484ull //  -2.77711081420602794433E2,
                                                    >(sqrself))
                               * sqrself,
                           self);
            }
        }
        /* origin: boost/simd/arch/common/simd/function/sinh.hpp */
        /*
         * ====================================================
         * copyright 2016 NumScale SAS
         *
         * Distributed under the Boost Software License, Version 1.0.
         * (See copy at http://boost.org/LICENSE_1_0.txt)
         * ====================================================
         */
        template <class A, class T>
        XSIMD_INLINE batch<T, A> sinh(batch<T, A> const& a, requires_arch<generic>) noexcept
        {
            using batch_type = batch<T, A>;
            batch_type half(0.5);
            batch_type x = abs(a);
            auto lt1 = x < batch_type(1.);
            batch_type bts = bitofsign(a);
            batch_type z(0.);
            if (any(lt1))
            {
                z = detail::sinh_kernel(x);
                if (all(lt1))
                    return z ^ bts;
            }
            auto test1 = x > (constants::maxlog<batch_type>() - constants::log_2<batch_type>());
            batch_type fac = select(test1, half, batch_type(1.));
            batch_type tmp = exp(x * fac);
            batch_type tmp1 = half * tmp;
            batch_type r = select(test1, tmp1 * tmp, tmp1 - half / tmp);
            return select(lt1, z, r) ^ bts;
        }
        template <class A, class T>
        XSIMD_INLINE batch<std::complex<T>, A> sinh(const batch<std::complex<T>, A>& z, requires_arch<generic>) noexcept
        {
            auto x = z.real();
            auto y = z.imag();
            return { sinh(x) * cos(y), cosh(x) * sin(y) };
        }

        // tan
        template <class A, class T>
        XSIMD_INLINE batch<T, A> tan(batch<T, A> const& self, requires_arch<generic>) noexcept
        {
            using batch_type = batch<T, A>;
            const batch_type x = abs(self);
            batch_type xr = constants::nan<batch_type>();
            const batch_type n = detail::trigo_reducer<batch_type>::reduce(x, xr);
            auto tmp = select(n >= batch_type(2.), batch_type(1.), batch_type(0.));
            auto swap_bit = fma(batch_type(-2.), tmp, n);
            auto test = (swap_bit == batch_type(0.));
            const batch_type y = detail::tan_eval(xr, test);
            return y ^ bitofsign(self);
        }
        template <class A, class T>
        XSIMD_INLINE batch<std::complex<T>, A> tan(batch<std::complex<T>, A> const& z, requires_arch<generic>) noexcept
        {
            using batch_type = batch<std::complex<T>, A>;
            using real_batch = typename batch_type::real_batch;
            real_batch d = cos(2 * z.real()) + cosh(2 * z.imag());
            batch_type winf(constants::infinity<real_batch>(), constants::infinity<real_batch>());
            real_batch wreal = sin(2 * z.real()) / d;
            real_batch wimag = sinh(2 * z.imag());
            batch_type wres = select(isinf(wimag), batch_type(wreal, real_batch(1.)), batch_type(wreal, wimag / d));
            return select(d == real_batch(0.), winf, wres);
        }

        // tanh
        namespace detail
        {
            /* origin: boost/simd/arch/common/detail/generic/tanh_kernel.hpp */
            /*
             * ====================================================
             * copyright 2016 NumScale SAS
             *
             * Distributed under the Boost Software License, Version 1.0.
             * (See copy at http://boost.org/LICENSE_1_0.txt)
             * ====================================================
             */
            template <class B>
            struct tanh_kernel;

            template <class A>
            struct tanh_kernel<batch<float, A>>
            {
                using batch_type = batch<float, A>;
                static XSIMD_INLINE batch_type tanh(const batch_type& x) noexcept
                {
                    batch_type sqrx = x * x;
                    return fma(detail::horner<batch_type,
                                              0xbeaaaa99, //    -3.33332819422E-1F
                                              0x3e088393, //    +1.33314422036E-1F
                                              0xbd5c1e2d, //    -5.37397155531E-2F
                                              0x3ca9134e, //    +2.06390887954E-2F
                                              0xbbbaf0ea //    -5.70498872745E-3F
                                              >(sqrx)
                                   * sqrx,
                               x, x);
                }

                static XSIMD_INLINE batch_type cotanh(const batch_type& x) noexcept
                {
                    return batch_type(1.) / tanh(x);
                }
            };

            template <class A>
            struct tanh_kernel<batch<double, A>>
            {
                using batch_type = batch<double, A>;
                static XSIMD_INLINE batch_type tanh(const batch_type& x) noexcept
                {
                    batch_type sqrx = x * x;
                    return fma(sqrx * p(sqrx) / q(sqrx), x, x);
                }

                static XSIMD_INLINE batch_type cotanh(const batch_type& x) noexcept
                {
                    batch_type sqrx = x * x;
                    batch_type qval = q(sqrx);
                    return qval / (x * fma(p(sqrx), sqrx, qval));
                }

                static XSIMD_INLINE batch_type p(const batch_type& x) noexcept
                {
                    return detail::horner<batch_type,
                                          0xc0993ac030580563, // -1.61468768441708447952E3
                                          0xc058d26a0e26682d, // -9.92877231001918586564E1,
                                          0xbfeedc5baafd6f4b // -9.64399179425052238628E-1
                                          >(x);
                }

                static XSIMD_INLINE batch_type q(const batch_type& x) noexcept
                {
                    return detail::horner1<batch_type,
                                           0x40b2ec102442040c, //  4.84406305325125486048E3
                                           0x40a176fa0e5535fa, //  2.23548839060100448583E3,
                                           0x405c33f28a581B86 //  1.12811678491632931402E2,
                                           >(x);
                }
            };

        }
        /* origin: boost/simd/arch/common/simd/function/tanh.hpp */
        /*
         * ====================================================
         * copyright 2016 NumScale SAS
         *
         * Distributed under the Boost Software License, Version 1.0.
         * (See copy at http://boost.org/LICENSE_1_0.txt)
         * ====================================================
         */
        template <class A, class T>
        XSIMD_INLINE batch<T, A> tanh(batch<T, A> const& self, requires_arch<generic>) noexcept
        {
            using batch_type = batch<T, A>;
            batch_type one(1.);
            batch_type x = abs(self);
            auto test = x < (batch_type(5.) / batch_type(8.));
            batch_type bts = bitofsign(self);
            batch_type z = one;
            if (any(test))
            {
                z = detail::tanh_kernel<batch_type>::tanh(x);
                if (all(test))
                    return z ^ bts;
            }
            batch_type r = fma(batch_type(-2.), one / (one + exp(x + x)), one);
            return select(test, z, r) ^ bts;
        }
        template <class A, class T>
        XSIMD_INLINE batch<std::complex<T>, A> tanh(const batch<std::complex<T>, A>& z, requires_arch<generic>) noexcept
        {
            using real_batch = typename batch<std::complex<T>, A>::real_batch;
            auto x = z.real();
            auto y = z.imag();
            real_batch two(2);
            auto d = cosh(two * x) + cos(two * y);
            return { sinh(two * x) / d, sin(two * y) / d };
        }

    }

}

#endif
