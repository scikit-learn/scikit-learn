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

#ifndef XSIMD_NUMERICAL_CONSTANT_HPP
#define XSIMD_NUMERICAL_CONSTANT_HPP

#include <limits>

#include "../types/xsimd_utils.hpp"

namespace xsimd
{

    namespace constants
    {

#define XSIMD_DEFINE_CONSTANT(NAME, SINGLE, DOUBLE) \
    template <class T>                              \
    XSIMD_INLINE T NAME() noexcept                  \
    {                                               \
        return T(NAME<typename T::value_type>());   \
    }                                               \
    template <>                                     \
    XSIMD_INLINE float NAME<float>() noexcept       \
    {                                               \
        return SINGLE;                              \
    }                                               \
    template <>                                     \
    XSIMD_INLINE double NAME<double>() noexcept     \
    {                                               \
        return DOUBLE;                              \
    }

#define XSIMD_DEFINE_CONSTANT_HEX(NAME, SINGLE, DOUBLE) \
    template <class T>                                  \
    XSIMD_INLINE T NAME() noexcept                      \
    {                                                   \
        return T(NAME<typename T::value_type>());       \
    }                                                   \
    template <>                                         \
    XSIMD_INLINE float NAME<float>() noexcept           \
    {                                                   \
        return bit_cast<float>((uint32_t)SINGLE);       \
    }                                                   \
    template <>                                         \
    XSIMD_INLINE double NAME<double>() noexcept         \
    {                                                   \
        return bit_cast<double>((uint64_t)DOUBLE);      \
    }

// Under fast-math, GCC might replace signmask (minus zero) by zero
#if defined(__FAST_MATH__) && defined(__GNUC__) && !defined(__clang__)
#pragma GCC push_options
#pragma GCC optimize("signed-zeros")
#endif
        XSIMD_DEFINE_CONSTANT(infinity, (std::numeric_limits<float>::infinity()), (std::numeric_limits<double>::infinity()))
        XSIMD_DEFINE_CONSTANT(invlog_2, 1.442695040888963407359924681001892137426645954152986f, 1.442695040888963407359924681001892137426645954152986)
        XSIMD_DEFINE_CONSTANT_HEX(invlog_2hi, 0x3fb8b000, 0x3ff7154765200000)
        XSIMD_DEFINE_CONSTANT_HEX(invlog_2lo, 0xb9389ad4, 0x3de705fc2eefa200)
        XSIMD_DEFINE_CONSTANT(invlog10_2, 3.32192809488736234787031942949f, 3.32192809488736234787031942949)
        XSIMD_DEFINE_CONSTANT_HEX(invpi, 0x3ea2f983, 0x3fd45f306dc9c883)
        XSIMD_DEFINE_CONSTANT(log_2, 0.6931471805599453094172321214581765680755001343602553f, 0.6931471805599453094172321214581765680755001343602553)
        XSIMD_DEFINE_CONSTANT_HEX(log_2hi, 0x3f318000, 0x3fe62e42fee00000)
        XSIMD_DEFINE_CONSTANT_HEX(log_2lo, 0xb95e8083, 0x3dea39ef35793c76)
        XSIMD_DEFINE_CONSTANT_HEX(log10_2hi, 0x3e9a0000, 0x3fd3440000000000)
        XSIMD_DEFINE_CONSTANT_HEX(log10_2lo, 0x39826a14, 0x3ed3509f79fef312)
        XSIMD_DEFINE_CONSTANT_HEX(logeps, 0xc17f1402, 0xc04205966f2b4f12)
        XSIMD_DEFINE_CONSTANT_HEX(logpi, 0x3f928682, 0x3ff250d048e7a1bd)
        XSIMD_DEFINE_CONSTANT_HEX(logsqrt2pi, 0x3f6b3f8e, 0x3fed67f1c864beb5)
        XSIMD_DEFINE_CONSTANT(maxflint, 16777216.0f, 9007199254740992.0)
        XSIMD_DEFINE_CONSTANT(maxlog, 88.3762626647949f, 709.78271289338400)
        XSIMD_DEFINE_CONSTANT(maxlog2, 127.0f, 1023.)
        XSIMD_DEFINE_CONSTANT(maxlog10, 38.23080825805664f, 308.2547155599167)
        XSIMD_DEFINE_CONSTANT_HEX(mediumpi, 0x43490fdb, 0x412921fb54442d18)
        XSIMD_DEFINE_CONSTANT(minlog, -88.3762626647949f, -708.3964185322641)
        XSIMD_DEFINE_CONSTANT(minlog2, -127.0f, -1023.)
        XSIMD_DEFINE_CONSTANT(minlog10, -37.89999771118164f, -308.2547155599167)
        XSIMD_DEFINE_CONSTANT(minusinfinity, (-infinity<float>()), (-infinity<double>()))
        XSIMD_DEFINE_CONSTANT_HEX(nan, 0xffffffff, 0xffffffffffffffff)
        XSIMD_DEFINE_CONSTANT_HEX(oneosqrteps, 0x453504f3, 0x4190000000000000)
        XSIMD_DEFINE_CONSTANT_HEX(oneotwoeps, 0x4a800000, 0x4320000000000000)
        XSIMD_DEFINE_CONSTANT_HEX(pi, 0x40490fdb, 0x400921fb54442d18)
        XSIMD_DEFINE_CONSTANT_HEX(pio_2lo, 0xb33bbd2e, 0x3c91a62633145c07)
        XSIMD_DEFINE_CONSTANT_HEX(pio_4lo, 0xb2bbbd2e, 0x3c81a62633145c07)
        XSIMD_DEFINE_CONSTANT_HEX(pio2, 0x3fc90fdb, 0x3ff921fb54442d18)
        XSIMD_DEFINE_CONSTANT_HEX(pio2_1, 0x3fc90f80, 0x3ff921fb54400000)
        XSIMD_DEFINE_CONSTANT_HEX(pio2_1t, 0x37354443, 0x3dd0b4611a626331)
        XSIMD_DEFINE_CONSTANT_HEX(pio2_2, 0x37354400, 0x3dd0b4611a600000)
        XSIMD_DEFINE_CONSTANT_HEX(pio2_2t, 0x2e85a308, 0x3ba3198a2e037073)
        XSIMD_DEFINE_CONSTANT_HEX(pio2_3, 0x2e85a300, 0x3ba3198a2e000000)
        XSIMD_DEFINE_CONSTANT_HEX(pio2_3t, 0x248d3132, 0x397b839a252049c1)
        XSIMD_DEFINE_CONSTANT_HEX(pio4, 0x3f490fdb, 0x3fe921fb54442d18)
        XSIMD_DEFINE_CONSTANT_HEX(signmask, 0x80000000, 0x8000000000000000)
        XSIMD_DEFINE_CONSTANT(smallestposval, std::numeric_limits<float>::min(), std::numeric_limits<double>::min())
        XSIMD_DEFINE_CONSTANT_HEX(sqrt_2pi, 0x40206c99, 0x40040d931ff62704)
        XSIMD_DEFINE_CONSTANT_HEX(sqrteps, 0x39b504f3, 0x3e50000000000000)
        XSIMD_DEFINE_CONSTANT_HEX(tanpio8, 0x3ed413cd, 0x3fda827999fcef31)
        XSIMD_DEFINE_CONSTANT_HEX(tan3pio8, 0x401a827a, 0x4003504f333f9de6)
        XSIMD_DEFINE_CONSTANT_HEX(twentypi, 0x427b53d1, 0x404f6a7a2955385e)
        XSIMD_DEFINE_CONSTANT_HEX(twoopi, 0x3f22f983, 0x3fe45f306dc9c883)
        XSIMD_DEFINE_CONSTANT(twotonmb, 8388608.0f, 4503599627370496.0)
        XSIMD_DEFINE_CONSTANT_HEX(twotonmbo3, 0x3ba14518, 0x3ed428a2f98d7286)
#if defined(__FAST_MATH__) && defined(__GNUC__) && !defined(__clang__)
#pragma GCC pop_options
#endif

#undef XSIMD_DEFINE_CONSTANT
#undef XSIMD_DEFINE_CONSTANT_HEX

        template <class T>
        constexpr T allbits() noexcept;

        template <class T>
        constexpr as_integer_t<T> mask1frexp() noexcept;

        template <class T>
        constexpr as_integer_t<T> mask2frexp() noexcept;

        template <class T>
        constexpr as_integer_t<T> maxexponent() noexcept;

        template <class T>
        constexpr as_integer_t<T> maxexponentm1() noexcept;

        template <class T>
        constexpr int32_t nmb() noexcept;

        template <class T>
        constexpr T zero() noexcept;

        template <class T>
        constexpr T minvalue() noexcept;

        template <class T>
        constexpr T maxvalue() noexcept;

        /**************************
         * allbits implementation *
         **************************/

        namespace detail
        {
            template <class T, bool = std::is_integral<T>::value>
            struct allbits_impl
            {
                static constexpr T get_value() noexcept
                {
                    return T(~0);
                }
            };

            template <class T>
            struct allbits_impl<T, false>
            {
                static constexpr T get_value() noexcept
                {
                    return nan<T>();
                }
            };
        }

        template <class T>
        XSIMD_INLINE constexpr T allbits() noexcept
        {
            return T(detail::allbits_impl<typename T::value_type>::get_value());
        }

        /*****************************
         * mask1frexp implementation *
         *****************************/

        template <class T>
        XSIMD_INLINE constexpr as_integer_t<T> mask1frexp() noexcept
        {
            return as_integer_t<T>(mask1frexp<typename T::value_type>());
        }

        template <>
        XSIMD_INLINE constexpr int32_t mask1frexp<float>() noexcept
        {
            return 0x7f800000;
        }

        template <>
        XSIMD_INLINE constexpr int64_t mask1frexp<double>() noexcept
        {
            return 0x7ff0000000000000;
        }

        /*****************************
         * mask2frexp implementation *
         *****************************/

        template <class T>
        XSIMD_INLINE constexpr as_integer_t<T> mask2frexp() noexcept
        {
            return as_integer_t<T>(mask2frexp<typename T::value_type>());
        }

        template <>
        XSIMD_INLINE constexpr int32_t mask2frexp<float>() noexcept
        {
            return 0x3f000000;
        }

        template <>
        XSIMD_INLINE constexpr int64_t mask2frexp<double>() noexcept
        {
            return 0x3fe0000000000000;
        }

        /******************************
         * maxexponent implementation *
         ******************************/

        template <class T>
        XSIMD_INLINE constexpr as_integer_t<T> maxexponent() noexcept
        {
            return as_integer_t<T>(maxexponent<typename T::value_type>());
        }

        template <>
        XSIMD_INLINE constexpr int32_t maxexponent<float>() noexcept
        {
            return 127;
        }

        template <>
        XSIMD_INLINE constexpr int64_t maxexponent<double>() noexcept
        {
            return 1023;
        }

        /******************************
         * maxexponent implementation *
         ******************************/

        template <class T>
        XSIMD_INLINE constexpr as_integer_t<T> maxexponentm1() noexcept
        {
            return as_integer_t<T>(maxexponentm1<typename T::value_type>());
        }

        template <>
        XSIMD_INLINE constexpr int32_t maxexponentm1<float>() noexcept
        {
            return 126;
        }

        template <>
        XSIMD_INLINE constexpr int64_t maxexponentm1<double>() noexcept
        {
            return 1022;
        }

        /**********************
         * nmb implementation *
         **********************/

        template <class T>
        XSIMD_INLINE constexpr int32_t nmb() noexcept
        {
            return nmb<typename T::value_type>();
        }

        template <>
        XSIMD_INLINE constexpr int32_t nmb<float>() noexcept
        {
            return 23;
        }

        template <>
        XSIMD_INLINE constexpr int32_t nmb<double>() noexcept
        {
            return 52;
        }

        /***********************
         * zero implementation *
         ***********************/

        template <class T>
        XSIMD_INLINE constexpr T zero() noexcept
        {
            return T(typename T::value_type(0));
        }

        /***************************
         * minvalue implementation *
         ***************************/

        namespace detail
        {
            template <class T>
            struct minvalue_impl
            {
                static constexpr T get_value() noexcept
                {
                    return std::numeric_limits<typename T::value_type>::min();
                }
            };

            template <class T>
            struct minvalue_common
            {
                static constexpr T get_value() noexcept
                {
                    return std::numeric_limits<T>::min();
                }
            };

            template <>
            struct minvalue_impl<int8_t> : minvalue_common<int8_t>
            {
            };
            template <>
            struct minvalue_impl<uint8_t> : minvalue_common<uint8_t>
            {
            };
            template <>
            struct minvalue_impl<int16_t> : minvalue_common<int16_t>
            {
            };
            template <>
            struct minvalue_impl<uint16_t> : minvalue_common<uint16_t>
            {
            };
            template <>
            struct minvalue_impl<int32_t> : minvalue_common<int32_t>
            {
            };
            template <>
            struct minvalue_impl<uint32_t> : minvalue_common<uint32_t>
            {
            };
            template <>
            struct minvalue_impl<int64_t> : minvalue_common<int64_t>
            {
            };
            template <>
            struct minvalue_impl<uint64_t> : minvalue_common<uint64_t>
            {
            };

            template <>
            struct minvalue_impl<float>
            {
                XSIMD_INLINE static float get_value() noexcept
                {
                    return bit_cast<float>((uint32_t)0xff7fffff);
                }
            };

            template <>
            struct minvalue_impl<double>
            {
                XSIMD_INLINE static double get_value() noexcept
                {
                    return bit_cast<double>((uint64_t)0xffefffffffffffff);
                }
            };
        }

        template <class T>
        constexpr T minvalue() noexcept
        {
            return T(detail::minvalue_impl<typename T::value_type>::get_value());
        }

        /***************************
         * maxvalue implementation *
         ***************************/

        template <class T>
        constexpr T maxvalue() noexcept
        {
            return T(std::numeric_limits<typename T::value_type>::max());
        }
    }

}

#endif
