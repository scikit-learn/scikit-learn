/***************************************************************************
 * Copyright (c) Johan Mabille, Sylvain Corlay, Wolf Vollprecht and         *
 * Martin Renou                                                             *
 * Copyright (c) QuantStack                                                 *
 * Copyright (c) Serge Guelton                                              *
 * Copyright (c) Anutosh Bhat                                               *
 *                                                                          *
 * Distributed under the terms of the BSD 3-Clause License.                 *
 *                                                                          *
 * The full license is in the file LICENSE, distributed with this software. *
 ****************************************************************************/

#ifndef XSIMD_WASM_HPP
#define XSIMD_WASM_HPP

#include <type_traits>

#include "../types/xsimd_wasm_register.hpp"

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
        template <class A, class T>
        XSIMD_INLINE batch<T, A> avg(batch<T, A> const&, batch<T, A> const&, requires_arch<generic>) noexcept;

        // abs
        template <class A, class T, typename std::enable_if<std::is_integral<T>::value && std::is_signed<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> abs(batch<T, A> const& self, requires_arch<wasm>) noexcept
        {
            XSIMD_IF_CONSTEXPR(sizeof(T) == 1)
            {
                return wasm_i8x16_abs(self);
            }
            else XSIMD_IF_CONSTEXPR(sizeof(T) == 2)
            {
                return wasm_i16x8_abs(self);
            }
            else XSIMD_IF_CONSTEXPR(sizeof(T) == 4)
            {
                return wasm_i32x4_abs(self);
            }
            else XSIMD_IF_CONSTEXPR(sizeof(T) == 8)
            {
                return wasm_i64x2_abs(self);
            }
            else
            {
                assert(false && "unsupported arch/op combination");
                return {};
            }
        }

        template <class A>
        XSIMD_INLINE batch<float, A> abs(batch<float, A> const& self, requires_arch<wasm>) noexcept
        {
            return wasm_f32x4_abs(self);
        }

        template <class A>
        XSIMD_INLINE batch<double, A> abs(batch<double, A> const& self, requires_arch<wasm>) noexcept
        {
            return wasm_f64x2_abs(self);
        }

        // add
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> add(batch<T, A> const& self, batch<T, A> const& other, requires_arch<wasm>) noexcept
        {
            XSIMD_IF_CONSTEXPR(sizeof(T) == 1)
            {
                return wasm_i8x16_add(self, other);
            }
            else XSIMD_IF_CONSTEXPR(sizeof(T) == 2)
            {
                return wasm_i16x8_add(self, other);
            }
            else XSIMD_IF_CONSTEXPR(sizeof(T) == 4)
            {
                return wasm_i32x4_add(self, other);
            }
            else XSIMD_IF_CONSTEXPR(sizeof(T) == 8)
            {
                return wasm_i64x2_add(self, other);
            }
            else
            {
                assert(false && "unsupported arch/op combination");
                return {};
            }
        }

        template <class A>
        XSIMD_INLINE batch<float, A> add(batch<float, A> const& self, batch<float, A> const& other, requires_arch<wasm>) noexcept
        {
            return wasm_f32x4_add(self, other);
        }

        template <class A>
        XSIMD_INLINE batch<double, A> add(batch<double, A> const& self, batch<double, A> const& other, requires_arch<wasm>) noexcept
        {
            return wasm_f64x2_add(self, other);
        }

        // avgr
        template <class A, class T, class = typename std::enable_if<std::is_unsigned<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> avgr(batch<T, A> const& self, batch<T, A> const& other, requires_arch<wasm>) noexcept
        {
            XSIMD_IF_CONSTEXPR(sizeof(T) == 1)
            {
                return wasm_u8x16_avgr(self, other);
            }
            else XSIMD_IF_CONSTEXPR(sizeof(T) == 2)
            {
                return wasm_u16x8_avgr(self, other);
            }
            else
            {
                return avgr(self, other, generic {});
            }
        }

        // avg
        template <class A, class T, class = typename std::enable_if<std::is_unsigned<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> avg(batch<T, A> const& self, batch<T, A> const& other, requires_arch<wasm>) noexcept
        {
            XSIMD_IF_CONSTEXPR(sizeof(T) == 1)
            {
                auto adj = ((self ^ other) << 7) >> 7;
                return avgr(self, other, A {}) - adj;
            }
            else XSIMD_IF_CONSTEXPR(sizeof(T) == 2)
            {
                auto adj = ((self ^ other) << 15) >> 15;
                return avgr(self, other, A {}) - adj;
            }
            else
            {
                return avg(self, other, generic {});
            }
        }

        // all
        template <class A>
        XSIMD_INLINE bool all(batch_bool<float, A> const& self, requires_arch<wasm>) noexcept
        {
            return wasm_i32x4_bitmask(self) == 0x0F;
        }
        template <class A>
        XSIMD_INLINE bool all(batch_bool<double, A> const& self, requires_arch<wasm>) noexcept
        {
            return wasm_i64x2_bitmask(self) == 0x03;
        }
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE bool all(batch_bool<T, A> const& self, requires_arch<wasm>) noexcept
        {
            return wasm_i8x16_bitmask(self) == 0xFFFF;
        }

        // any
        template <class A>
        XSIMD_INLINE bool any(batch_bool<float, A> const& self, requires_arch<wasm>) noexcept
        {
            return wasm_i32x4_bitmask(self) != 0;
        }
        template <class A>
        XSIMD_INLINE bool any(batch_bool<double, A> const& self, requires_arch<wasm>) noexcept
        {
            return wasm_i64x2_bitmask(self) != 0;
        }
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE bool any(batch_bool<T, A> const& self, requires_arch<wasm>) noexcept
        {
            return wasm_i8x16_bitmask(self) != 0;
        }

        // batch_bool_cast
        template <class A, class T_out, class T_in>
        XSIMD_INLINE batch_bool<T_out, A> batch_bool_cast(batch_bool<T_in, A> const& self, batch_bool<T_out, A> const&, requires_arch<wasm>) noexcept
        {
            return { bitwise_cast<T_out>(batch<T_in, A>(self.data)).data };
        }

        // bitwise_and
        template <class A, class T>
        XSIMD_INLINE batch<T, A> bitwise_and(batch<T, A> const& self, batch<T, A> const& other, requires_arch<wasm>) noexcept
        {
            return wasm_v128_and(self, other);
        }

        template <class A, class T>
        XSIMD_INLINE batch_bool<T, A> bitwise_and(batch_bool<T, A> const& self, batch_bool<T, A> const& other, requires_arch<wasm>) noexcept
        {
            return wasm_v128_and(self, other);
        }

        // bitwise_andnot
        template <class A, class T>
        XSIMD_INLINE batch<T, A> bitwise_andnot(batch<T, A> const& self, batch<T, A> const& other, requires_arch<wasm>) noexcept
        {
            return wasm_v128_andnot(self, other);
        }

        template <class A, class T>
        XSIMD_INLINE batch_bool<T, A> bitwise_andnot(batch_bool<T, A> const& self, batch_bool<T, A> const& other, requires_arch<wasm>) noexcept
        {
            return wasm_v128_andnot(self, other);
        }

        // bitwise_cast
        template <class A, class T, class Tp>
        XSIMD_INLINE batch<Tp, A> bitwise_cast(batch<T, A> const& self, batch<Tp, A> const&, requires_arch<wasm>) noexcept
        {
            return batch<Tp, A>(self.data);
        }

        // bitwise_or
        template <class A, class T>
        XSIMD_INLINE batch<T, A> bitwise_or(batch<T, A> const& self, batch<T, A> const& other, requires_arch<wasm>) noexcept
        {
            return wasm_v128_or(self, other);
        }

        template <class A, class T>
        XSIMD_INLINE batch_bool<T, A> bitwise_or(batch_bool<T, A> const& self, batch_bool<T, A> const& other, requires_arch<wasm>) noexcept
        {
            return wasm_v128_or(self, other);
        }

        // bitwise_lshift
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> bitwise_lshift(batch<T, A> const& self, int32_t other, requires_arch<wasm>) noexcept
        {
            XSIMD_IF_CONSTEXPR(sizeof(T) == 1)
            {
                return wasm_i8x16_shl(self, other);
            }
            else XSIMD_IF_CONSTEXPR(sizeof(T) == 2)
            {
                return wasm_i16x8_shl(self, other);
            }
            else XSIMD_IF_CONSTEXPR(sizeof(T) == 4)
            {
                return wasm_i32x4_shl(self, other);
            }
            else XSIMD_IF_CONSTEXPR(sizeof(T) == 8)
            {
                return wasm_i64x2_shl(self, other);
            }
            else
            {
                assert(false && "unsupported arch/op combination");
                return {};
            }
        }

        // bitwise_rshift
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> bitwise_rshift(batch<T, A> const& self, int32_t other, requires_arch<wasm>) noexcept
        {
            if (std::is_signed<T>::value)
            {
                XSIMD_IF_CONSTEXPR(sizeof(T) == 1)
                {
                    return wasm_i8x16_shr(self, other);
                }
                else XSIMD_IF_CONSTEXPR(sizeof(T) == 2)
                {
                    return wasm_i16x8_shr(self, other);
                }
                else XSIMD_IF_CONSTEXPR(sizeof(T) == 4)
                {
                    return wasm_i32x4_shr(self, other);
                }
                else XSIMD_IF_CONSTEXPR(sizeof(T) == 8)
                {
                    return wasm_i64x2_shr(self, other);
                }
                else
                {
                    assert(false && "unsupported arch/op combination");
                    return {};
                }
            }
            else
            {
                XSIMD_IF_CONSTEXPR(sizeof(T) == 1)
                {
                    return wasm_u8x16_shr(self, other);
                }
                else XSIMD_IF_CONSTEXPR(sizeof(T) == 2)
                {
                    return wasm_u16x8_shr(self, other);
                }
                else XSIMD_IF_CONSTEXPR(sizeof(T) == 4)
                {
                    return wasm_u32x4_shr(self, other);
                }
                else XSIMD_IF_CONSTEXPR(sizeof(T) == 8)
                {
                    return wasm_u64x2_shr(self, other);
                }
                else
                {
                    assert(false && "unsupported arch/op combination");
                    return {};
                }
            }
        }

        // bitwise_not
        template <class A, class T>
        XSIMD_INLINE batch<T, A> bitwise_not(batch<T, A> const& self, requires_arch<wasm>) noexcept
        {
            return wasm_v128_not(self);
        }

        template <class A, class T>
        XSIMD_INLINE batch_bool<T, A> bitwise_not(batch_bool<T, A> const& self, requires_arch<wasm>) noexcept
        {
            return wasm_v128_not(self);
        }

        // bitwise_xor
        template <class A, class T>
        XSIMD_INLINE batch<T, A> bitwise_xor(batch<T, A> const& self, batch<T, A> const& other, requires_arch<wasm>) noexcept
        {
            return wasm_v128_xor(self, other);
        }

        template <class A, class T>
        XSIMD_INLINE batch_bool<T, A> bitwise_xor(batch_bool<T, A> const& self, batch_bool<T, A> const& other, requires_arch<wasm>) noexcept
        {
            return wasm_v128_xor(self, other);
        }

        // broadcast
        template <class A>
        batch<float, A> XSIMD_INLINE broadcast(float val, requires_arch<wasm>) noexcept
        {
            return wasm_f32x4_splat(val);
        }
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> broadcast(T val, requires_arch<wasm>) noexcept
        {
            XSIMD_IF_CONSTEXPR(sizeof(T) == 1)
            {
                return wasm_i8x16_splat(val);
            }
            else XSIMD_IF_CONSTEXPR(sizeof(T) == 2)
            {
                return wasm_i16x8_splat(val);
            }
            else XSIMD_IF_CONSTEXPR(sizeof(T) == 4)
            {
                return wasm_i32x4_splat(val);
            }
            else XSIMD_IF_CONSTEXPR(sizeof(T) == 8)
            {
                return wasm_i64x2_splat(val);
            }
            else
            {
                assert(false && "unsupported arch/op combination");
                return {};
            }
        }
        template <class A>
        XSIMD_INLINE batch<double, A> broadcast(double val, requires_arch<wasm>) noexcept
        {
            return wasm_f64x2_splat(val);
        }

        // ceil
        template <class A>
        XSIMD_INLINE batch<float, A> ceil(batch<float, A> const& self, requires_arch<wasm>) noexcept
        {
            return wasm_f32x4_ceil(self);
        }
        template <class A>
        XSIMD_INLINE batch<double, A> ceil(batch<double, A> const& self, requires_arch<wasm>) noexcept
        {
            return wasm_f64x2_ceil(self);
        }

        // div
        template <class A>
        XSIMD_INLINE batch<float, A> div(batch<float, A> const& self, batch<float, A> const& other, requires_arch<wasm>) noexcept
        {
            return wasm_f32x4_div(self, other);
        }
        template <class A>
        XSIMD_INLINE batch<double, A> div(batch<double, A> const& self, batch<double, A> const& other, requires_arch<wasm>) noexcept
        {
            return wasm_f64x2_div(self, other);
        }

        // eq
        template <class A>
        XSIMD_INLINE batch_bool<float, A> eq(batch<float, A> const& self, batch<float, A> const& other, requires_arch<wasm>) noexcept
        {
            return wasm_f32x4_eq(self, other);
        }
        template <class A>
        XSIMD_INLINE batch_bool<float, A> eq(batch_bool<float, A> const& self, batch_bool<float, A> const& other, requires_arch<wasm>) noexcept
        {
            return wasm_i32x4_eq(self, other);
        }
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch_bool<T, A> eq(batch<T, A> const& self, batch<T, A> const& other, requires_arch<wasm>) noexcept
        {
            XSIMD_IF_CONSTEXPR(sizeof(T) == 1)
            {
                return wasm_i8x16_eq(self, other);
            }
            else XSIMD_IF_CONSTEXPR(sizeof(T) == 2)
            {
                return wasm_i16x8_eq(self, other);
            }
            else XSIMD_IF_CONSTEXPR(sizeof(T) == 4)
            {
                return wasm_i32x4_eq(self, other);
            }
            else XSIMD_IF_CONSTEXPR(sizeof(T) == 8)
            {
                return wasm_i64x2_eq(self, other);
            }
            else
            {
                assert(false && "unsupported arch/op combination");
                return {};
            }
        }
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch_bool<T, A> eq(batch_bool<T, A> const& self, batch_bool<T, A> const& other, requires_arch<wasm>) noexcept
        {
            XSIMD_IF_CONSTEXPR(sizeof(T) == 1)
            {
                return wasm_i8x16_eq(self, other);
            }
            else XSIMD_IF_CONSTEXPR(sizeof(T) == 2)
            {
                return wasm_i16x8_eq(self, other);
            }
            else XSIMD_IF_CONSTEXPR(sizeof(T) == 4)
            {
                return wasm_i32x4_eq(self, other);
            }
            else XSIMD_IF_CONSTEXPR(sizeof(T) == 8)
            {
                return wasm_i64x2_eq(self, other);
            }
            else
            {
                assert(false && "unsupported arch/op combination");
                return {};
            }
        }
        template <class A>
        XSIMD_INLINE batch_bool<double, A> eq(batch<double, A> const& self, batch<double, A> const& other, requires_arch<wasm>) noexcept
        {
            return wasm_f64x2_eq(self, other);
        }
        template <class A>
        XSIMD_INLINE batch_bool<double, A> eq(batch_bool<double, A> const& self, batch_bool<double, A> const& other, requires_arch<wasm>) noexcept
        {
            return wasm_i64x2_eq(self, other);
        }

        // fast_cast
        namespace detail
        {
            template <class A>
            XSIMD_INLINE batch<float, A> fast_cast(batch<int32_t, A> const& self, batch<float, A> const&, requires_arch<wasm>) noexcept
            {
                return wasm_f32x4_convert_i32x4(self);
            }

            template <class A>
            XSIMD_INLINE batch<double, A> fast_cast(batch<uint64_t, A> const& x, batch<double, A> const&, requires_arch<wasm>) noexcept
            {
                // from https://stackoverflow.com/questions/41144668/how-to-efficiently-perform-double-int64-conversions-with-sse-avx
                // adapted to wasm
                v128_t xH = wasm_u64x2_shr(x, 32);
                xH = wasm_v128_or(xH, wasm_f64x2_splat(19342813113834066795298816.)); //  2^84
                v128_t mask = wasm_i16x8_make(0xFFFF, 0xFFFF, 0x0000, 0x0000, 0xFFFF, 0xFFFF, 0x0000, 0x0000);
                v128_t xL = wasm_v128_or(wasm_v128_and(mask, x), wasm_v128_andnot(wasm_f64x2_splat(0x0010000000000000), mask)); //  2^52
                v128_t f = wasm_f64x2_sub(xH, wasm_f64x2_splat(19342813118337666422669312.)); //  2^84 + 2^52
                return wasm_f64x2_add(f, xL);
            }

            template <class A>
            XSIMD_INLINE batch<double, A> fast_cast(batch<int64_t, A> const& x, batch<double, A> const&, requires_arch<wasm>) noexcept
            {
                // from https://stackoverflow.com/questions/41144668/how-to-efficiently-perform-double-int64-conversions-with-sse-avx
                // adapted to wasm
                v128_t xH = wasm_i32x4_shr(x, 16);
                xH = wasm_v128_and(xH, wasm_i16x8_make(0x0000, 0x0000, 0xFFFF, 0xFFFF, 0x0000, 0x0000, 0xFFFF, 0xFFFF));
                xH = wasm_i64x2_add(xH, wasm_f64x2_splat(442721857769029238784.)); //  3*2^67
                v128_t mask = wasm_i16x8_make(0xFFFF, 0xFFFF, 0xFFFF, 0x0000, 0xFFFF, 0xFFFF, 0xFFFF, 0x0000);
                v128_t xL = wasm_v128_or(wasm_v128_and(mask, x), wasm_v128_andnot(wasm_f64x2_splat(0x0010000000000000), mask)); //  2^52
                v128_t f = wasm_f64x2_sub(xH, wasm_f64x2_splat(442726361368656609280.)); //  3*2^67 + 2^52
                return wasm_f64x2_add(f, xL);
            }

            template <class A>
            XSIMD_INLINE batch<int32_t, A> fast_cast(batch<float, A> const& self, batch<int32_t, A> const&, requires_arch<wasm>) noexcept
            {
                return wasm_i32x4_make(
                    static_cast<int32_t>(wasm_f32x4_extract_lane(self, 0)),
                    static_cast<int32_t>(wasm_f32x4_extract_lane(self, 1)),
                    static_cast<int32_t>(wasm_f32x4_extract_lane(self, 2)),
                    static_cast<int32_t>(wasm_f32x4_extract_lane(self, 3)));
            }
        }

        // floor
        template <class A>
        XSIMD_INLINE batch<float, A> floor(batch<float, A> const& self, requires_arch<wasm>) noexcept
        {
            return wasm_f32x4_floor(self);
        }

        template <class A>
        XSIMD_INLINE batch<double, A> floor(batch<double, A> const& self, requires_arch<wasm>) noexcept
        {
            return wasm_f64x2_floor(self);
        }

        // from_mask
        template <class A>
        XSIMD_INLINE batch_bool<float, A> from_mask(batch_bool<float, A> const&, uint64_t mask, requires_arch<wasm>) noexcept
        {
            alignas(A::alignment()) static const uint32_t lut[][4] = {
                { 0x00000000, 0x00000000, 0x00000000, 0x00000000 },
                { 0xFFFFFFFF, 0x00000000, 0x00000000, 0x00000000 },
                { 0x00000000, 0xFFFFFFFF, 0x00000000, 0x00000000 },
                { 0xFFFFFFFF, 0xFFFFFFFF, 0x00000000, 0x00000000 },
                { 0x00000000, 0x00000000, 0xFFFFFFFF, 0x00000000 },
                { 0xFFFFFFFF, 0x00000000, 0xFFFFFFFF, 0x00000000 },
                { 0x00000000, 0xFFFFFFFF, 0xFFFFFFFF, 0x00000000 },
                { 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0x00000000 },
                { 0x00000000, 0x00000000, 0x00000000, 0xFFFFFFFF },
                { 0xFFFFFFFF, 0x00000000, 0x00000000, 0xFFFFFFFF },
                { 0x00000000, 0xFFFFFFFF, 0x00000000, 0xFFFFFFFF },
                { 0xFFFFFFFF, 0xFFFFFFFF, 0x00000000, 0xFFFFFFFF },
                { 0x00000000, 0x00000000, 0xFFFFFFFF, 0xFFFFFFFF },
                { 0xFFFFFFFF, 0x00000000, 0xFFFFFFFF, 0xFFFFFFFF },
                { 0x00000000, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF },
                { 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF },
            };
            assert(!(mask & ~0xFul) && "inbound mask");
            return wasm_v128_load((const v128_t*)lut[mask]);
        }
        template <class A>
        XSIMD_INLINE batch_bool<double, A> from_mask(batch_bool<double, A> const&, uint64_t mask, requires_arch<wasm>) noexcept
        {
            alignas(A::alignment()) static const uint64_t lut[][4] = {
                { 0x0000000000000000ul, 0x0000000000000000ul },
                { 0xFFFFFFFFFFFFFFFFul, 0x0000000000000000ul },
                { 0x0000000000000000ul, 0xFFFFFFFFFFFFFFFFul },
                { 0xFFFFFFFFFFFFFFFFul, 0xFFFFFFFFFFFFFFFFul },
            };
            assert(!(mask & ~0x3ul) && "inbound mask");
            return wasm_v128_load((const v128_t*)lut[mask]);
        }
        template <class T, class A, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch_bool<T, A> from_mask(batch_bool<T, A> const&, uint64_t mask, requires_arch<wasm>) noexcept
        {
            alignas(A::alignment()) static const uint64_t lut64[] = {
                0x0000000000000000,
                0x000000000000FFFF,
                0x00000000FFFF0000,
                0x00000000FFFFFFFF,
                0x0000FFFF00000000,
                0x0000FFFF0000FFFF,
                0x0000FFFFFFFF0000,
                0x0000FFFFFFFFFFFF,
                0xFFFF000000000000,
                0xFFFF00000000FFFF,
                0xFFFF0000FFFF0000,
                0xFFFF0000FFFFFFFF,
                0xFFFFFFFF00000000,
                0xFFFFFFFF0000FFFF,
                0xFFFFFFFFFFFF0000,
                0xFFFFFFFFFFFFFFFF,
            };
            alignas(A::alignment()) static const uint32_t lut32[] = {
                0x00000000,
                0x000000FF,
                0x0000FF00,
                0x0000FFFF,
                0x00FF0000,
                0x00FF00FF,
                0x00FFFF00,
                0x00FFFFFF,
                0xFF000000,
                0xFF0000FF,
                0xFF00FF00,
                0xFF00FFFF,
                0xFFFF0000,
                0xFFFF00FF,
                0xFFFFFF00,
                0xFFFFFFFF,
            };
            alignas(A::alignment()) static const uint32_t lut16[][4] = {
                { 0x00000000, 0x00000000, 0x00000000, 0x00000000 },
                { 0xFFFFFFFF, 0x00000000, 0x00000000, 0x00000000 },
                { 0x00000000, 0xFFFFFFFF, 0x00000000, 0x00000000 },
                { 0xFFFFFFFF, 0xFFFFFFFF, 0x00000000, 0x00000000 },
                { 0x00000000, 0x00000000, 0xFFFFFFFF, 0x00000000 },
                { 0xFFFFFFFF, 0x00000000, 0xFFFFFFFF, 0x00000000 },
                { 0x00000000, 0xFFFFFFFF, 0xFFFFFFFF, 0x00000000 },
                { 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0x00000000 },
                { 0x00000000, 0x00000000, 0x00000000, 0xFFFFFFFF },
                { 0xFFFFFFFF, 0x00000000, 0x00000000, 0xFFFFFFFF },
                { 0x00000000, 0xFFFFFFFF, 0x00000000, 0xFFFFFFFF },
                { 0xFFFFFFFF, 0xFFFFFFFF, 0x00000000, 0xFFFFFFFF },
                { 0x00000000, 0x00000000, 0xFFFFFFFF, 0xFFFFFFFF },
                { 0xFFFFFFFF, 0x00000000, 0xFFFFFFFF, 0xFFFFFFFF },
                { 0x00000000, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF },
                { 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF },
            };
            alignas(A::alignment()) static const uint64_t lut8[][4] = {
                { 0x0000000000000000ul, 0x0000000000000000ul },
                { 0xFFFFFFFFFFFFFFFFul, 0x0000000000000000ul },
                { 0x0000000000000000ul, 0xFFFFFFFFFFFFFFFFul },
                { 0xFFFFFFFFFFFFFFFFul, 0xFFFFFFFFFFFFFFFFul },
            };
            XSIMD_IF_CONSTEXPR(sizeof(T) == 1)
            {
                assert(!(mask & ~0xFFFF) && "inbound mask");
                return wasm_i32x4_make(lut32[mask & 0xF], lut32[(mask >> 4) & 0xF], lut32[(mask >> 8) & 0xF], lut32[mask >> 12]);
            }
            else XSIMD_IF_CONSTEXPR(sizeof(T) == 2)
            {
                assert(!(mask & ~0xFF) && "inbound mask");
                return wasm_i64x2_make(lut64[mask & 0xF], lut64[mask >> 4]);
            }
            else XSIMD_IF_CONSTEXPR(sizeof(T) == 4)
            {
                assert(!(mask & ~0xFul) && "inbound mask");
                return wasm_v128_load((const v128_t*)lut16[mask]);
            }
            else XSIMD_IF_CONSTEXPR(sizeof(T) == 8)
            {
                assert(!(mask & ~0x3ul) && "inbound mask");
                return wasm_v128_load((const v128_t*)lut8[mask]);
            }
        }

        // ge
        template <class A>
        XSIMD_INLINE batch_bool<float, A> ge(batch<float, A> const& self, batch<float, A> const& other, requires_arch<wasm>) noexcept
        {
            return wasm_f32x4_ge(self, other);
        }
        template <class A>
        XSIMD_INLINE batch_bool<double, A> ge(batch<double, A> const& self, batch<double, A> const& other, requires_arch<wasm>) noexcept
        {
            return wasm_f64x2_ge(self, other);
        }

        // gt
        template <class A>
        XSIMD_INLINE batch_bool<float, A> gt(batch<float, A> const& self, batch<float, A> const& other, requires_arch<wasm>) noexcept
        {
            return wasm_f32x4_gt(self, other);
        }
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch_bool<T, A> gt(batch<T, A> const& self, batch<T, A> const& other, requires_arch<wasm>) noexcept
        {
            if (std::is_signed<T>::value)
            {
                XSIMD_IF_CONSTEXPR(sizeof(T) == 1)
                {
                    return wasm_i8x16_gt(self, other);
                }
                else XSIMD_IF_CONSTEXPR(sizeof(T) == 2)
                {
                    return wasm_i16x8_gt(self, other);
                }
                else XSIMD_IF_CONSTEXPR(sizeof(T) == 4)
                {
                    return wasm_i32x4_gt(self, other);
                }
                else XSIMD_IF_CONSTEXPR(sizeof(T) == 8)
                {
                    return wasm_i64x2_gt(self, other);
                }
                else
                {
                    assert(false && "unsupported arch/op combination");
                    return {};
                }
            }
            else
            {
                XSIMD_IF_CONSTEXPR(sizeof(T) == 1)
                {
                    return wasm_u8x16_gt(self, other);
                }
                else XSIMD_IF_CONSTEXPR(sizeof(T) == 2)
                {
                    return wasm_u16x8_gt(self, other);
                }
                else XSIMD_IF_CONSTEXPR(sizeof(T) == 4)
                {
                    return wasm_u32x4_gt(self, other);
                }
                else
                {
                    return gt(self, other, generic {});
                }
            }
        }

        template <class A>
        XSIMD_INLINE batch_bool<double, A> gt(batch<double, A> const& self, batch<double, A> const& other, requires_arch<wasm>) noexcept
        {
            return wasm_f64x2_gt(self, other);
        }

        // haddp
        template <class A>
        XSIMD_INLINE batch<float, A> haddp(batch<float, A> const* row, requires_arch<wasm>) noexcept
        {
            v128_t tmp0 = wasm_i32x4_shuffle(row[0], row[1], 0, 4, 1, 5);
            v128_t tmp1 = wasm_i32x4_shuffle(row[0], row[1], 2, 6, 3, 7);
            v128_t tmp2 = wasm_i32x4_shuffle(row[2], row[3], 2, 6, 3, 7);
            tmp0 = wasm_f32x4_add(tmp0, tmp1);
            tmp1 = wasm_i32x4_shuffle(row[2], row[3], 0, 4, 1, 5);
            tmp1 = wasm_f32x4_add(tmp1, tmp2);
            tmp2 = wasm_i32x4_shuffle(tmp1, tmp0, 6, 7, 2, 3);
            tmp0 = wasm_i32x4_shuffle(tmp0, tmp1, 0, 1, 4, 5);
            return wasm_f32x4_add(tmp0, tmp2);
        }
        template <class A>
        XSIMD_INLINE batch<double, A> haddp(batch<double, A> const* row, requires_arch<wasm>) noexcept
        {
            return wasm_f64x2_add(wasm_i64x2_shuffle(row[0], row[1], 0, 2),
                                  wasm_i64x2_shuffle(row[0], row[1], 1, 3));
        }

        // insert
        template <class A, size_t I>
        XSIMD_INLINE batch<float, A> insert(batch<float, A> const& self, float val, index<I> pos, requires_arch<wasm>) noexcept
        {
            return wasm_f32x4_replace_lane(self, pos, val);
        }
        template <class A, class T, size_t I, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> insert(batch<T, A> const& self, T val, index<I> pos, requires_arch<wasm>) noexcept
        {
            if (std::is_signed<T>::value)
            {
                XSIMD_IF_CONSTEXPR(sizeof(T) == 1)
                {
                    return wasm_i8x16_replace_lane(self, pos, val);
                }
                else XSIMD_IF_CONSTEXPR(sizeof(T) == 2)
                {
                    return wasm_i16x8_replace_lane(self, pos, val);
                }
                else XSIMD_IF_CONSTEXPR(sizeof(T) == 4)
                {
                    return wasm_i32x4_replace_lane(self, pos, val);
                }
                else XSIMD_IF_CONSTEXPR(sizeof(T) == 8)
                {
                    return wasm_i64x2_replace_lane(self, pos, val);
                }
                else
                {
                    assert(false && "unsupported arch/op combination");
                    return {};
                }
            }
            else
            {
                XSIMD_IF_CONSTEXPR(sizeof(T) == 1)
                {
                    return wasm_u8x16_replace_lane(self, pos, val);
                }
                else XSIMD_IF_CONSTEXPR(sizeof(T) == 2)
                {
                    return wasm_u16x8_replace_lane(self, pos, val);
                }
                else XSIMD_IF_CONSTEXPR(sizeof(T) == 4)
                {
                    return wasm_u32x4_replace_lane(self, pos, val);
                }
                else XSIMD_IF_CONSTEXPR(sizeof(T) == 8)
                {
                    return wasm_u64x2_replace_lane(self, pos, val);
                }
                else
                {
                    assert(false && "unsupported arch/op combination");
                    return {};
                }
            }
        }

        template <class A, size_t I>
        XSIMD_INLINE batch<double, A> insert(batch<double, A> const& self, double val, index<I> pos, requires_arch<wasm>) noexcept
        {
            return wasm_f64x2_replace_lane(self, pos, val);
        }

        // isnan
        template <class A>
        XSIMD_INLINE batch_bool<float, A> isnan(batch<float, A> const& self, requires_arch<wasm>) noexcept
        {
            return wasm_v128_or(wasm_f32x4_ne(self, self), wasm_f32x4_ne(self, self));
        }
        template <class A>
        XSIMD_INLINE batch_bool<double, A> isnan(batch<double, A> const& self, requires_arch<wasm>) noexcept
        {
            return wasm_v128_or(wasm_f64x2_ne(self, self), wasm_f64x2_ne(self, self));
        }

        // le
        template <class A>
        XSIMD_INLINE batch_bool<float, A> le(batch<float, A> const& self, batch<float, A> const& other, requires_arch<wasm>) noexcept
        {
            return wasm_f32x4_le(self, other);
        }
        template <class A>
        XSIMD_INLINE batch_bool<double, A> le(batch<double, A> const& self, batch<double, A> const& other, requires_arch<wasm>) noexcept
        {
            return wasm_f64x2_le(self, other);
        }

        // load_aligned
        template <class A>
        XSIMD_INLINE batch<float, A> load_aligned(float const* mem, convert<float>, requires_arch<wasm>) noexcept
        {
            return wasm_v128_load(mem);
        }
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> load_aligned(T const* mem, convert<T>, requires_arch<wasm>) noexcept
        {
            return wasm_v128_load((v128_t const*)mem);
        }
        template <class A>
        XSIMD_INLINE batch<double, A> load_aligned(double const* mem, convert<double>, requires_arch<wasm>) noexcept
        {
            return wasm_v128_load(mem);
        }

        // load_complex
        namespace detail
        {
            template <class A>
            XSIMD_INLINE batch<std::complex<float>, A> load_complex(batch<float, A> const& hi, batch<float, A> const& lo, requires_arch<wasm>) noexcept
            {
                return { wasm_i32x4_shuffle(hi, lo, 0, 2, 4, 6), wasm_i32x4_shuffle(hi, lo, 1, 3, 5, 7) };
            }
            template <class A>
            XSIMD_INLINE batch<std::complex<double>, A> load_complex(batch<double, A> const& hi, batch<double, A> const& lo, requires_arch<wasm>) noexcept
            {
                return { wasm_i64x2_shuffle(hi, lo, 0, 2), wasm_i64x2_shuffle(hi, lo, 1, 3) };
            }
        }

        // load_unaligned
        template <class A>
        XSIMD_INLINE batch<float, A> load_unaligned(float const* mem, convert<float>, requires_arch<wasm>) noexcept
        {
            return wasm_v128_load(mem);
        }
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> load_unaligned(T const* mem, convert<T>, requires_arch<wasm>) noexcept
        {
            return wasm_v128_load((v128_t const*)mem);
        }
        template <class A>
        XSIMD_INLINE batch<double, A> load_unaligned(double const* mem, convert<double>, requires_arch<wasm>) noexcept
        {
            return wasm_v128_load(mem);
        }

        // lt
        template <class A>
        XSIMD_INLINE batch_bool<float, A> lt(batch<float, A> const& self, batch<float, A> const& other, requires_arch<wasm>) noexcept
        {
            return wasm_f32x4_lt(self, other);
        }
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch_bool<T, A> lt(batch<T, A> const& self, batch<T, A> const& other, requires_arch<wasm>) noexcept
        {
            if (std::is_signed<T>::value)
            {
                XSIMD_IF_CONSTEXPR(sizeof(T) == 1)
                {
                    return wasm_i8x16_lt(self, other);
                }
                else XSIMD_IF_CONSTEXPR(sizeof(T) == 2)
                {
                    return wasm_i16x8_lt(self, other);
                }
                else XSIMD_IF_CONSTEXPR(sizeof(T) == 4)
                {
                    return wasm_i32x4_lt(self, other);
                }
                else XSIMD_IF_CONSTEXPR(sizeof(T) == 8)
                {
                    return wasm_i64x2_lt(self, other);
                }
                else
                {
                    assert(false && "unsupported arch/op combination");
                    return {};
                }
            }
            else
            {
                XSIMD_IF_CONSTEXPR(sizeof(T) == 1)
                {
                    return wasm_u8x16_lt(self, other);
                }
                else XSIMD_IF_CONSTEXPR(sizeof(T) == 2)
                {
                    return wasm_u16x8_lt(self, other);
                }
                else XSIMD_IF_CONSTEXPR(sizeof(T) == 4)
                {
                    return wasm_u32x4_lt(self, other);
                }
                else XSIMD_IF_CONSTEXPR(sizeof(T) == 8)
                {
                    auto xself = wasm_v128_xor(self, wasm_i64x2_splat(std::numeric_limits<int64_t>::lowest()));
                    auto xother = wasm_v128_xor(other, wasm_i64x2_splat(std::numeric_limits<int64_t>::lowest()));
                    v128_t tmp1 = wasm_i64x2_sub(xself, xother);
                    v128_t tmp2 = wasm_v128_xor(xself, xother);
                    v128_t tmp3 = wasm_v128_andnot(xself, xother);
                    v128_t tmp4 = wasm_v128_andnot(tmp1, tmp2);
                    v128_t tmp5 = wasm_v128_or(tmp3, tmp4);
                    v128_t tmp6 = wasm_i32x4_shr(tmp5, 31);
                    return wasm_i32x4_shuffle(tmp6, wasm_i32x4_splat(0), 1, 1, 3, 3);
                }
                else
                {
                    assert(false && "unsupported arch/op combination");
                    return {};
                }
            }
        }

        template <class A>
        XSIMD_INLINE batch_bool<double, A> lt(batch<double, A> const& self, batch<double, A> const& other, requires_arch<wasm>) noexcept
        {
            return wasm_f64x2_lt(self, other);
        }

        // mask
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE uint64_t mask(batch_bool<T, A> const& self, requires_arch<wasm>) noexcept
        {
            XSIMD_IF_CONSTEXPR(sizeof(T) == 1)
            {
                return wasm_i8x16_bitmask(self);
            }
            else XSIMD_IF_CONSTEXPR(sizeof(T) == 2)
            {
                return wasm_i16x8_bitmask(self);
            }
            else XSIMD_IF_CONSTEXPR(sizeof(T) == 4)
            {
                return wasm_i32x4_bitmask(self);
            }
            else XSIMD_IF_CONSTEXPR(sizeof(T) == 8)
            {
                return wasm_i64x2_bitmask(self);
            }
            else
            {
                assert(false && "unsupported arch/op combination");
                return {};
            }
        }
        template <class A>
        XSIMD_INLINE uint64_t mask(batch_bool<float, A> const& self, requires_arch<wasm>) noexcept
        {
            return wasm_i32x4_bitmask(self);
        }

        template <class A>
        XSIMD_INLINE uint64_t mask(batch_bool<double, A> const& self, requires_arch<wasm>) noexcept
        {
            return wasm_i64x2_bitmask(self);
        }

        // max
        template <class A>
        XSIMD_INLINE batch<float, A> max(batch<float, A> const& self, batch<float, A> const& other, requires_arch<wasm>) noexcept
        {
            return wasm_f32x4_pmax(self, other);
        }
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> max(batch<T, A> const& self, batch<T, A> const& other, requires_arch<wasm>) noexcept
        {
            return select(self > other, self, other);
        }
        template <class A>
        XSIMD_INLINE batch<double, A> max(batch<double, A> const& self, batch<double, A> const& other, requires_arch<wasm>) noexcept
        {
            return wasm_f64x2_pmax(self, other);
        }

        // min
        template <class A>
        XSIMD_INLINE batch<float, A> min(batch<float, A> const& self, batch<float, A> const& other, requires_arch<wasm>) noexcept
        {
            return wasm_f32x4_pmin(self, other);
        }
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> min(batch<T, A> const& self, batch<T, A> const& other, requires_arch<wasm>) noexcept
        {
            return select(self <= other, self, other);
        }
        template <class A>
        XSIMD_INLINE batch<double, A> min(batch<double, A> const& self, batch<double, A> const& other, requires_arch<wasm>) noexcept
        {
            return wasm_f64x2_pmin(self, other);
        }

        // mul
        template <class A>
        XSIMD_INLINE batch<float, A> mul(batch<float, A> const& self, batch<float, A> const& other, requires_arch<wasm>) noexcept
        {
            return wasm_f32x4_mul(self, other);
        }
        template <class A>
        XSIMD_INLINE batch<double, A> mul(batch<double, A> const& self, batch<double, A> const& other, requires_arch<wasm>) noexcept
        {
            return wasm_f64x2_mul(self, other);
        }

        // neg
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> neg(batch<T, A> const& self, requires_arch<wasm>) noexcept
        {
            XSIMD_IF_CONSTEXPR(sizeof(T) == 1)
            {
                return wasm_i8x16_neg(self);
            }
            else XSIMD_IF_CONSTEXPR(sizeof(T) == 2)
            {
                return wasm_i16x8_neg(self);
            }
            else XSIMD_IF_CONSTEXPR(sizeof(T) == 4)
            {
                return wasm_i32x4_neg(self);
            }
            else XSIMD_IF_CONSTEXPR(sizeof(T) == 8)
            {
                return wasm_i64x2_neg(self);
            }
            else
            {
                assert(false && "unsupported arch/op combination");
                return {};
            }
        }

        template <class A>
        XSIMD_INLINE batch<float, A> neg(batch<float, A> const& self, requires_arch<wasm>) noexcept
        {
            return wasm_f32x4_neg(self);
        }

        template <class A>
        XSIMD_INLINE batch<double, A> neg(batch<double, A> const& self, requires_arch<wasm>) noexcept
        {
            return wasm_f64x2_neg(self);
        }

        // neq
        template <class A>
        XSIMD_INLINE batch_bool<float, A> neq(batch<float, A> const& self, batch<float, A> const& other, requires_arch<wasm>) noexcept
        {
            return wasm_f32x4_ne(self, other);
        }
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch_bool<T, A> neq(batch<T, A> const& self, batch<T, A> const& other, requires_arch<wasm>) noexcept
        {
            return ~(self == other);
        }
        template <class A>
        XSIMD_INLINE batch_bool<float, A> neq(batch_bool<float, A> const& self, batch_bool<float, A> const& other, requires_arch<wasm>) noexcept
        {
            return wasm_f32x4_ne(self, other);
        }
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch_bool<T, A> neq(batch_bool<T, A> const& self, batch_bool<T, A> const& other, requires_arch<wasm>) noexcept
        {
            return ~(self == other);
        }

        template <class A>
        XSIMD_INLINE batch_bool<double, A> neq(batch<double, A> const& self, batch<double, A> const& other, requires_arch<wasm>) noexcept
        {
            return wasm_f64x2_ne(self, other);
        }
        template <class A>
        XSIMD_INLINE batch_bool<double, A> neq(batch_bool<double, A> const& self, batch_bool<double, A> const& other, requires_arch<wasm>) noexcept
        {
            return wasm_f64x2_ne(self, other);
        }

        // reciprocal
        template <class A>
        XSIMD_INLINE batch<float, A> reciprocal(batch<float, A> const& self, requires_arch<wasm>) noexcept
        {
            v128_t one = wasm_f32x4_splat(1.0f);
            return wasm_f32x4_div(one, self);
        }
        template <class A>
        XSIMD_INLINE batch<double, A> reciprocal(batch<double, A> const& self, requires_arch<wasm>) noexcept
        {
            v128_t one = wasm_f64x2_splat(1.0);
            return wasm_f64x2_div(one, self);
        }

        // reduce_add
        template <class A>
        XSIMD_INLINE float reduce_add(batch<float, A> const& self, requires_arch<wasm>) noexcept
        {
            v128_t tmp0 = wasm_f32x4_add(self, wasm_i32x4_shuffle(self, self, 6, 7, 2, 3));
            v128_t tmp1 = wasm_i32x4_shuffle(tmp0, tmp0, 1, 0, 4, 4);
            v128_t tmp2 = wasm_f32x4_add(tmp0, tmp1);
            v128_t tmp3 = wasm_i32x4_shuffle(tmp0, tmp2, 4, 1, 2, 3);
            return wasm_f32x4_extract_lane(tmp3, 0);
        }
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE T reduce_add(batch<T, A> const& self, requires_arch<wasm>) noexcept
        {
            XSIMD_IF_CONSTEXPR(sizeof(T) == 4)
            {
                v128_t tmp0 = wasm_i32x4_shuffle(self, wasm_i32x4_splat(0), 2, 3, 0, 0);
                v128_t tmp1 = wasm_i32x4_add(self, tmp0);
                v128_t tmp2 = wasm_i32x4_shuffle(tmp1, wasm_i32x4_splat(0), 1, 0, 0, 0);
                v128_t tmp3 = wasm_i32x4_add(tmp1, tmp2);
                return wasm_i32x4_extract_lane(tmp3, 0);
            }
            else XSIMD_IF_CONSTEXPR(sizeof(T) == 8)
            {
                v128_t tmp0 = wasm_i32x4_shuffle(self, wasm_i32x4_splat(0), 2, 3, 0, 0);
                v128_t tmp1 = wasm_i64x2_add(self, tmp0);
                return wasm_i64x2_extract_lane(tmp1, 0);
            }
            else
            {
                return hadd(self, generic {});
            }
        }
        template <class A>
        XSIMD_INLINE double reduce_add(batch<double, A> const& self, requires_arch<wasm>) noexcept
        {
            v128_t tmp0 = wasm_i64x2_shuffle(self, self, 1, 3);
            v128_t tmp1 = wasm_f64x2_add(self, tmp0);
            v128_t tmp2 = wasm_i64x2_shuffle(tmp0, tmp1, 2, 1);
            return wasm_f64x2_extract_lane(tmp2, 0);
        }

        // rsqrt
        template <class A>
        XSIMD_INLINE batch<float, A> rsqrt(batch<float, A> const& self, requires_arch<wasm>) noexcept
        {
            v128_t one = wasm_f32x4_splat(1.0f);
            return wasm_f32x4_div(one, wasm_f32x4_sqrt(self));
        }
        template <class A>
        XSIMD_INLINE batch<double, A> rsqrt(batch<double, A> const& self, requires_arch<wasm>) noexcept
        {
            v128_t one = wasm_f64x2_splat(1.0);
            return wasm_f64x2_div(one, wasm_f64x2_sqrt(self));
        }

        // slide_left
        template <size_t N, class A, class T>
        XSIMD_INLINE batch<T, A> slide_left(batch<T, A> const& x, requires_arch<wasm>) noexcept
        {
            return wasm_i8x16_shuffle(
                wasm_i64x2_const(0, 0), x, ((N) & 0xF0) ? 0 : 16 - ((N) & 0xF),
                ((N) & 0xF0) ? 0 : 17 - ((N) & 0xF), ((N) & 0xF0) ? 0 : 18 - ((N) & 0xF),
                ((N) & 0xF0) ? 0 : 19 - ((N) & 0xF), ((N) & 0xF0) ? 0 : 20 - ((N) & 0xF),
                ((N) & 0xF0) ? 0 : 21 - ((N) & 0xF), ((N) & 0xF0) ? 0 : 22 - ((N) & 0xF),
                ((N) & 0xF0) ? 0 : 23 - ((N) & 0xF), ((N) & 0xF0) ? 0 : 24 - ((N) & 0xF),
                ((N) & 0xF0) ? 0 : 25 - ((N) & 0xF), ((N) & 0xF0) ? 0 : 26 - ((N) & 0xF),
                ((N) & 0xF0) ? 0 : 27 - ((N) & 0xF), ((N) & 0xF0) ? 0 : 28 - ((N) & 0xF),
                ((N) & 0xF0) ? 0 : 29 - ((N) & 0xF), ((N) & 0xF0) ? 0 : 30 - ((N) & 0xF),
                ((N) & 0xF0) ? 0 : 31 - ((N) & 0xF));
        }

        // slide_right
        template <size_t N, class A, class T>
        XSIMD_INLINE batch<T, A> slide_right(batch<T, A> const& x, requires_arch<wasm>) noexcept
        {
            return wasm_i8x16_shuffle(
                x, wasm_i64x2_const(0, 0), ((N) & 0xF0) ? 16 : ((N) & 0xF) + 0,
                ((N) & 0xF0) ? 16 : ((N) & 0xF) + 1, ((N) & 0xF0) ? 16 : ((N) & 0xF) + 2,
                ((N) & 0xF0) ? 16 : ((N) & 0xF) + 3, ((N) & 0xF0) ? 16 : ((N) & 0xF) + 4,
                ((N) & 0xF0) ? 16 : ((N) & 0xF) + 5, ((N) & 0xF0) ? 16 : ((N) & 0xF) + 6,
                ((N) & 0xF0) ? 16 : ((N) & 0xF) + 7, ((N) & 0xF0) ? 16 : ((N) & 0xF) + 8,
                ((N) & 0xF0) ? 16 : ((N) & 0xF) + 9, ((N) & 0xF0) ? 16 : ((N) & 0xF) + 10,
                ((N) & 0xF0) ? 16 : ((N) & 0xF) + 11, ((N) & 0xF0) ? 16 : ((N) & 0xF) + 12,
                ((N) & 0xF0) ? 16 : ((N) & 0xF) + 13, ((N) & 0xF0) ? 16 : ((N) & 0xF) + 14,
                ((N) & 0xF0) ? 16 : ((N) & 0xF) + 15);
        }

        // sadd
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> sadd(batch<T, A> const& self, batch<T, A> const& other, requires_arch<wasm>) noexcept
        {
            if (std::is_signed<T>::value)
            {
                XSIMD_IF_CONSTEXPR(sizeof(T) == 1)
                {
                    return wasm_i8x16_add_sat(self, other);
                }
                else XSIMD_IF_CONSTEXPR(sizeof(T) == 2)
                {
                    return wasm_i16x8_add_sat(self, other);
                }
                else
                {
                    return sadd(self, other, generic {});
                }
            }
            else
            {
                XSIMD_IF_CONSTEXPR(sizeof(T) == 1)
                {
                    return wasm_u8x16_add_sat(self, other);
                }
                else XSIMD_IF_CONSTEXPR(sizeof(T) == 2)
                {
                    return wasm_u16x8_add_sat(self, other);
                }
                else
                {
                    return sadd(self, other, generic {});
                }
            }
        }

        // select
        template <class A>
        XSIMD_INLINE batch<float, A> select(batch_bool<float, A> const& cond, batch<float, A> const& true_br, batch<float, A> const& false_br, requires_arch<wasm>) noexcept
        {
            return wasm_v128_or(wasm_v128_and(cond, true_br), wasm_v128_andnot(false_br, cond));
        }

        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> select(batch_bool<T, A> const& cond, batch<T, A> const& true_br, batch<T, A> const& false_br, requires_arch<wasm>) noexcept
        {
            return wasm_v128_or(wasm_v128_and(cond, true_br), wasm_v128_andnot(false_br, cond));
        }
        template <class A, class T, bool... Values, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> select(batch_bool_constant<T, A, Values...> const&, batch<T, A> const& true_br, batch<T, A> const& false_br, requires_arch<wasm>) noexcept
        {
            return select(batch_bool<T, A> { Values... }, true_br, false_br, wasm {});
        }
        template <class A>
        XSIMD_INLINE batch<double, A> select(batch_bool<double, A> const& cond, batch<double, A> const& true_br, batch<double, A> const& false_br, requires_arch<wasm>) noexcept
        {
            return wasm_v128_or(wasm_v128_and(cond, true_br), wasm_v128_andnot(false_br, cond));
        }

        // shuffle
        template <class A, class ITy, ITy I0, ITy I1, ITy I2, ITy I3>
        XSIMD_INLINE batch<float, A> shuffle(batch<float, A> const& x, batch<float, A> const& y, batch_constant<ITy, A, I0, I1, I2, I3>, requires_arch<wasm>) noexcept
        {
            return wasm_i32x4_shuffle(x, y, I0, I1, I2, I3);
        }

        template <class A, class ITy, ITy I0, ITy I1>
        XSIMD_INLINE batch<double, A> shuffle(batch<double, A> const& x, batch<double, A> const& y, batch_constant<ITy, A, I0, I1>, requires_arch<wasm>) noexcept
        {
            return wasm_i64x2_shuffle(x, y, I0, I1);
        }

        // set
        template <class A, class... Values>
        XSIMD_INLINE batch<float, A> set(batch<float, A> const&, requires_arch<wasm>, Values... values) noexcept
        {
            static_assert(sizeof...(Values) == batch<float, A>::size, "consistent init");
            return wasm_f32x4_make(values...);
        }

        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> set(batch<T, A> const&, requires_arch<wasm>, T v0, T v1) noexcept
        {
            return wasm_i64x2_make(v0, v1);
        }

        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> set(batch<T, A> const&, requires_arch<wasm>, T v0, T v1, T v2, T v3) noexcept
        {
            return wasm_i32x4_make(v0, v1, v2, v3);
        }

        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> set(batch<T, A> const&, requires_arch<wasm>, T v0, T v1, T v2, T v3, T v4, T v5, T v6, T v7) noexcept
        {
            return wasm_i16x8_make(v0, v1, v2, v3, v4, v5, v6, v7);
        }

        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> set(batch<T, A> const&, requires_arch<wasm>, T v0, T v1, T v2, T v3, T v4, T v5, T v6, T v7, T v8, T v9, T v10, T v11, T v12, T v13, T v14, T v15) noexcept
        {
            return wasm_i8x16_make(v0, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, v15);
        }

        template <class A, class... Values>
        XSIMD_INLINE batch<double, A> set(batch<double, A> const&, requires_arch<wasm>, Values... values) noexcept
        {
            static_assert(sizeof...(Values) == batch<double, A>::size, "consistent init");
            return wasm_f64x2_make(values...);
        }

        template <class A, class T, class... Values, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch_bool<T, A> set(batch_bool<T, A> const&, requires_arch<wasm>, Values... values) noexcept
        {
            return set(batch<T, A>(), A {}, static_cast<T>(values ? -1LL : 0LL)...).data;
        }

        template <class A, class... Values>
        XSIMD_INLINE batch_bool<float, A> set(batch_bool<float, A> const&, requires_arch<wasm>, Values... values) noexcept
        {
            static_assert(sizeof...(Values) == batch_bool<float, A>::size, "consistent init");
            return set(batch<int32_t, A>(), A {}, static_cast<int32_t>(values ? -1LL : 0LL)...).data;
        }

        template <class A, class... Values>
        XSIMD_INLINE batch_bool<double, A> set(batch_bool<double, A> const&, requires_arch<wasm>, Values... values) noexcept
        {
            static_assert(sizeof...(Values) == batch_bool<double, A>::size, "consistent init");
            return set(batch<int64_t, A>(), A {}, static_cast<int64_t>(values ? -1LL : 0LL)...).data;
        }

        // ssub
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> ssub(batch<T, A> const& self, batch<T, A> const& other, requires_arch<wasm>) noexcept
        {
            if (std::is_signed<T>::value)
            {
                XSIMD_IF_CONSTEXPR(sizeof(T) == 1)
                {
                    return wasm_i8x16_sub_sat(self, other);
                }
                else XSIMD_IF_CONSTEXPR(sizeof(T) == 2)
                {
                    return wasm_i16x8_sub_sat(self, other);
                }
                else
                {
                    return ssub(self, other, generic {});
                }
            }
            else
            {
                XSIMD_IF_CONSTEXPR(sizeof(T) == 1)
                {
                    return wasm_u8x16_sub_sat(self, other);
                }
                else XSIMD_IF_CONSTEXPR(sizeof(T) == 2)
                {
                    return wasm_u16x8_sub_sat(self, other);
                }
                else
                {
                    return ssub(self, other, generic {});
                }
            }
        }

        // store_aligned
        template <class A>
        XSIMD_INLINE void store_aligned(float* mem, batch<float, A> const& self, requires_arch<wasm>) noexcept
        {
            return wasm_v128_store(mem, self);
        }
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE void store_aligned(T* mem, batch<T, A> const& self, requires_arch<wasm>) noexcept
        {
            return wasm_v128_store((v128_t*)mem, self);
        }
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE void store_aligned(T* mem, batch_bool<T, A> const& self, requires_arch<wasm>) noexcept
        {
            return wasm_v128_store((v128_t*)mem, self);
        }
        template <class A>
        XSIMD_INLINE void store_aligned(double* mem, batch<double, A> const& self, requires_arch<wasm>) noexcept
        {
            return wasm_v128_store(mem, self);
        }

        // store_complex
        namespace detail
        {
            // complex_low
            template <class A>
            XSIMD_INLINE batch<float, A> complex_low(batch<std::complex<float>, A> const& self, requires_arch<wasm>) noexcept
            {
                return wasm_i32x4_shuffle(self.real(), self.imag(), 0, 4, 1, 5);
            }
            // complex_high
            template <class A>
            XSIMD_INLINE batch<float, A> complex_high(batch<std::complex<float>, A> const& self, requires_arch<wasm>) noexcept
            {
                return wasm_i32x4_shuffle(self.real(), self.imag(), 2, 6, 3, 7);
            }
            template <class A>
            XSIMD_INLINE batch<double, A> complex_low(batch<std::complex<double>, A> const& self, requires_arch<wasm>) noexcept
            {
                return wasm_i64x2_shuffle(self.real(), self.imag(), 0, 2);
            }
            template <class A>
            XSIMD_INLINE batch<double, A> complex_high(batch<std::complex<double>, A> const& self, requires_arch<wasm>) noexcept
            {
                return wasm_i64x2_shuffle(self.real(), self.imag(), 1, 3);
            }
        }

        // store_unaligned
        template <class A>
        XSIMD_INLINE void store_unaligned(float* mem, batch<float, A> const& self, requires_arch<wasm>) noexcept
        {
            return wasm_v128_store(mem, self);
        }
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE void store_unaligned(T* mem, batch<T, A> const& self, requires_arch<wasm>) noexcept
        {
            return wasm_v128_store((v128_t*)mem, self);
        }
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE void store_unaligned(T* mem, batch_bool<T, A> const& self, requires_arch<wasm>) noexcept
        {
            return wasm_v128_store((v128_t*)mem, self);
        }
        template <class A>
        XSIMD_INLINE void store_unaligned(double* mem, batch<double, A> const& self, requires_arch<wasm>) noexcept
        {
            return wasm_v128_store(mem, self);
        }

        // sub
        template <class A>
        XSIMD_INLINE batch<float, A> sub(batch<float, A> const& self, batch<float, A> const& other, requires_arch<wasm>) noexcept
        {
            return wasm_f32x4_sub(self, other);
        }
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> sub(batch<T, A> const& self, batch<T, A> const& other, requires_arch<wasm>) noexcept
        {
            XSIMD_IF_CONSTEXPR(sizeof(T) == 1)
            {
                return wasm_i8x16_sub(self, other);
            }
            else XSIMD_IF_CONSTEXPR(sizeof(T) == 2)
            {
                return wasm_i16x8_sub(self, other);
            }
            else XSIMD_IF_CONSTEXPR(sizeof(T) == 4)
            {
                return wasm_i32x4_sub(self, other);
            }
            else XSIMD_IF_CONSTEXPR(sizeof(T) == 8)
            {
                return wasm_i64x2_sub(self, other);
            }
            else
            {
                assert(false && "unsupported arch/op combination");
                return {};
            }
        }
        template <class A>
        XSIMD_INLINE batch<double, A> sub(batch<double, A> const& self, batch<double, A> const& other, requires_arch<wasm>) noexcept
        {
            return wasm_f64x2_sub(self, other);
        }

        // sqrt
        template <class A>
        XSIMD_INLINE batch<float, A> sqrt(batch<float, A> const& val, requires_arch<wasm>) noexcept
        {
            return wasm_f32x4_sqrt(val);
        }
        template <class A>
        XSIMD_INLINE batch<double, A> sqrt(batch<double, A> const& val, requires_arch<wasm>) noexcept
        {
            return wasm_f64x2_sqrt(val);
        }

        // swizzle
        template <class A, uint32_t V0, uint32_t V1, uint32_t V2, uint32_t V3>
        XSIMD_INLINE batch<float, A> swizzle(batch<float, A> const& self, batch_constant<uint32_t, A, V0, V1, V2, V3>, requires_arch<wasm>) noexcept
        {
            return wasm_i32x4_shuffle(self, self, V0, V1, V2, V3);
        }

        template <class A, uint64_t V0, uint64_t V1>
        XSIMD_INLINE batch<double, A> swizzle(batch<double, A> const& self, batch_constant<uint64_t, A, V0, V1>, requires_arch<wasm>) noexcept
        {
            return wasm_i64x2_shuffle(self, self, V0, V1);
        }

        template <class A, uint64_t V0, uint64_t V1>
        XSIMD_INLINE batch<uint64_t, A> swizzle(batch<uint64_t, A> const& self, batch_constant<uint64_t, A, V0, V1>, requires_arch<wasm>) noexcept
        {
            return wasm_i64x2_shuffle(self, self, V0, V1);
        }

        template <class A, uint64_t V0, uint64_t V1>
        XSIMD_INLINE batch<int64_t, A> swizzle(batch<int64_t, A> const& self, batch_constant<uint64_t, A, V0, V1> mask, requires_arch<wasm>) noexcept
        {
            return bitwise_cast<int64_t>(swizzle(bitwise_cast<uint64_t>(self), mask, wasm {}));
        }

        template <class A, uint32_t V0, uint32_t V1, uint32_t V2, uint32_t V3>
        XSIMD_INLINE batch<uint32_t, A> swizzle(batch<uint32_t, A> const& self, batch_constant<uint32_t, A, V0, V1, V2, V3>, requires_arch<wasm>) noexcept
        {
            return wasm_i32x4_shuffle(self, self, V0, V1, V2, V3);
        }

        template <class A, uint32_t V0, uint32_t V1, uint32_t V2, uint32_t V3>
        XSIMD_INLINE batch<int32_t, A> swizzle(batch<int32_t, A> const& self, batch_constant<uint32_t, A, V0, V1, V2, V3> mask, requires_arch<wasm>) noexcept
        {
            return bitwise_cast<int32_t>(swizzle(bitwise_cast<uint32_t>(self), mask, wasm {}));
        }

        template <class A, uint16_t V0, uint16_t V1, uint16_t V2, uint16_t V3, uint16_t V4, uint16_t V5, uint16_t V6, uint16_t V7>
        XSIMD_INLINE batch<uint16_t, A> swizzle(batch<uint16_t, A> const& self, batch_constant<uint16_t, A, V0, V1, V2, V3, V4, V5, V6, V7>, requires_arch<wasm>) noexcept
        {
            return wasm_i16x8_shuffle(self, self, V0, V1, V2, V3, V4, V5, V6, V7);
        }

        template <class A, uint16_t V0, uint16_t V1, uint16_t V2, uint16_t V3, uint16_t V4, uint16_t V5, uint16_t V6, uint16_t V7>
        XSIMD_INLINE batch<int16_t, A> swizzle(batch<int16_t, A> const& self, batch_constant<uint16_t, A, V0, V1, V2, V3, V4, V5, V6, V7> mask, requires_arch<wasm>) noexcept
        {
            return bitwise_cast<int16_t>(swizzle(bitwise_cast<uint16_t>(self), mask, wasm {}));
        }

        template <class A, uint8_t V0, uint8_t V1, uint8_t V2, uint8_t V3, uint8_t V4, uint8_t V5, uint8_t V6, uint8_t V7,
                  uint8_t V8, uint8_t V9, uint8_t V10, uint8_t V11, uint8_t V12, uint8_t V13, uint8_t V14, uint8_t V15>
        XSIMD_INLINE batch<uint8_t, A> swizzle(batch<uint8_t, A> const& self, batch_constant<uint8_t, A, V0, V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13, V14, V15>, requires_arch<wasm>) noexcept
        {
            return wasm_i8x16_shuffle(self, self, V0, V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13, V14, V15);
        }

        template <class A, uint8_t V0, uint8_t V1, uint8_t V2, uint8_t V3, uint8_t V4, uint8_t V5, uint8_t V6, uint8_t V7,
                  uint8_t V8, uint8_t V9, uint8_t V10, uint8_t V11, uint8_t V12, uint8_t V13, uint8_t V14, uint8_t V15>
        XSIMD_INLINE batch<int8_t, A> swizzle(batch<int8_t, A> const& self, batch_constant<uint8_t, A, V0, V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13, V14, V15> mask, requires_arch<wasm>) noexcept
        {
            return bitwise_cast<int8_t>(swizzle(bitwise_cast<uint8_t>(self), mask, wasm {}));
        }

        // trunc
        template <class A>
        XSIMD_INLINE batch<float, A> trunc(batch<float, A> const& self, requires_arch<wasm>) noexcept
        {
            return wasm_f32x4_trunc(self);
        }
        template <class A>
        XSIMD_INLINE batch<double, A> trunc(batch<double, A> const& self, requires_arch<wasm>) noexcept
        {
            return wasm_f64x2_trunc(self);
        }

        // zip_hi
        template <class A>
        XSIMD_INLINE batch<float, A> zip_hi(batch<float, A> const& self, batch<float, A> const& other, requires_arch<wasm>) noexcept
        {
            return wasm_i32x4_shuffle(self, other, 2, 6, 3, 7);
        }
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> zip_hi(batch<T, A> const& self, batch<T, A> const& other, requires_arch<wasm>) noexcept
        {
            XSIMD_IF_CONSTEXPR(sizeof(T) == 1)
            {
                return wasm_i8x16_shuffle(self, other, 8, 24, 9, 25, 10, 26, 11, 27, 12, 28, 13, 29, 14, 30, 15, 31);
            }
            else XSIMD_IF_CONSTEXPR(sizeof(T) == 2)
            {
                return wasm_i16x8_shuffle(self, other, 4, 12, 5, 13, 6, 14, 7, 15);
            }
            else XSIMD_IF_CONSTEXPR(sizeof(T) == 4)
            {
                return wasm_i32x4_shuffle(self, other, 2, 6, 3, 7);
            }
            else XSIMD_IF_CONSTEXPR(sizeof(T) == 8)
            {
                return wasm_i64x2_shuffle(self, other, 1, 3);
            }
            else
            {
                assert(false && "unsupported arch/op combination");
                return {};
            }
        }
        template <class A>
        XSIMD_INLINE batch<double, A> zip_hi(batch<double, A> const& self, batch<double, A> const& other, requires_arch<wasm>) noexcept
        {
            return wasm_i64x2_shuffle(self, other, 1, 3);
        }

        // zip_lo
        template <class A>
        XSIMD_INLINE batch<float, A> zip_lo(batch<float, A> const& self, batch<float, A> const& other, requires_arch<wasm>) noexcept
        {
            return wasm_i32x4_shuffle(self, other, 0, 4, 1, 5);
        }
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> zip_lo(batch<T, A> const& self, batch<T, A> const& other, requires_arch<wasm>) noexcept
        {
            XSIMD_IF_CONSTEXPR(sizeof(T) == 1)
            {
                return wasm_i8x16_shuffle(self, other, 0, 16, 1, 17, 2, 18, 3, 19, 4, 20, 5, 21, 6, 22, 7, 23);
            }
            else XSIMD_IF_CONSTEXPR(sizeof(T) == 2)
            {
                return wasm_i16x8_shuffle(self, other, 0, 8, 1, 9, 2, 10, 3, 11);
            }
            else XSIMD_IF_CONSTEXPR(sizeof(T) == 4)
            {
                return wasm_i32x4_shuffle(self, other, 0, 4, 1, 5);
            }
            else XSIMD_IF_CONSTEXPR(sizeof(T) == 8)
            {
                return wasm_i64x2_shuffle(self, other, 0, 2);
            }
            else
            {
                assert(false && "unsupported arch/op combination");
                return {};
            }
        }
        template <class A>
        XSIMD_INLINE batch<double, A> zip_lo(batch<double, A> const& self, batch<double, A> const& other, requires_arch<wasm>) noexcept
        {
            return wasm_i64x2_shuffle(self, other, 0, 2);
        }
    }
}

#endif
