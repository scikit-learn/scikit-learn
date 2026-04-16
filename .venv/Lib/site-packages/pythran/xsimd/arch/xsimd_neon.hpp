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

#ifndef XSIMD_NEON_HPP
#define XSIMD_NEON_HPP

#include <algorithm>
#include <complex>
#include <tuple>
#include <type_traits>

#include "../types/xsimd_neon_register.hpp"
#include "../types/xsimd_utils.hpp"

// Wrap intrinsics so we can pass them as function pointers
// - OP: intrinsics name prefix, e.g., vorrq
// - RT: type traits to deduce intrinsics return types
#define WRAP_BINARY_UINT_EXCLUDING_64(OP, RT)                                     \
    namespace wrap                                                                \
    {                                                                             \
        XSIMD_INLINE RT<uint8x16_t> OP##_u8(uint8x16_t a, uint8x16_t b) noexcept  \
        {                                                                         \
            return ::OP##_u8(a, b);                                               \
        }                                                                         \
        XSIMD_INLINE RT<uint16x8_t> OP##_u16(uint16x8_t a, uint16x8_t b) noexcept \
        {                                                                         \
            return ::OP##_u16(a, b);                                              \
        }                                                                         \
        XSIMD_INLINE RT<uint32x4_t> OP##_u32(uint32x4_t a, uint32x4_t b) noexcept \
        {                                                                         \
            return ::OP##_u32(a, b);                                              \
        }                                                                         \
    }

#define WRAP_BINARY_INT_EXCLUDING_64(OP, RT)                                   \
    WRAP_BINARY_UINT_EXCLUDING_64(OP, RT)                                      \
    namespace wrap                                                             \
    {                                                                          \
        XSIMD_INLINE RT<int8x16_t> OP##_s8(int8x16_t a, int8x16_t b) noexcept  \
        {                                                                      \
            return ::OP##_s8(a, b);                                            \
        }                                                                      \
        XSIMD_INLINE RT<int16x8_t> OP##_s16(int16x8_t a, int16x8_t b) noexcept \
        {                                                                      \
            return ::OP##_s16(a, b);                                           \
        }                                                                      \
        XSIMD_INLINE RT<int32x4_t> OP##_s32(int32x4_t a, int32x4_t b) noexcept \
        {                                                                      \
            return ::OP##_s32(a, b);                                           \
        }                                                                      \
    }

#define WRAP_BINARY_INT(OP, RT)                                                   \
    WRAP_BINARY_INT_EXCLUDING_64(OP, RT)                                          \
    namespace wrap                                                                \
    {                                                                             \
        XSIMD_INLINE RT<uint64x2_t> OP##_u64(uint64x2_t a, uint64x2_t b) noexcept \
        {                                                                         \
            return ::OP##_u64(a, b);                                              \
        }                                                                         \
        XSIMD_INLINE RT<int64x2_t> OP##_s64(int64x2_t a, int64x2_t b) noexcept    \
        {                                                                         \
            return ::OP##_s64(a, b);                                              \
        }                                                                         \
    }

#define WRAP_BINARY_FLOAT(OP, RT)                                                    \
    namespace wrap                                                                   \
    {                                                                                \
        XSIMD_INLINE RT<float32x4_t> OP##_f32(float32x4_t a, float32x4_t b) noexcept \
        {                                                                            \
            return ::OP##_f32(a, b);                                                 \
        }                                                                            \
    }

#define WRAP_UNARY_INT_EXCLUDING_64(OP)                         \
    namespace wrap                                              \
    {                                                           \
        XSIMD_INLINE uint8x16_t OP##_u8(uint8x16_t a) noexcept  \
        {                                                       \
            return ::OP##_u8(a);                                \
        }                                                       \
        XSIMD_INLINE int8x16_t OP##_s8(int8x16_t a) noexcept    \
        {                                                       \
            return ::OP##_s8(a);                                \
        }                                                       \
        XSIMD_INLINE uint16x8_t OP##_u16(uint16x8_t a) noexcept \
        {                                                       \
            return ::OP##_u16(a);                               \
        }                                                       \
        XSIMD_INLINE int16x8_t OP##_s16(int16x8_t a) noexcept   \
        {                                                       \
            return ::OP##_s16(a);                               \
        }                                                       \
        XSIMD_INLINE uint32x4_t OP##_u32(uint32x4_t a) noexcept \
        {                                                       \
            return ::OP##_u32(a);                               \
        }                                                       \
        XSIMD_INLINE int32x4_t OP##_s32(int32x4_t a) noexcept   \
        {                                                       \
            return ::OP##_s32(a);                               \
        }                                                       \
    }

#define WRAP_UNARY_INT(OP)                                      \
    WRAP_UNARY_INT_EXCLUDING_64(OP)                             \
    namespace wrap                                              \
    {                                                           \
        XSIMD_INLINE uint64x2_t OP##_u64(uint64x2_t a) noexcept \
        {                                                       \
            return ::OP##_u64(a);                               \
        }                                                       \
        XSIMD_INLINE int64x2_t OP##_s64(int64x2_t a) noexcept   \
        {                                                       \
            return ::OP##_s64(a);                               \
        }                                                       \
    }

#define WRAP_UNARY_FLOAT(OP)                                      \
    namespace wrap                                                \
    {                                                             \
        XSIMD_INLINE float32x4_t OP##_f32(float32x4_t a) noexcept \
        {                                                         \
            return ::OP##_f32(a);                                 \
        }                                                         \
    }

// Dummy identity caster to ease coding
XSIMD_INLINE uint8x16_t vreinterpretq_u8_u8(uint8x16_t arg) noexcept { return arg; }
XSIMD_INLINE int8x16_t vreinterpretq_s8_s8(int8x16_t arg) noexcept { return arg; }
XSIMD_INLINE uint16x8_t vreinterpretq_u16_u16(uint16x8_t arg) noexcept { return arg; }
XSIMD_INLINE int16x8_t vreinterpretq_s16_s16(int16x8_t arg) noexcept { return arg; }
XSIMD_INLINE uint32x4_t vreinterpretq_u32_u32(uint32x4_t arg) noexcept { return arg; }
XSIMD_INLINE int32x4_t vreinterpretq_s32_s32(int32x4_t arg) noexcept { return arg; }
XSIMD_INLINE uint64x2_t vreinterpretq_u64_u64(uint64x2_t arg) noexcept { return arg; }
XSIMD_INLINE int64x2_t vreinterpretq_s64_s64(int64x2_t arg) noexcept { return arg; }
XSIMD_INLINE float32x4_t vreinterpretq_f32_f32(float32x4_t arg) noexcept { return arg; }

namespace xsimd
{
    template <typename T, class A, bool... Values>
    struct batch_bool_constant;

    namespace kernel
    {
        using namespace types;

        namespace detail
        {
            template <template <class> class return_type, class... T>
            struct neon_dispatcher_base
            {
                struct unary
                {
                    using container_type = std::tuple<return_type<T> (*)(T)...>;
                    const container_type m_func;

                    template <class U>
                    return_type<U> apply(U rhs) const noexcept
                    {
                        using func_type = return_type<U> (*)(U);
                        auto func = xsimd::detail::get<func_type>(m_func);
                        return func(rhs);
                    }
                };

                struct binary
                {
                    using container_type = std::tuple<return_type<T> (*)(T, T)...>;
                    const container_type m_func;

                    template <class U>
                    return_type<U> apply(U lhs, U rhs) const noexcept
                    {
                        using func_type = return_type<U> (*)(U, U);
                        auto func = xsimd::detail::get<func_type>(m_func);
                        return func(lhs, rhs);
                    }
                };
            };

            /***************************
             *  arithmetic dispatchers *
             ***************************/

            template <class T>
            using identity_return_type = T;

            template <class... T>
            struct neon_dispatcher_impl : neon_dispatcher_base<identity_return_type, T...>
            {
            };

            using neon_dispatcher = neon_dispatcher_impl<uint8x16_t, int8x16_t,
                                                         uint16x8_t, int16x8_t,
                                                         uint32x4_t, int32x4_t,
                                                         uint64x2_t, int64x2_t,
                                                         float32x4_t>;

            using excluding_int64_dispatcher = neon_dispatcher_impl<uint8x16_t, int8x16_t,
                                                                    uint16x8_t, int16x8_t,
                                                                    uint32x4_t, int32x4_t,
                                                                    float32x4_t>;

            using excluding_int64f32_dispatcher = neon_dispatcher_impl<uint8x16_t, int8x16_t,
                                                                       uint16x8_t, int16x8_t,
                                                                       uint32x4_t, int32x4_t>;

            /**************************
             * comparison dispatchers *
             **************************/

            template <class T>
            struct comp_return_type_impl;

            template <>
            struct comp_return_type_impl<uint8x16_t>
            {
                using type = uint8x16_t;
            };

            template <>
            struct comp_return_type_impl<int8x16_t>
            {
                using type = uint8x16_t;
            };

            template <>
            struct comp_return_type_impl<uint16x8_t>
            {
                using type = uint16x8_t;
            };

            template <>
            struct comp_return_type_impl<int16x8_t>
            {
                using type = uint16x8_t;
            };

            template <>
            struct comp_return_type_impl<uint32x4_t>
            {
                using type = uint32x4_t;
            };

            template <>
            struct comp_return_type_impl<int32x4_t>
            {
                using type = uint32x4_t;
            };

            template <>
            struct comp_return_type_impl<uint64x2_t>
            {
                using type = uint64x2_t;
            };

            template <>
            struct comp_return_type_impl<int64x2_t>
            {
                using type = uint64x2_t;
            };

            template <>
            struct comp_return_type_impl<float32x4_t>
            {
                using type = uint32x4_t;
            };

            template <class T>
            using comp_return_type = typename comp_return_type_impl<T>::type;

            template <class... T>
            struct neon_comp_dispatcher_impl : neon_dispatcher_base<comp_return_type, T...>
            {
            };

            using excluding_int64_comp_dispatcher = neon_comp_dispatcher_impl<uint8x16_t, int8x16_t,
                                                                              uint16x8_t, int16x8_t,
                                                                              uint32x4_t, int32x4_t,
                                                                              float32x4_t>;

            /**************************************
             * enabling / disabling metafunctions *
             **************************************/

            template <class T>
            using enable_neon_type_t = typename std::enable_if<std::is_integral<T>::value || std::is_same<T, float>::value,
                                                               int>::type;

            template <class T>
            using exclude_int64_neon_t
                = typename std::enable_if<(std::is_integral<T>::value && sizeof(T) != 8) || std::is_same<T, float>::value, int>::type;
        }

        /*************
         * broadcast *
         *************/

        template <class A, class T, detail::enable_sized_unsigned_t<T, 1> = 0>
        XSIMD_INLINE batch<T, A> broadcast(T val, requires_arch<neon>) noexcept
        {
            return vdupq_n_u8(uint8_t(val));
        }

        template <class A, class T, detail::enable_sized_signed_t<T, 1> = 0>
        XSIMD_INLINE batch<T, A> broadcast(T val, requires_arch<neon>) noexcept
        {
            return vdupq_n_s8(int8_t(val));
        }

        template <class A, class T, detail::enable_sized_unsigned_t<T, 2> = 0>
        XSIMD_INLINE batch<T, A> broadcast(T val, requires_arch<neon>) noexcept
        {
            return vdupq_n_u16(uint16_t(val));
        }

        template <class A, class T, detail::enable_sized_signed_t<T, 2> = 0>
        XSIMD_INLINE batch<T, A> broadcast(T val, requires_arch<neon>) noexcept
        {
            return vdupq_n_s16(int16_t(val));
        }

        template <class A, class T, detail::enable_sized_unsigned_t<T, 4> = 0>
        XSIMD_INLINE batch<T, A> broadcast(T val, requires_arch<neon>) noexcept
        {
            return vdupq_n_u32(uint32_t(val));
        }

        template <class A, class T, detail::enable_sized_signed_t<T, 4> = 0>
        XSIMD_INLINE batch<T, A> broadcast(T val, requires_arch<neon>) noexcept
        {
            return vdupq_n_s32(int32_t(val));
        }

        template <class A, class T, detail::enable_sized_unsigned_t<T, 8> = 0>
        XSIMD_INLINE batch<T, A> broadcast(T val, requires_arch<neon>) noexcept
        {
            return vdupq_n_u64(uint64_t(val));
        }

        template <class A, class T, detail::enable_sized_signed_t<T, 8> = 0>
        XSIMD_INLINE batch<T, A> broadcast(T val, requires_arch<neon>) noexcept
        {
            return vdupq_n_s64(int64_t(val));
        }

        template <class A>
        XSIMD_INLINE batch<float, A> broadcast(float val, requires_arch<neon>) noexcept
        {
            return vdupq_n_f32(val);
        }

        /*******
         * set *
         *******/

        template <class A, class T, class... Args, detail::enable_integral_t<T> = 0>
        XSIMD_INLINE batch<T, A> set(batch<T, A> const&, requires_arch<neon>, Args... args) noexcept
        {
            return xsimd::types::detail::neon_vector_type<T> { args... };
        }

        template <class A, class T, class... Args, detail::enable_integral_t<T> = 0>
        XSIMD_INLINE batch_bool<T, A> set(batch_bool<T, A> const&, requires_arch<neon>, Args... args) noexcept
        {
            using register_type = typename batch_bool<T, A>::register_type;
            using unsigned_type = as_unsigned_integer_t<T>;
            return register_type { static_cast<unsigned_type>(args ? -1LL : 0LL)... };
        }

        template <class A>
        XSIMD_INLINE batch<float, A> set(batch<float, A> const&, requires_arch<neon>, float f0, float f1, float f2, float f3) noexcept
        {
            return float32x4_t { f0, f1, f2, f3 };
        }

        template <class A>
        XSIMD_INLINE batch<std::complex<float>, A> set(batch<std::complex<float>, A> const&, requires_arch<neon>,
                                                       std::complex<float> c0, std::complex<float> c1,
                                                       std::complex<float> c2, std::complex<float> c3) noexcept
        {
            return batch<std::complex<float>, A>(float32x4_t { c0.real(), c1.real(), c2.real(), c3.real() },
                                                 float32x4_t { c0.imag(), c1.imag(), c2.imag(), c3.imag() });
        }

        template <class A, class... Args>
        XSIMD_INLINE batch_bool<float, A> set(batch_bool<float, A> const&, requires_arch<neon>, Args... args) noexcept
        {
            using register_type = typename batch_bool<float, A>::register_type;
            using unsigned_type = as_unsigned_integer_t<float>;
            return register_type { static_cast<unsigned_type>(args ? -1LL : 0LL)... };
        }

        /*************
         * from_bool *
         *************/

        template <class A, class T, detail::enable_sized_unsigned_t<T, 1> = 0>
        XSIMD_INLINE batch<T, A> from_bool(batch_bool<T, A> const& arg, requires_arch<neon>) noexcept
        {
            return vandq_u8(arg, vdupq_n_u8(1));
        }

        template <class A, class T, detail::enable_sized_signed_t<T, 1> = 0>
        XSIMD_INLINE batch<T, A> from_bool(batch_bool<T, A> const& arg, requires_arch<neon>) noexcept
        {
            return vandq_s8(reinterpret_cast<int8x16_t>(arg.data), vdupq_n_s8(1));
        }

        template <class A, class T, detail::enable_sized_unsigned_t<T, 2> = 0>
        XSIMD_INLINE batch<T, A> from_bool(batch_bool<T, A> const& arg, requires_arch<neon>) noexcept
        {
            return vandq_u16(arg, vdupq_n_u16(1));
        }

        template <class A, class T, detail::enable_sized_signed_t<T, 2> = 0>
        XSIMD_INLINE batch<T, A> from_bool(batch_bool<T, A> const& arg, requires_arch<neon>) noexcept
        {
            return vandq_s16(reinterpret_cast<int16x8_t>(arg.data), vdupq_n_s16(1));
        }

        template <class A, class T, detail::enable_sized_unsigned_t<T, 4> = 0>
        XSIMD_INLINE batch<T, A> from_bool(batch_bool<T, A> const& arg, requires_arch<neon>) noexcept
        {
            return vandq_u32(arg, vdupq_n_u32(1));
        }

        template <class A, class T, detail::enable_sized_signed_t<T, 4> = 0>
        XSIMD_INLINE batch<T, A> from_bool(batch_bool<T, A> const& arg, requires_arch<neon>) noexcept
        {
            return vandq_s32(reinterpret_cast<int32x4_t>(arg.data), vdupq_n_s32(1));
        }

        template <class A, class T, detail::enable_sized_unsigned_t<T, 8> = 0>
        XSIMD_INLINE batch<T, A> from_bool(batch_bool<T, A> const& arg, requires_arch<neon>) noexcept
        {
            return vandq_u64(arg, vdupq_n_u64(1));
        }

        template <class A, class T, detail::enable_sized_signed_t<T, 8> = 0>
        XSIMD_INLINE batch<T, A> from_bool(batch_bool<T, A> const& arg, requires_arch<neon>) noexcept
        {
            return vandq_s64(reinterpret_cast<int64x2_t>(arg.data), vdupq_n_s64(1));
        }

        template <class A>
        XSIMD_INLINE batch<float, A> from_bool(batch_bool<float, A> const& arg, requires_arch<neon>) noexcept
        {
            return vreinterpretq_f32_u32(vandq_u32(arg, vreinterpretq_u32_f32(vdupq_n_f32(1.f))));
        }

        /********
         * load *
         ********/

        // It is not possible to use a call to A::alignment() here, so use an
        // immediate instead.
#if defined(__clang__) || defined(__GNUC__)
#define xsimd_aligned_load(inst, type, expr) inst((type)__builtin_assume_aligned(expr, 16))
#elif defined(_MSC_VER)
#define xsimd_aligned_load(inst, type, expr) inst##_ex((type)expr, 128)
#else
#define xsimd_aligned_load(inst, type, expr) inst((type)expr)
#endif

        template <class A, class T, detail::enable_sized_unsigned_t<T, 1> = 0>
        XSIMD_INLINE batch<T, A> load_aligned(T const* src, convert<T>, requires_arch<neon>) noexcept
        {
            return xsimd_aligned_load(vld1q_u8, uint8_t*, src);
        }

        template <class A, class T, detail::enable_sized_signed_t<T, 1> = 0>
        XSIMD_INLINE batch<T, A> load_aligned(T const* src, convert<T>, requires_arch<neon>) noexcept
        {
            return xsimd_aligned_load(vld1q_s8, int8_t*, src);
        }

        template <class A, class T, detail::enable_sized_unsigned_t<T, 2> = 0>
        XSIMD_INLINE batch<T, A> load_aligned(T const* src, convert<T>, requires_arch<neon>) noexcept
        {
            return xsimd_aligned_load(vld1q_u16, uint16_t*, src);
        }
        template <class A, class T, detail::enable_sized_signed_t<T, 2> = 0>
        XSIMD_INLINE batch<T, A> load_aligned(T const* src, convert<T>, requires_arch<neon>) noexcept
        {
            return xsimd_aligned_load(vld1q_s16, int16_t*, src);
        }
        template <class A, class T, detail::enable_sized_unsigned_t<T, 4> = 0>
        XSIMD_INLINE batch<T, A> load_aligned(T const* src, convert<T>, requires_arch<neon>) noexcept
        {
            return xsimd_aligned_load(vld1q_u32, uint32_t*, src);
        }
        template <class A, class T, detail::enable_sized_signed_t<T, 4> = 0>
        XSIMD_INLINE batch<T, A> load_aligned(T const* src, convert<T>, requires_arch<neon>) noexcept
        {
            return xsimd_aligned_load(vld1q_s32, int32_t*, src);
        }
        template <class A, class T, detail::enable_sized_unsigned_t<T, 8> = 0>
        XSIMD_INLINE batch<T, A> load_aligned(T const* src, convert<T>, requires_arch<neon>) noexcept
        {
            return xsimd_aligned_load(vld1q_u64, uint64_t*, src);
        }
        template <class A, class T, detail::enable_sized_signed_t<T, 8> = 0>
        XSIMD_INLINE batch<T, A> load_aligned(T const* src, convert<T>, requires_arch<neon>) noexcept
        {
            return xsimd_aligned_load(vld1q_s64, int64_t*, src);
        }

        template <class A>
        XSIMD_INLINE batch<float, A> load_aligned(float const* src, convert<float>, requires_arch<neon>) noexcept
        {
            return xsimd_aligned_load(vld1q_f32, float*, src);
        }

#undef xsimd_aligned_load

        template <class A, class T, detail::enable_sized_unsigned_t<T, 1> = 0>
        XSIMD_INLINE batch<T, A> load_unaligned(T const* src, convert<T>, requires_arch<neon>) noexcept
        {
            return vld1q_u8((uint8_t*)src);
        }

        template <class A, class T, detail::enable_sized_signed_t<T, 1> = 0>
        XSIMD_INLINE batch<T, A> load_unaligned(T const* src, convert<T>, requires_arch<neon>) noexcept
        {
            return vld1q_s8((int8_t*)src);
        }

        template <class A, class T, detail::enable_sized_unsigned_t<T, 2> = 0>
        XSIMD_INLINE batch<T, A> load_unaligned(T const* src, convert<T>, requires_arch<neon>) noexcept
        {
            return vld1q_u16((uint16_t*)src);
        }
        template <class A, class T, detail::enable_sized_signed_t<T, 2> = 0>
        XSIMD_INLINE batch<T, A> load_unaligned(T const* src, convert<T>, requires_arch<neon>) noexcept
        {
            return vld1q_s16((int16_t*)src);
        }
        template <class A, class T, detail::enable_sized_unsigned_t<T, 4> = 0>
        XSIMD_INLINE batch<T, A> load_unaligned(T const* src, convert<T>, requires_arch<neon>) noexcept
        {
            return vld1q_u32((uint32_t*)src);
        }
        template <class A, class T, detail::enable_sized_signed_t<T, 4> = 0>
        XSIMD_INLINE batch<T, A> load_unaligned(T const* src, convert<T>, requires_arch<neon>) noexcept
        {
            return vld1q_s32((int32_t*)src);
        }
        template <class A, class T, detail::enable_sized_unsigned_t<T, 8> = 0>
        XSIMD_INLINE batch<T, A> load_unaligned(T const* src, convert<T>, requires_arch<neon>) noexcept
        {
            return vld1q_u64((uint64_t*)src);
        }
        template <class A, class T, detail::enable_sized_signed_t<T, 8> = 0>
        XSIMD_INLINE batch<T, A> load_unaligned(T const* src, convert<T>, requires_arch<neon>) noexcept
        {
            return vld1q_s64((int64_t*)src);
        }

        template <class A>
        XSIMD_INLINE batch<float, A> load_unaligned(float const* src, convert<float>, requires_arch<neon>) noexcept
        {
            return vld1q_f32(src);
        }

        /*********
         * store *
         *********/

        template <class A, class T, detail::enable_sized_unsigned_t<T, 1> = 0>
        XSIMD_INLINE void store_aligned(T* dst, batch<T, A> const& src, requires_arch<neon>) noexcept
        {
            vst1q_u8((uint8_t*)dst, src);
        }

        template <class A, class T, detail::enable_sized_signed_t<T, 1> = 0>
        XSIMD_INLINE void store_aligned(T* dst, batch<T, A> const& src, requires_arch<neon>) noexcept
        {
            vst1q_s8((int8_t*)dst, src);
        }

        template <class A, class T, detail::enable_sized_unsigned_t<T, 2> = 0>
        XSIMD_INLINE void store_aligned(T* dst, batch<T, A> const& src, requires_arch<neon>) noexcept
        {
            vst1q_u16((uint16_t*)dst, src);
        }

        template <class A, class T, detail::enable_sized_signed_t<T, 2> = 0>
        XSIMD_INLINE void store_aligned(T* dst, batch<T, A> const& src, requires_arch<neon>) noexcept
        {
            vst1q_s16((int16_t*)dst, src);
        }

        template <class A, class T, detail::enable_sized_unsigned_t<T, 4> = 0>
        XSIMD_INLINE void store_aligned(T* dst, batch<T, A> const& src, requires_arch<neon>) noexcept
        {
            vst1q_u32((uint32_t*)dst, src);
        }

        template <class A, class T, detail::enable_sized_signed_t<T, 4> = 0>
        XSIMD_INLINE void store_aligned(T* dst, batch<T, A> const& src, requires_arch<neon>) noexcept
        {
            vst1q_s32((int32_t*)dst, src);
        }

        template <class A, class T, detail::enable_sized_unsigned_t<T, 8> = 0>
        XSIMD_INLINE void store_aligned(T* dst, batch<T, A> const& src, requires_arch<neon>) noexcept
        {
            vst1q_u64((uint64_t*)dst, src);
        }

        template <class A, class T, detail::enable_sized_signed_t<T, 8> = 0>
        XSIMD_INLINE void store_aligned(T* dst, batch<T, A> const& src, requires_arch<neon>) noexcept
        {
            vst1q_s64((int64_t*)dst, src);
        }

        template <class A>
        XSIMD_INLINE void store_aligned(float* dst, batch<float, A> const& src, requires_arch<neon>) noexcept
        {
            vst1q_f32(dst, src);
        }

        template <class A, class T>
        XSIMD_INLINE void store_unaligned(T* dst, batch<T, A> const& src, requires_arch<neon>) noexcept
        {
            store_aligned<A>(dst, src, A {});
        }

        /****************
         * load_complex *
         ****************/

        template <class A>
        XSIMD_INLINE batch<std::complex<float>, A> load_complex_aligned(std::complex<float> const* mem, convert<std::complex<float>>, requires_arch<neon>) noexcept
        {
            using real_batch = batch<float, A>;
            const float* buf = reinterpret_cast<const float*>(mem);
            float32x4x2_t tmp = vld2q_f32(buf);
            real_batch real = tmp.val[0],
                       imag = tmp.val[1];
            return batch<std::complex<float>, A> { real, imag };
        }

        template <class A>
        XSIMD_INLINE batch<std::complex<float>, A> load_complex_unaligned(std::complex<float> const* mem, convert<std::complex<float>> cvt, requires_arch<neon>) noexcept
        {
            return load_complex_aligned<A>(mem, cvt, A {});
        }

        /*****************
         * store_complex *
         *****************/

        template <class A>
        XSIMD_INLINE void store_complex_aligned(std::complex<float>* dst, batch<std::complex<float>, A> const& src, requires_arch<neon>) noexcept
        {
            float32x4x2_t tmp;
            tmp.val[0] = src.real();
            tmp.val[1] = src.imag();
            float* buf = reinterpret_cast<float*>(dst);
            vst2q_f32(buf, tmp);
        }

        template <class A>
        XSIMD_INLINE void store_complex_unaligned(std::complex<float>* dst, batch<std::complex<float>, A> const& src, requires_arch<neon>) noexcept
        {
            store_complex_aligned(dst, src, A {});
        }

        /*******
         * neg *
         *******/

        template <class A, class T, detail::enable_sized_unsigned_t<T, 1> = 0>
        XSIMD_INLINE batch<T, A> neg(batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            return vreinterpretq_u8_s8(vnegq_s8(vreinterpretq_s8_u8(rhs)));
        }

        template <class A, class T, detail::enable_sized_signed_t<T, 1> = 0>
        XSIMD_INLINE batch<T, A> neg(batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            return vnegq_s8(rhs);
        }

        template <class A, class T, detail::enable_sized_unsigned_t<T, 2> = 0>
        XSIMD_INLINE batch<T, A> neg(batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            return vreinterpretq_u16_s16(vnegq_s16(vreinterpretq_s16_u16(rhs)));
        }

        template <class A, class T, detail::enable_sized_signed_t<T, 2> = 0>
        XSIMD_INLINE batch<T, A> neg(batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            return vnegq_s16(rhs);
        }

        template <class A, class T, detail::enable_sized_unsigned_t<T, 4> = 0>
        XSIMD_INLINE batch<T, A> neg(batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            return vreinterpretq_u32_s32(vnegq_s32(vreinterpretq_s32_u32(rhs)));
        }

        template <class A, class T, detail::enable_sized_signed_t<T, 4> = 0>
        XSIMD_INLINE batch<T, A> neg(batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            return vnegq_s32(rhs);
        }

        template <class A, class T, detail::enable_sized_unsigned_t<T, 8> = 0>
        XSIMD_INLINE batch<T, A> neg(batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            return batch<T, A> { -rhs.get(0), -rhs.get(1) };
        }

        template <class A, class T, detail::enable_sized_signed_t<T, 8> = 0>
        XSIMD_INLINE batch<T, A> neg(batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            return batch<T, A> { -rhs.get(0), -rhs.get(1) };
        }

        template <class A>
        XSIMD_INLINE batch<float, A> neg(batch<float, A> const& rhs, requires_arch<neon>) noexcept
        {
            return vnegq_f32(rhs);
        }

        /*******
         * add *
         *******/

        WRAP_BINARY_INT(vaddq, detail::identity_return_type)
        WRAP_BINARY_FLOAT(vaddq, detail::identity_return_type)

        template <class A, class T, detail::enable_neon_type_t<T> = 0>
        XSIMD_INLINE batch<T, A> add(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            using register_type = typename batch<T, A>::register_type;
            const detail::neon_dispatcher::binary dispatcher = {
                std::make_tuple(wrap::vaddq_u8, wrap::vaddq_s8, wrap::vaddq_u16, wrap::vaddq_s16,
                                wrap::vaddq_u32, wrap::vaddq_s32, wrap::vaddq_u64, wrap::vaddq_s64,
                                wrap::vaddq_f32)
            };
            return dispatcher.apply(register_type(lhs), register_type(rhs));
        }

        /*******
         * avg *
         *******/

        WRAP_BINARY_UINT_EXCLUDING_64(vhaddq, detail::identity_return_type)

        template <class A, class T, class = typename std::enable_if<(std::is_unsigned<T>::value && sizeof(T) != 8), void>::type>
        XSIMD_INLINE batch<T, A> avg(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            using register_type = typename batch<T, A>::register_type;
            const detail::neon_dispatcher_impl<uint8x16_t, uint16x8_t, uint32x4_t>::binary dispatcher = {
                std::make_tuple(wrap::vhaddq_u8, wrap::vhaddq_u16, wrap::vhaddq_u32)
            };
            return dispatcher.apply(register_type(lhs), register_type(rhs));
        }

        /********
         * avgr *
         ********/

        WRAP_BINARY_UINT_EXCLUDING_64(vrhaddq, detail::identity_return_type)

        template <class A, class T, class = typename std::enable_if<(std::is_unsigned<T>::value && sizeof(T) != 8), void>::type>
        XSIMD_INLINE batch<T, A> avgr(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            using register_type = typename batch<T, A>::register_type;
            const detail::neon_dispatcher_impl<uint8x16_t, uint16x8_t, uint32x4_t>::binary dispatcher = {
                std::make_tuple(wrap::vrhaddq_u8, wrap::vrhaddq_u16, wrap::vrhaddq_u32)
            };
            return dispatcher.apply(register_type(lhs), register_type(rhs));
        }

        /********
         * sadd *
         ********/

        WRAP_BINARY_INT(vqaddq, detail::identity_return_type)

        template <class A, class T, detail::enable_neon_type_t<T> = 0>
        XSIMD_INLINE batch<T, A> sadd(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            using register_type = typename batch<T, A>::register_type;
            const detail::neon_dispatcher::binary dispatcher = {
                std::make_tuple(wrap::vqaddq_u8, wrap::vqaddq_s8, wrap::vqaddq_u16, wrap::vqaddq_s16,
                                wrap::vqaddq_u32, wrap::vqaddq_s32, wrap::vqaddq_u64, wrap::vqaddq_s64,
                                wrap::vaddq_f32)
            };
            return dispatcher.apply(register_type(lhs), register_type(rhs));
        }

        /*******
         * sub *
         *******/

        WRAP_BINARY_INT(vsubq, detail::identity_return_type)
        WRAP_BINARY_FLOAT(vsubq, detail::identity_return_type)

        template <class A, class T, detail::enable_neon_type_t<T> = 0>
        XSIMD_INLINE batch<T, A> sub(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            using register_type = typename batch<T, A>::register_type;
            const detail::neon_dispatcher::binary dispatcher = {
                std::make_tuple(wrap::vsubq_u8, wrap::vsubq_s8, wrap::vsubq_u16, wrap::vsubq_s16,
                                wrap::vsubq_u32, wrap::vsubq_s32, wrap::vsubq_u64, wrap::vsubq_s64,
                                wrap::vsubq_f32)
            };
            return dispatcher.apply(register_type(lhs), register_type(rhs));
        }

        /********
         * ssub *
         ********/

        WRAP_BINARY_INT(vqsubq, detail::identity_return_type)

        template <class A, class T, detail::enable_neon_type_t<T> = 0>
        XSIMD_INLINE batch<T, A> ssub(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            using register_type = typename batch<T, A>::register_type;
            const detail::neon_dispatcher::binary dispatcher = {
                std::make_tuple(wrap::vqsubq_u8, wrap::vqsubq_s8, wrap::vqsubq_u16, wrap::vqsubq_s16,
                                wrap::vqsubq_u32, wrap::vqsubq_s32, wrap::vqsubq_u64, wrap::vqsubq_s64,
                                wrap::vsubq_f32)
            };
            return dispatcher.apply(register_type(lhs), register_type(rhs));
        }

        /*******
         * mul *
         *******/

        WRAP_BINARY_INT_EXCLUDING_64(vmulq, detail::identity_return_type)
        WRAP_BINARY_FLOAT(vmulq, detail::identity_return_type)

        template <class A, class T, detail::exclude_int64_neon_t<T> = 0>
        XSIMD_INLINE batch<T, A> mul(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            using register_type = typename batch<T, A>::register_type;
            const detail::excluding_int64_dispatcher::binary dispatcher = {
                std::make_tuple(wrap::vmulq_u8, wrap::vmulq_s8, wrap::vmulq_u16, wrap::vmulq_s16,
                                wrap::vmulq_u32, wrap::vmulq_s32, wrap::vmulq_f32)
            };
            return dispatcher.apply(register_type(lhs), register_type(rhs));
        }

        /*******
         * div *
         *******/

#if defined(XSIMD_FAST_INTEGER_DIVISION)
        template <class A, class T, detail::enable_sized_signed_t<T, 4> = 0>
        XSIMD_INLINE batch<T, A> div(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            return vcvtq_s32_f32(vcvtq_f32_s32(lhs) / vcvtq_f32_s32(rhs));
        }

        template <class A, class T, detail::enable_sized_unsigned_t<T, 4> = 0>
        XSIMD_INLINE batch<T, A> div(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            return vcvtq_u32_f32(vcvtq_f32_u32(lhs) / vcvtq_f32_u32(rhs));
        }
#endif

        template <class A>
        XSIMD_INLINE batch<float, A> div(batch<float, A> const& lhs, batch<float, A> const& rhs, requires_arch<neon>) noexcept
        {
            // from stackoverflow & https://projectne10.github.io/Ne10/doc/NE10__divc_8neon_8c_source.html
            // get an initial estimate of 1/b.
            float32x4_t rcp = reciprocal(rhs);

            // use a couple Newton-Raphson steps to refine the estimate.  Depending on your
            // application's accuracy requirements, you may be able to get away with only
            // one refinement (instead of the two used here).  Be sure to test!
            rcp = vmulq_f32(vrecpsq_f32(rhs, rcp), rcp);
            rcp = vmulq_f32(vrecpsq_f32(rhs, rcp), rcp);

            // and finally, compute a / b = a * (1 / b)
            return vmulq_f32(lhs, rcp);
        }

        /******
         * eq *
         ******/

        WRAP_BINARY_INT_EXCLUDING_64(vceqq, detail::comp_return_type)
        WRAP_BINARY_FLOAT(vceqq, detail::comp_return_type)

        template <class A, class T, detail::exclude_int64_neon_t<T> = 0>
        XSIMD_INLINE batch_bool<T, A> eq(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            using register_type = typename batch<T, A>::register_type;
            const detail::excluding_int64_comp_dispatcher::binary dispatcher = {
                std::make_tuple(wrap::vceqq_u8, wrap::vceqq_s8, wrap::vceqq_u16, wrap::vceqq_s16,
                                wrap::vceqq_u32, wrap::vceqq_s32, wrap::vceqq_f32)
            };
            return dispatcher.apply(register_type(lhs), register_type(rhs));
        }

        template <class A, class T, detail::exclude_int64_neon_t<T> = 0>
        XSIMD_INLINE batch_bool<T, A> eq(batch_bool<T, A> const& lhs, batch_bool<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            using register_type = typename batch_bool<T, A>::register_type;
            using dispatcher_type = detail::neon_comp_dispatcher_impl<uint8x16_t, uint16x8_t, uint32x4_t>::binary;
            const dispatcher_type dispatcher = {
                std::make_tuple(wrap::vceqq_u8, wrap::vceqq_u16, wrap::vceqq_u32)
            };
            return dispatcher.apply(register_type(lhs), register_type(rhs));
        }

        template <class A, class T, detail::enable_sized_integral_t<T, 8> = 0>
        XSIMD_INLINE batch_bool<T, A> eq(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            return batch_bool<T, A>({ lhs.get(0) == rhs.get(0), lhs.get(1) == rhs.get(1) });
        }

        template <class A, class T, detail::enable_sized_integral_t<T, 8> = 0>
        XSIMD_INLINE batch_bool<T, A> eq(batch_bool<T, A> const& lhs, batch_bool<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            return batch_bool<T, A>({ lhs.get(0) == rhs.get(0), lhs.get(1) == rhs.get(1) });
        }

        /*************
         * fast_cast *
         *************/

        namespace detail
        {
            template <class A>
            XSIMD_INLINE batch<float, A> fast_cast(batch<int32_t, A> const& self, batch<float, A> const&, requires_arch<neon>) noexcept
            {
                return vcvtq_f32_s32(self);
            }

            template <class A>
            XSIMD_INLINE batch<float, A> fast_cast(batch<uint32_t, A> const& self, batch<float, A> const&, requires_arch<neon>) noexcept
            {
                return vcvtq_f32_u32(self);
            }

            template <class A>
            XSIMD_INLINE batch<int32_t, A> fast_cast(batch<float, A> const& self, batch<int32_t, A> const&, requires_arch<neon>) noexcept
            {
                return vcvtq_s32_f32(self);
            }

            template <class A>
            XSIMD_INLINE batch<uint32_t, A> fast_cast(batch<float, A> const& self, batch<uint32_t, A> const&, requires_arch<neon>) noexcept
            {
                return vcvtq_u32_f32(self);
            }

        }

        /******
         * lt *
         ******/

        WRAP_BINARY_INT_EXCLUDING_64(vcltq, detail::comp_return_type)
        WRAP_BINARY_FLOAT(vcltq, detail::comp_return_type)

        template <class A, class T, detail::exclude_int64_neon_t<T> = 0>
        XSIMD_INLINE batch_bool<T, A> lt(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            using register_type = typename batch<T, A>::register_type;
            const detail::excluding_int64_comp_dispatcher::binary dispatcher = {
                std::make_tuple(wrap::vcltq_u8, wrap::vcltq_s8, wrap::vcltq_u16, wrap::vcltq_s16,
                                wrap::vcltq_u32, wrap::vcltq_s32, wrap::vcltq_f32)
            };
            return dispatcher.apply(register_type(lhs), register_type(rhs));
        }

        template <class A, class T, detail::enable_sized_integral_t<T, 8> = 0>
        XSIMD_INLINE batch_bool<T, A> lt(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            return batch_bool<T, A>({ lhs.get(0) < rhs.get(0), lhs.get(1) < rhs.get(1) });
        }

        /******
         * le *
         ******/

        WRAP_BINARY_INT_EXCLUDING_64(vcleq, detail::comp_return_type)
        WRAP_BINARY_FLOAT(vcleq, detail::comp_return_type)

        template <class A, class T, detail::exclude_int64_neon_t<T> = 0>
        XSIMD_INLINE batch_bool<T, A> le(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            using register_type = typename batch<T, A>::register_type;
            const detail::excluding_int64_comp_dispatcher::binary dispatcher = {
                std::make_tuple(wrap::vcleq_u8, wrap::vcleq_s8, wrap::vcleq_u16, wrap::vcleq_s16,
                                wrap::vcleq_u32, wrap::vcleq_s32, wrap::vcleq_f32)
            };
            return dispatcher.apply(register_type(lhs), register_type(rhs));
        }

        template <class A, class T, detail::enable_sized_integral_t<T, 8> = 0>
        XSIMD_INLINE batch_bool<T, A> le(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            return batch_bool<T, A>({ lhs.get(0) <= rhs.get(0), lhs.get(1) <= rhs.get(1) });
        }

        /******
         * gt *
         ******/

        WRAP_BINARY_INT_EXCLUDING_64(vcgtq, detail::comp_return_type)
        WRAP_BINARY_FLOAT(vcgtq, detail::comp_return_type)

        template <class A, class T, detail::exclude_int64_neon_t<T> = 0>
        XSIMD_INLINE batch_bool<T, A> gt(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            using register_type = typename batch<T, A>::register_type;
            const detail::excluding_int64_comp_dispatcher::binary dispatcher = {
                std::make_tuple(wrap::vcgtq_u8, wrap::vcgtq_s8, wrap::vcgtq_u16, wrap::vcgtq_s16,
                                wrap::vcgtq_u32, wrap::vcgtq_s32, wrap::vcgtq_f32)
            };
            return dispatcher.apply(register_type(lhs), register_type(rhs));
        }

        template <class A, class T, detail::enable_sized_integral_t<T, 8> = 0>
        XSIMD_INLINE batch_bool<T, A> gt(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            return batch_bool<T, A>({ lhs.get(0) > rhs.get(0), lhs.get(1) > rhs.get(1) });
        }

        /******
         * ge *
         ******/

        WRAP_BINARY_INT_EXCLUDING_64(vcgeq, detail::comp_return_type)
        WRAP_BINARY_FLOAT(vcgeq, detail::comp_return_type)

        template <class A, class T, detail::exclude_int64_neon_t<T> = 0>
        XSIMD_INLINE batch_bool<T, A> ge(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            using register_type = typename batch<T, A>::register_type;
            const detail::excluding_int64_comp_dispatcher::binary dispatcher = {
                std::make_tuple(wrap::vcgeq_u8, wrap::vcgeq_s8, wrap::vcgeq_u16, wrap::vcgeq_s16,
                                wrap::vcgeq_u32, wrap::vcgeq_s32, wrap::vcgeq_f32)
            };
            return dispatcher.apply(register_type(lhs), register_type(rhs));
        }

        template <class A, class T, detail::enable_sized_integral_t<T, 8> = 0>
        XSIMD_INLINE batch_bool<T, A> ge(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            return batch_bool<T, A>({ lhs.get(0) >= rhs.get(0), lhs.get(1) >= rhs.get(1) });
        }

        /*******************
         * batch_bool_cast *
         *******************/

        template <class A, class T_out, class T_in>
        XSIMD_INLINE batch_bool<T_out, A> batch_bool_cast(batch_bool<T_in, A> const& self, batch_bool<T_out, A> const&, requires_arch<neon>) noexcept
        {
            using register_type = typename batch_bool<T_out, A>::register_type;
            return register_type(self);
        }

        /***************
         * bitwise_and *
         ***************/

        WRAP_BINARY_INT(vandq, detail::identity_return_type)

        namespace detail
        {
            XSIMD_INLINE float32x4_t bitwise_and_f32(float32x4_t lhs, float32x4_t rhs) noexcept
            {
                return vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(lhs),
                                                       vreinterpretq_u32_f32(rhs)));
            }

            template <class V>
            V bitwise_and_neon(V const& lhs, V const& rhs)
            {
                const neon_dispatcher::binary dispatcher = {
                    std::make_tuple(wrap::vandq_u8, wrap::vandq_s8, wrap::vandq_u16, wrap::vandq_s16,
                                    wrap::vandq_u32, wrap::vandq_s32, wrap::vandq_u64, wrap::vandq_s64,
                                    bitwise_and_f32)
                };
                return dispatcher.apply(lhs, rhs);
            }
        }

        template <class A, class T, detail::enable_neon_type_t<T> = 0>
        XSIMD_INLINE batch<T, A> bitwise_and(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            using register_type = typename batch<T, A>::register_type;
            return detail::bitwise_and_neon(register_type(lhs), register_type(rhs));
        }

        template <class A, class T, detail::enable_neon_type_t<T> = 0>
        XSIMD_INLINE batch_bool<T, A> bitwise_and(batch_bool<T, A> const& lhs, batch_bool<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            using register_type = typename batch_bool<T, A>::register_type;
            return detail::bitwise_and_neon(register_type(lhs), register_type(rhs));
        }

        /**************
         * bitwise_or *
         **************/

        WRAP_BINARY_INT(vorrq, detail::identity_return_type)

        namespace detail
        {
            XSIMD_INLINE float32x4_t bitwise_or_f32(float32x4_t lhs, float32x4_t rhs) noexcept
            {
                return vreinterpretq_f32_u32(vorrq_u32(vreinterpretq_u32_f32(lhs),
                                                       vreinterpretq_u32_f32(rhs)));
            }

            template <class V>
            XSIMD_INLINE V bitwise_or_neon(V const& lhs, V const& rhs) noexcept
            {
                const neon_dispatcher::binary dispatcher = {
                    std::make_tuple(wrap::vorrq_u8, wrap::vorrq_s8, wrap::vorrq_u16, wrap::vorrq_s16,
                                    wrap::vorrq_u32, wrap::vorrq_s32, wrap::vorrq_u64, wrap::vorrq_s64,
                                    bitwise_or_f32)
                };
                return dispatcher.apply(lhs, rhs);
            }
        }

        template <class A, class T, detail::enable_neon_type_t<T> = 0>
        XSIMD_INLINE batch<T, A> bitwise_or(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            using register_type = typename batch<T, A>::register_type;
            return detail::bitwise_or_neon(register_type(lhs), register_type(rhs));
        }

        template <class A, class T, detail::enable_neon_type_t<T> = 0>
        XSIMD_INLINE batch_bool<T, A> bitwise_or(batch_bool<T, A> const& lhs, batch_bool<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            using register_type = typename batch_bool<T, A>::register_type;
            return detail::bitwise_or_neon(register_type(lhs), register_type(rhs));
        }

        /***************
         * bitwise_xor *
         ***************/

        WRAP_BINARY_INT(veorq, detail::identity_return_type)

        namespace detail
        {
            XSIMD_INLINE float32x4_t bitwise_xor_f32(float32x4_t lhs, float32x4_t rhs) noexcept
            {
                return vreinterpretq_f32_u32(veorq_u32(vreinterpretq_u32_f32(lhs),
                                                       vreinterpretq_u32_f32(rhs)));
            }

            template <class V>
            XSIMD_INLINE V bitwise_xor_neon(V const& lhs, V const& rhs) noexcept
            {
                const neon_dispatcher::binary dispatcher = {
                    std::make_tuple(wrap::veorq_u8, wrap::veorq_s8, wrap::veorq_u16, wrap::veorq_s16,
                                    wrap::veorq_u32, wrap::veorq_s32, wrap::veorq_u64, wrap::veorq_s64,
                                    bitwise_xor_f32)
                };
                return dispatcher.apply(lhs, rhs);
            }
        }

        template <class A, class T, detail::enable_neon_type_t<T> = 0>
        XSIMD_INLINE batch<T, A> bitwise_xor(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            using register_type = typename batch<T, A>::register_type;
            return detail::bitwise_xor_neon(register_type(lhs), register_type(rhs));
        }

        template <class A, class T, detail::enable_neon_type_t<T> = 0>
        XSIMD_INLINE batch_bool<T, A> bitwise_xor(batch_bool<T, A> const& lhs, batch_bool<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            using register_type = typename batch_bool<T, A>::register_type;
            return detail::bitwise_xor_neon(register_type(lhs), register_type(rhs));
        }

        /*******
         * neq *
         *******/

        template <class A, class T>
        XSIMD_INLINE batch_bool<T, A> neq(batch_bool<T, A> const& lhs, batch_bool<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            return bitwise_xor(lhs, rhs, A {});
        }

        /***************
         * bitwise_not *
         ***************/

        WRAP_UNARY_INT_EXCLUDING_64(vmvnq)

        namespace detail
        {
            XSIMD_INLINE int64x2_t bitwise_not_s64(int64x2_t arg) noexcept
            {
                return vreinterpretq_s64_s32(vmvnq_s32(vreinterpretq_s32_s64(arg)));
            }

            XSIMD_INLINE uint64x2_t bitwise_not_u64(uint64x2_t arg) noexcept
            {
                return vreinterpretq_u64_u32(vmvnq_u32(vreinterpretq_u32_u64(arg)));
            }

            XSIMD_INLINE float32x4_t bitwise_not_f32(float32x4_t arg) noexcept
            {
                return vreinterpretq_f32_u32(vmvnq_u32(vreinterpretq_u32_f32(arg)));
            }

            template <class V>
            XSIMD_INLINE V bitwise_not_neon(V const& arg) noexcept
            {
                const neon_dispatcher::unary dispatcher = {
                    std::make_tuple(wrap::vmvnq_u8, wrap::vmvnq_s8, wrap::vmvnq_u16, wrap::vmvnq_s16,
                                    wrap::vmvnq_u32, wrap::vmvnq_s32,
                                    bitwise_not_u64, bitwise_not_s64,
                                    bitwise_not_f32)
                };
                return dispatcher.apply(arg);
            }
        }

        template <class A, class T, detail::enable_neon_type_t<T> = 0>
        XSIMD_INLINE batch<T, A> bitwise_not(batch<T, A> const& arg, requires_arch<neon>) noexcept
        {
            using register_type = typename batch<T, A>::register_type;
            return detail::bitwise_not_neon(register_type(arg));
        }

        template <class A, class T, detail::enable_neon_type_t<T> = 0>
        XSIMD_INLINE batch_bool<T, A> bitwise_not(batch_bool<T, A> const& arg, requires_arch<neon>) noexcept
        {
            using register_type = typename batch_bool<T, A>::register_type;
            return detail::bitwise_not_neon(register_type(arg));
        }

        /******************
         * bitwise_andnot *
         ******************/

        WRAP_BINARY_INT(vbicq, detail::identity_return_type)

        namespace detail
        {
            XSIMD_INLINE float32x4_t bitwise_andnot_f32(float32x4_t lhs, float32x4_t rhs) noexcept
            {
                return vreinterpretq_f32_u32(vbicq_u32(vreinterpretq_u32_f32(lhs), vreinterpretq_u32_f32(rhs)));
            }

            template <class V>
            XSIMD_INLINE V bitwise_andnot_neon(V const& lhs, V const& rhs) noexcept
            {
                const detail::neon_dispatcher::binary dispatcher = {
                    std::make_tuple(wrap::vbicq_u8, wrap::vbicq_s8, wrap::vbicq_u16, wrap::vbicq_s16,
                                    wrap::vbicq_u32, wrap::vbicq_s32, wrap::vbicq_u64, wrap::vbicq_s64,
                                    bitwise_andnot_f32)
                };
                return dispatcher.apply(lhs, rhs);
            }
        }

        template <class A, class T, detail::enable_neon_type_t<T> = 0>
        XSIMD_INLINE batch<T, A> bitwise_andnot(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            using register_type = typename batch<T, A>::register_type;
            return detail::bitwise_andnot_neon(register_type(lhs), register_type(rhs));
        }

        template <class A, class T, detail::enable_neon_type_t<T> = 0>
        XSIMD_INLINE batch_bool<T, A> bitwise_andnot(batch_bool<T, A> const& lhs, batch_bool<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            using register_type = typename batch_bool<T, A>::register_type;
            return detail::bitwise_andnot_neon(register_type(lhs), register_type(rhs));
        }

        /*******
         * min *
         *******/

        WRAP_BINARY_INT_EXCLUDING_64(vminq, detail::identity_return_type)
        WRAP_BINARY_FLOAT(vminq, detail::identity_return_type)

        template <class A, class T, detail::exclude_int64_neon_t<T> = 0>
        XSIMD_INLINE batch<T, A> min(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            using register_type = typename batch<T, A>::register_type;
            const detail::excluding_int64_dispatcher::binary dispatcher = {
                std::make_tuple(wrap::vminq_u8, wrap::vminq_s8, wrap::vminq_u16, wrap::vminq_s16,
                                wrap::vminq_u32, wrap::vminq_s32, wrap::vminq_f32)
            };
            return dispatcher.apply(register_type(lhs), register_type(rhs));
        }

        template <class A, class T, detail::enable_sized_integral_t<T, 8> = 0>
        XSIMD_INLINE batch<T, A> min(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            return { std::min(lhs.get(0), rhs.get(0)), std::min(lhs.get(1), rhs.get(1)) };
        }

        /*******
         * max *
         *******/

        WRAP_BINARY_INT_EXCLUDING_64(vmaxq, detail::identity_return_type)
        WRAP_BINARY_FLOAT(vmaxq, detail::identity_return_type)

        template <class A, class T, detail::exclude_int64_neon_t<T> = 0>
        XSIMD_INLINE batch<T, A> max(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            using register_type = typename batch<T, A>::register_type;
            const detail::excluding_int64_dispatcher::binary dispatcher = {
                std::make_tuple(wrap::vmaxq_u8, wrap::vmaxq_s8, wrap::vmaxq_u16, wrap::vmaxq_s16,
                                wrap::vmaxq_u32, wrap::vmaxq_s32, wrap::vmaxq_f32)
            };
            return dispatcher.apply(register_type(lhs), register_type(rhs));
        }

        template <class A, class T, detail::enable_sized_integral_t<T, 8> = 0>
        XSIMD_INLINE batch<T, A> max(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            return { std::max(lhs.get(0), rhs.get(0)), std::max(lhs.get(1), rhs.get(1)) };
        }

        /*******
         * abs *
         *******/

        namespace wrap
        {
            XSIMD_INLINE int8x16_t vabsq_s8(int8x16_t a) noexcept { return ::vabsq_s8(a); }
            XSIMD_INLINE int16x8_t vabsq_s16(int16x8_t a) noexcept { return ::vabsq_s16(a); }
            XSIMD_INLINE int32x4_t vabsq_s32(int32x4_t a) noexcept { return ::vabsq_s32(a); }
        }
        WRAP_UNARY_FLOAT(vabsq)

        namespace detail
        {
            XSIMD_INLINE uint8x16_t abs_u8(uint8x16_t arg) noexcept
            {
                return arg;
            }

            XSIMD_INLINE uint16x8_t abs_u16(uint16x8_t arg) noexcept
            {
                return arg;
            }

            XSIMD_INLINE uint32x4_t abs_u32(uint32x4_t arg) noexcept
            {
                return arg;
            }
        }

        template <class A, class T, detail::exclude_int64_neon_t<T> = 0>
        XSIMD_INLINE batch<T, A> abs(batch<T, A> const& arg, requires_arch<neon>) noexcept
        {
            using register_type = typename batch<T, A>::register_type;
            const detail::excluding_int64_dispatcher::unary dispatcher = {
                std::make_tuple(detail::abs_u8, wrap::vabsq_s8, detail::abs_u16, wrap::vabsq_s16,
                                detail::abs_u32, wrap::vabsq_s32, wrap::vabsq_f32)
            };
            return dispatcher.apply(register_type(arg));
        }

        /********
         * rsqrt *
         ********/

        template <class A>
        XSIMD_INLINE batch<float, A> rsqrt(batch<float, A> const& arg, requires_arch<neon>) noexcept
        {
            return vrsqrteq_f32(arg);
        }

        /********
         * sqrt *
         ********/

        template <class A>
        XSIMD_INLINE batch<float, A> sqrt(batch<float, A> const& arg, requires_arch<neon>) noexcept
        {
            batch<float, A> sqrt_reciprocal = vrsqrteq_f32(arg);
            // one iter
            sqrt_reciprocal = sqrt_reciprocal * batch<float, A>(vrsqrtsq_f32(arg * sqrt_reciprocal, sqrt_reciprocal));
            batch<float, A> sqrt_approx = arg * sqrt_reciprocal * batch<float, A>(vrsqrtsq_f32(arg * sqrt_reciprocal, sqrt_reciprocal));
            batch<float, A> zero(0.f);
            return select(arg == zero, zero, sqrt_approx);
        }

        /********************
         * Fused operations *
         ********************/

#ifdef __ARM_FEATURE_FMA
        template <class A>
        XSIMD_INLINE batch<float, A> fma(batch<float, A> const& x, batch<float, A> const& y, batch<float, A> const& z, requires_arch<neon>) noexcept
        {
            return vfmaq_f32(z, x, y);
        }

        template <class A>
        XSIMD_INLINE batch<float, A> fms(batch<float, A> const& x, batch<float, A> const& y, batch<float, A> const& z, requires_arch<neon>) noexcept
        {
            return vfmaq_f32(-z, x, y);
        }
#endif

        /*********
         * haddp *
         *********/

        template <class A>
        XSIMD_INLINE batch<float, A> haddp(const batch<float, A>* row, requires_arch<neon>) noexcept
        {
            // row = (a,b,c,d)
            float32x2_t tmp1, tmp2, tmp3;
            // tmp1 = (a0 + a2, a1 + a3)
            tmp1 = vpadd_f32(vget_low_f32(row[0]), vget_high_f32(row[0]));
            // tmp2 = (b0 + b2, b1 + b3)
            tmp2 = vpadd_f32(vget_low_f32(row[1]), vget_high_f32(row[1]));
            // tmp1 = (a0..3, b0..3)
            tmp1 = vpadd_f32(tmp1, tmp2);
            // tmp2 = (c0 + c2, c1 + c3)
            tmp2 = vpadd_f32(vget_low_f32(row[2]), vget_high_f32(row[2]));
            // tmp3 = (d0 + d2, d1 + d3)
            tmp3 = vpadd_f32(vget_low_f32(row[3]), vget_high_f32(row[3]));
            // tmp1 = (c0..3, d0..3)
            tmp2 = vpadd_f32(tmp2, tmp3);
            // return = (a0..3, b0..3, c0..3, d0..3)
            return vcombine_f32(tmp1, tmp2);
        }

        /**************
         * reciprocal *
         **************/

        template <class A>
        XSIMD_INLINE batch<float, A>
        reciprocal(const batch<float, A>& x,
                   kernel::requires_arch<neon>) noexcept
        {
            return vrecpeq_f32(x);
        }

        /**********
         * insert *
         **********/

        template <class A, class T, size_t I, detail::enable_sized_unsigned_t<T, 1> = 0>
        XSIMD_INLINE batch<T, A> insert(batch<T, A> const& self, T val, index<I>, requires_arch<neon>) noexcept
        {
            return vsetq_lane_u8(val, self, I);
        }

        template <class A, class T, size_t I, detail::enable_sized_signed_t<T, 1> = 0>
        XSIMD_INLINE batch<T, A> insert(batch<T, A> const& self, T val, index<I>, requires_arch<neon>) noexcept
        {
            return vsetq_lane_s8(val, self, I);
        }

        template <class A, class T, size_t I, detail::enable_sized_unsigned_t<T, 2> = 0>
        XSIMD_INLINE batch<T, A> insert(batch<T, A> const& self, T val, index<I>, requires_arch<neon>) noexcept
        {
            return vsetq_lane_u16(val, self, I);
        }

        template <class A, class T, size_t I, detail::enable_sized_signed_t<T, 2> = 0>
        XSIMD_INLINE batch<int16_t, A> insert(batch<int16_t, A> const& self, int16_t val, index<I>, requires_arch<neon>) noexcept
        {
            return vsetq_lane_s16(val, self, I);
        }

        template <class A, class T, size_t I, detail::enable_sized_unsigned_t<T, 4> = 0>
        XSIMD_INLINE batch<T, A> insert(batch<T, A> const& self, T val, index<I>, requires_arch<neon>) noexcept
        {
            return vsetq_lane_u32(val, self, I);
        }

        template <class A, class T, size_t I, detail::enable_sized_signed_t<T, 4> = 0>
        XSIMD_INLINE batch<T, A> insert(batch<T, A> const& self, T val, index<I>, requires_arch<neon>) noexcept
        {
            return vsetq_lane_s32(val, self, I);
        }

        template <class A, class T, size_t I, detail::enable_sized_unsigned_t<T, 8> = 0>
        XSIMD_INLINE batch<T, A> insert(batch<T, A> const& self, T val, index<I>, requires_arch<neon>) noexcept
        {
            return vsetq_lane_u64(val, self, I);
        }

        template <class A, class T, size_t I, detail::enable_sized_signed_t<T, 8> = 0>
        XSIMD_INLINE batch<T, A> insert(batch<T, A> const& self, T val, index<I>, requires_arch<neon>) noexcept
        {
            return vsetq_lane_s64(val, self, I);
        }

        template <class A, size_t I>
        XSIMD_INLINE batch<float, A> insert(batch<float, A> const& self, float val, index<I>, requires_arch<neon>) noexcept
        {
            return vsetq_lane_f32(val, self, I);
        }

        /********************
         * nearbyint_as_int *
         *******************/

        template <class A>
        XSIMD_INLINE batch<int32_t, A> nearbyint_as_int(batch<float, A> const& self,
                                                        requires_arch<neon>) noexcept
        {
            /* origin: https://github.com/DLTcollab/sse2neon/blob/cad518a93b326f0f644b7972d488d04eaa2b0475/sse2neon.h#L4028-L4047 */
            //  Contributors to this work are:
            //   John W. Ratcliff <jratcliffscarab@gmail.com>
            //   Brandon Rowlett <browlett@nvidia.com>
            //   Ken Fast <kfast@gdeb.com>
            //   Eric van Beurden <evanbeurden@nvidia.com>
            //   Alexander Potylitsin <apotylitsin@nvidia.com>
            //   Hasindu Gamaarachchi <hasindu2008@gmail.com>
            //   Jim Huang <jserv@biilabs.io>
            //   Mark Cheng <marktwtn@biilabs.io>
            //   Malcolm James MacLeod <malcolm@gulden.com>
            //   Devin Hussey (easyaspi314) <husseydevin@gmail.com>
            //   Sebastian Pop <spop@amazon.com>
            //   Developer Ecosystem Engineering <DeveloperEcosystemEngineering@apple.com>
            //   Danila Kutenin <danilak@google.com>
            //   François Turban (JishinMaster) <francois.turban@gmail.com>
            //   Pei-Hsuan Hung <afcidk@gmail.com>
            //   Yang-Hao Yuan <yanghau@biilabs.io>
            //   Syoyo Fujita <syoyo@lighttransport.com>
            //   Brecht Van Lommel <brecht@blender.org>

            /*
             * sse2neon is freely redistributable under the MIT License.
             *
             * Permission is hereby granted, free of charge, to any person obtaining a copy
             * of this software and associated documentation files (the "Software"), to deal
             * in the Software without restriction, including without limitation the rights
             * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
             * copies of the Software, and to permit persons to whom the Software is
             * furnished to do so, subject to the following conditions:
             *
             * The above copyright notice and this permission notice shall be included in
             * all copies or substantial portions of the Software.
             *
             * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
             * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
             * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
             * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
             * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
             * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
             * SOFTWARE.
             */

            const auto signmask = vdupq_n_u32(0x80000000);
            const auto half = vbslq_f32(signmask, self,
                                        vdupq_n_f32(0.5f)); /* +/- 0.5 */
            const auto r_normal = vcvtq_s32_f32(vaddq_f32(
                self, half)); /* round to integer: [a + 0.5]*/
            const auto r_trunc = vcvtq_s32_f32(self); /* truncate to integer: [a] */
            const auto plusone = vreinterpretq_s32_u32(vshrq_n_u32(
                vreinterpretq_u32_s32(vnegq_s32(r_trunc)), 31)); /* 1 or 0 */
            const auto r_even = vbicq_s32(vaddq_s32(r_trunc, plusone),
                                          vdupq_n_s32(1)); /* ([a] + {0,1}) & ~1 */
            const auto delta = vsubq_f32(
                self,
                vcvtq_f32_s32(r_trunc)); /* compute delta: delta = (a - [a]) */
            const auto is_delta_half = vceqq_f32(delta, half); /* delta == +/- 0.5 */
            return vbslq_s32(is_delta_half, r_even, r_normal);
        }

        /**************
         * reduce_add *
         **************/

        namespace detail
        {
            template <class T, class A, class V>
            XSIMD_INLINE T sum_batch(V const& arg) noexcept
            {
                T res = T(0);
                for (std::size_t i = 0; i < batch<T, A>::size; ++i)
                {
                    res += arg[i];
                }
                return res;
            }
        }

        template <class A, class T, detail::enable_sized_unsigned_t<T, 1> = 0>
        XSIMD_INLINE typename batch<T, A>::value_type reduce_add(batch<T, A> const& arg, requires_arch<neon>) noexcept
        {
            uint8x8_t tmp = vpadd_u8(vget_low_u8(arg), vget_high_u8(arg));
            tmp = vpadd_u8(tmp, tmp);
            tmp = vpadd_u8(tmp, tmp);
            tmp = vpadd_u8(tmp, tmp);
            return vget_lane_u8(tmp, 0);
        }

        template <class A, class T, detail::enable_sized_signed_t<T, 1> = 0>
        XSIMD_INLINE typename batch<T, A>::value_type reduce_add(batch<T, A> const& arg, requires_arch<neon>) noexcept
        {
            int8x8_t tmp = vpadd_s8(vget_low_s8(arg), vget_high_s8(arg));
            tmp = vpadd_s8(tmp, tmp);
            tmp = vpadd_s8(tmp, tmp);
            tmp = vpadd_s8(tmp, tmp);
            return vget_lane_s8(tmp, 0);
        }

        template <class A, class T, detail::enable_sized_unsigned_t<T, 2> = 0>
        XSIMD_INLINE typename batch<T, A>::value_type reduce_add(batch<T, A> const& arg, requires_arch<neon>) noexcept
        {
            uint16x4_t tmp = vpadd_u16(vget_low_u16(arg), vget_high_u16(arg));
            tmp = vpadd_u16(tmp, tmp);
            tmp = vpadd_u16(tmp, tmp);
            return vget_lane_u16(tmp, 0);
        }

        template <class A, class T, detail::enable_sized_signed_t<T, 2> = 0>
        XSIMD_INLINE typename batch<T, A>::value_type reduce_add(batch<T, A> const& arg, requires_arch<neon>) noexcept
        {
            int16x4_t tmp = vpadd_s16(vget_low_s16(arg), vget_high_s16(arg));
            tmp = vpadd_s16(tmp, tmp);
            tmp = vpadd_s16(tmp, tmp);
            return vget_lane_s16(tmp, 0);
        }

        template <class A, class T, detail::enable_sized_unsigned_t<T, 4> = 0>
        XSIMD_INLINE typename batch<T, A>::value_type reduce_add(batch<T, A> const& arg, requires_arch<neon>) noexcept
        {
            uint32x2_t tmp = vpadd_u32(vget_low_u32(arg), vget_high_u32(arg));
            tmp = vpadd_u32(tmp, tmp);
            return vget_lane_u32(tmp, 0);
        }

        template <class A, class T, detail::enable_sized_signed_t<T, 4> = 0>
        XSIMD_INLINE typename batch<T, A>::value_type reduce_add(batch<T, A> const& arg, requires_arch<neon>) noexcept
        {
            int32x2_t tmp = vpadd_s32(vget_low_s32(arg), vget_high_s32(arg));
            tmp = vpadd_s32(tmp, tmp);
            return vget_lane_s32(tmp, 0);
        }

        template <class A, class T, detail::enable_sized_integral_t<T, 8> = 0>
        XSIMD_INLINE typename batch<T, A>::value_type reduce_add(batch<T, A> const& arg, requires_arch<neon>) noexcept
        {
            return arg.get(0) + arg.get(1);
        }

        template <class A>
        XSIMD_INLINE float reduce_add(batch<float, A> const& arg, requires_arch<neon>) noexcept
        {
            float32x2_t tmp = vpadd_f32(vget_low_f32(arg), vget_high_f32(arg));
            tmp = vpadd_f32(tmp, tmp);
            return vget_lane_f32(tmp, 0);
        }

        /**************
         * reduce_max *
         **************/

        // Using generic implementation because ARM doe snot provide intrinsics
        // for this operation

        /**************
         * reduce_min *
         **************/

        // Using generic implementation because ARM doe snot provide intrinsics
        // for this operation

        /**********
         * select *
         **********/

        namespace wrap
        {
            XSIMD_INLINE uint8x16_t vbslq_u8(uint8x16_t a, uint8x16_t b, uint8x16_t c) noexcept { return ::vbslq_u8(a, b, c); }
            XSIMD_INLINE int8x16_t vbslq_s8(uint8x16_t a, int8x16_t b, int8x16_t c) noexcept { return ::vbslq_s8(a, b, c); }
            XSIMD_INLINE uint16x8_t vbslq_u16(uint16x8_t a, uint16x8_t b, uint16x8_t c) noexcept { return ::vbslq_u16(a, b, c); }
            XSIMD_INLINE int16x8_t vbslq_s16(uint16x8_t a, int16x8_t b, int16x8_t c) noexcept { return ::vbslq_s16(a, b, c); }
            XSIMD_INLINE uint32x4_t vbslq_u32(uint32x4_t a, uint32x4_t b, uint32x4_t c) noexcept { return ::vbslq_u32(a, b, c); }
            XSIMD_INLINE int32x4_t vbslq_s32(uint32x4_t a, int32x4_t b, int32x4_t c) noexcept { return ::vbslq_s32(a, b, c); }
            XSIMD_INLINE uint64x2_t vbslq_u64(uint64x2_t a, uint64x2_t b, uint64x2_t c) noexcept { return ::vbslq_u64(a, b, c); }
            XSIMD_INLINE int64x2_t vbslq_s64(uint64x2_t a, int64x2_t b, int64x2_t c) noexcept { return ::vbslq_s64(a, b, c); }
            XSIMD_INLINE float32x4_t vbslq_f32(uint32x4_t a, float32x4_t b, float32x4_t c) noexcept { return ::vbslq_f32(a, b, c); }
        }

        namespace detail
        {
            template <class... T>
            struct neon_select_dispatcher_impl
            {
                using container_type = std::tuple<T (*)(comp_return_type<T>, T, T)...>;
                const container_type m_func;

                template <class U>
                U apply(comp_return_type<U> cond, U lhs, U rhs) const noexcept
                {
                    using func_type = U (*)(comp_return_type<U>, U, U);
                    auto func = xsimd::detail::get<func_type>(m_func);
                    return func(cond, lhs, rhs);
                }
            };

            using neon_select_dispatcher = neon_select_dispatcher_impl<uint8x16_t, int8x16_t,
                                                                       uint16x8_t, int16x8_t,
                                                                       uint32x4_t, int32x4_t,
                                                                       uint64x2_t, int64x2_t,
                                                                       float32x4_t>;
        }

        template <class A, class T, detail::enable_neon_type_t<T> = 0>
        XSIMD_INLINE batch<T, A> select(batch_bool<T, A> const& cond, batch<T, A> const& a, batch<T, A> const& b, requires_arch<neon>) noexcept
        {
            using bool_register_type = typename batch_bool<T, A>::register_type;
            using register_type = typename batch<T, A>::register_type;
            const detail::neon_select_dispatcher dispatcher = {
                std::make_tuple(wrap::vbslq_u8, wrap::vbslq_s8, wrap::vbslq_u16, wrap::vbslq_s16,
                                wrap::vbslq_u32, wrap::vbslq_s32, wrap::vbslq_u64, wrap::vbslq_s64,
                                wrap::vbslq_f32)
            };
            return dispatcher.apply(bool_register_type(cond), register_type(a), register_type(b));
        }

        template <class A, class T, bool... b, detail::enable_neon_type_t<T> = 0>
        XSIMD_INLINE batch<T, A> select(batch_bool_constant<T, A, b...> const&, batch<T, A> const& true_br, batch<T, A> const& false_br, requires_arch<neon>) noexcept
        {
            return select(batch_bool<T, A> { b... }, true_br, false_br, neon {});
        }

        /**********
         * zip_lo *
         **********/

        template <class A, class T, detail::enable_sized_unsigned_t<T, 1> = 0>
        XSIMD_INLINE batch<T, A> zip_lo(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            uint8x8x2_t tmp = vzip_u8(vget_low_u8(lhs), vget_low_u8(rhs));
            return vcombine_u8(tmp.val[0], tmp.val[1]);
        }

        template <class A, class T, detail::enable_sized_signed_t<T, 1> = 0>
        XSIMD_INLINE batch<T, A> zip_lo(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            int8x8x2_t tmp = vzip_s8(vget_low_s8(lhs), vget_low_s8(rhs));
            return vcombine_s8(tmp.val[0], tmp.val[1]);
        }

        template <class A, class T, detail::enable_sized_unsigned_t<T, 2> = 0>
        XSIMD_INLINE batch<T, A> zip_lo(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            uint16x4x2_t tmp = vzip_u16(vget_low_u16(lhs), vget_low_u16(rhs));
            return vcombine_u16(tmp.val[0], tmp.val[1]);
        }

        template <class A, class T, detail::enable_sized_signed_t<T, 2> = 0>
        XSIMD_INLINE batch<T, A> zip_lo(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            int16x4x2_t tmp = vzip_s16(vget_low_s16(lhs), vget_low_s16(rhs));
            return vcombine_s16(tmp.val[0], tmp.val[1]);
        }

        template <class A, class T, detail::enable_sized_unsigned_t<T, 4> = 0>
        XSIMD_INLINE batch<T, A> zip_lo(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            uint32x2x2_t tmp = vzip_u32(vget_low_u32(lhs), vget_low_u32(rhs));
            return vcombine_u32(tmp.val[0], tmp.val[1]);
        }

        template <class A, class T, detail::enable_sized_signed_t<T, 4> = 0>
        XSIMD_INLINE batch<T, A> zip_lo(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            int32x2x2_t tmp = vzip_s32(vget_low_s32(lhs), vget_low_s32(rhs));
            return vcombine_s32(tmp.val[0], tmp.val[1]);
        }

        template <class A, class T, detail::enable_sized_unsigned_t<T, 8> = 0>
        XSIMD_INLINE batch<T, A> zip_lo(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            return vcombine_u64(vget_low_u64(lhs), vget_low_u64(rhs));
        }

        template <class A, class T, detail::enable_sized_signed_t<T, 8> = 0>
        XSIMD_INLINE batch<T, A> zip_lo(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            return vcombine_s64(vget_low_s64(lhs), vget_low_s64(rhs));
        }

        template <class A>
        XSIMD_INLINE batch<float, A> zip_lo(batch<float, A> const& lhs, batch<float, A> const& rhs, requires_arch<neon>) noexcept
        {
            float32x2x2_t tmp = vzip_f32(vget_low_f32(lhs), vget_low_f32(rhs));
            return vcombine_f32(tmp.val[0], tmp.val[1]);
        }

        /**********
         * zip_hi *
         **********/

        template <class A, class T, detail::enable_sized_unsigned_t<T, 1> = 0>
        XSIMD_INLINE batch<T, A> zip_hi(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            uint8x8x2_t tmp = vzip_u8(vget_high_u8(lhs), vget_high_u8(rhs));
            return vcombine_u8(tmp.val[0], tmp.val[1]);
        }

        template <class A, class T, detail::enable_sized_signed_t<T, 1> = 0>
        XSIMD_INLINE batch<T, A> zip_hi(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            int8x8x2_t tmp = vzip_s8(vget_high_s8(lhs), vget_high_s8(rhs));
            return vcombine_s8(tmp.val[0], tmp.val[1]);
        }

        template <class A, class T, detail::enable_sized_unsigned_t<T, 2> = 0>
        XSIMD_INLINE batch<T, A> zip_hi(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            uint16x4x2_t tmp = vzip_u16(vget_high_u16(lhs), vget_high_u16(rhs));
            return vcombine_u16(tmp.val[0], tmp.val[1]);
        }

        template <class A, class T, detail::enable_sized_signed_t<T, 2> = 0>
        XSIMD_INLINE batch<T, A> zip_hi(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            int16x4x2_t tmp = vzip_s16(vget_high_s16(lhs), vget_high_s16(rhs));
            return vcombine_s16(tmp.val[0], tmp.val[1]);
        }

        template <class A, class T, detail::enable_sized_unsigned_t<T, 4> = 0>
        XSIMD_INLINE batch<T, A> zip_hi(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            uint32x2x2_t tmp = vzip_u32(vget_high_u32(lhs), vget_high_u32(rhs));
            return vcombine_u32(tmp.val[0], tmp.val[1]);
        }

        template <class A, class T, detail::enable_sized_signed_t<T, 4> = 0>
        XSIMD_INLINE batch<T, A> zip_hi(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            int32x2x2_t tmp = vzip_s32(vget_high_s32(lhs), vget_high_s32(rhs));
            return vcombine_s32(tmp.val[0], tmp.val[1]);
        }

        template <class A, class T, detail::enable_sized_unsigned_t<T, 8> = 0>
        XSIMD_INLINE batch<T, A> zip_hi(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            return vcombine_u64(vget_high_u64(lhs), vget_high_u64(rhs));
        }

        template <class A, class T, detail::enable_sized_signed_t<T, 8> = 0>
        XSIMD_INLINE batch<T, A> zip_hi(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            return vcombine_s64(vget_high_s64(lhs), vget_high_s64(rhs));
        }

        template <class A>
        XSIMD_INLINE batch<float, A> zip_hi(batch<float, A> const& lhs, batch<float, A> const& rhs, requires_arch<neon>) noexcept
        {
            float32x2x2_t tmp = vzip_f32(vget_high_f32(lhs), vget_high_f32(rhs));
            return vcombine_f32(tmp.val[0], tmp.val[1]);
        }

        /****************
         * extract_pair *
         ****************/

        namespace detail
        {
            template <class A, class T>
            XSIMD_INLINE batch<T, A> extract_pair(batch<T, A> const&, batch<T, A> const& /*rhs*/, std::size_t, ::xsimd::detail::index_sequence<>) noexcept
            {
                assert(false && "extract_pair out of bounds");
                return batch<T, A> {};
            }

            template <class A, class T, size_t I, size_t... Is, detail::enable_sized_unsigned_t<T, 1> = 0>
            XSIMD_INLINE batch<T, A> extract_pair(batch<T, A> const& lhs, batch<T, A> const& rhs, std::size_t n, ::xsimd::detail::index_sequence<I, Is...>) noexcept
            {
                if (n == I)
                {
                    return vextq_u8(rhs, lhs, I);
                }
                else
                {
                    return extract_pair(lhs, rhs, n, ::xsimd::detail::index_sequence<Is...>());
                }
            }

            template <class A, class T, size_t I, size_t... Is, detail::enable_sized_signed_t<T, 1> = 0>
            XSIMD_INLINE batch<T, A> extract_pair(batch<T, A> const& lhs, batch<T, A> const& rhs, std::size_t n, ::xsimd::detail::index_sequence<I, Is...>) noexcept
            {
                if (n == I)
                {
                    return vextq_s8(rhs, lhs, I);
                }
                else
                {
                    return extract_pair(lhs, rhs, n, ::xsimd::detail::index_sequence<Is...>());
                }
            }

            template <class A, class T, size_t I, size_t... Is, detail::enable_sized_unsigned_t<T, 2> = 0>
            XSIMD_INLINE batch<T, A> extract_pair(batch<T, A> const& lhs, batch<T, A> const& rhs, std::size_t n, ::xsimd::detail::index_sequence<I, Is...>) noexcept
            {
                if (n == I)
                {
                    return vextq_u16(rhs, lhs, I);
                }
                else
                {
                    return extract_pair(lhs, rhs, n, ::xsimd::detail::index_sequence<Is...>());
                }
            }

            template <class A, class T, size_t I, size_t... Is, detail::enable_sized_signed_t<T, 2> = 0>
            XSIMD_INLINE batch<T, A> extract_pair(batch<T, A> const& lhs, batch<T, A> const& rhs, std::size_t n, ::xsimd::detail::index_sequence<I, Is...>) noexcept
            {
                if (n == I)
                {
                    return vextq_s16(rhs, lhs, I);
                }
                else
                {
                    return extract_pair(lhs, rhs, n, ::xsimd::detail::index_sequence<Is...>());
                }
            }

            template <class A, class T, size_t I, size_t... Is, detail::enable_sized_unsigned_t<T, 4> = 0>
            XSIMD_INLINE batch<T, A> extract_pair(batch<T, A> const& lhs, batch<T, A> const& rhs, std::size_t n, ::xsimd::detail::index_sequence<I, Is...>) noexcept
            {
                if (n == I)
                {
                    return vextq_u32(rhs, lhs, I);
                }
                else
                {
                    return extract_pair(lhs, rhs, n, ::xsimd::detail::index_sequence<Is...>());
                }
            }

            template <class A, class T, size_t I, size_t... Is, detail::enable_sized_signed_t<T, 4> = 0>
            XSIMD_INLINE batch<T, A> extract_pair(batch<T, A> const& lhs, batch<T, A> const& rhs, std::size_t n, ::xsimd::detail::index_sequence<I, Is...>) noexcept
            {
                if (n == I)
                {
                    return vextq_s32(rhs, lhs, I);
                }
                else
                {
                    return extract_pair(lhs, rhs, n, ::xsimd::detail::index_sequence<Is...>());
                }
            }

            template <class A, class T, size_t I, size_t... Is, detail::enable_sized_unsigned_t<T, 8> = 0>
            XSIMD_INLINE batch<T, A> extract_pair(batch<T, A> const& lhs, batch<T, A> const& rhs, std::size_t n, ::xsimd::detail::index_sequence<I, Is...>) noexcept
            {
                if (n == I)
                {
                    return vextq_u64(rhs, lhs, I);
                }
                else
                {
                    return extract_pair(lhs, rhs, n, ::xsimd::detail::index_sequence<Is...>());
                }
            }

            template <class A, class T, size_t I, size_t... Is, detail::enable_sized_signed_t<T, 8> = 0>
            XSIMD_INLINE batch<T, A> extract_pair(batch<T, A> const& lhs, batch<T, A> const& rhs, std::size_t n, ::xsimd::detail::index_sequence<I, Is...>) noexcept
            {
                if (n == I)
                {
                    return vextq_s64(rhs, lhs, I);
                }
                else
                {
                    return extract_pair(lhs, rhs, n, ::xsimd::detail::index_sequence<Is...>());
                }
            }

            template <class A, size_t I, size_t... Is>
            XSIMD_INLINE batch<float, A> extract_pair(batch<float, A> const& lhs, batch<float, A> const& rhs, std::size_t n, ::xsimd::detail::index_sequence<I, Is...>) noexcept
            {
                if (n == I)
                {
                    return vextq_f32(rhs, lhs, I);
                }
                else
                {
                    return extract_pair(lhs, rhs, n, ::xsimd::detail::index_sequence<Is...>());
                }
            }

            template <class A, class T, size_t... Is>
            XSIMD_INLINE batch<T, A> extract_pair_impl(batch<T, A> const& lhs, batch<T, A> const& rhs, std::size_t n, ::xsimd::detail::index_sequence<0, Is...>) noexcept
            {
                if (n == 0)
                {
                    return rhs;
                }
                else
                {
                    return extract_pair(lhs, rhs, n, ::xsimd::detail::index_sequence<Is...>());
                }
            }
        }

        template <class A, class T>
        XSIMD_INLINE batch<T, A> extract_pair(batch<T, A> const& lhs, batch<T, A> const& rhs, std::size_t n, requires_arch<neon>) noexcept
        {
            constexpr std::size_t size = batch<T, A>::size;
            assert(n < size && "index in bounds");
            return detail::extract_pair_impl(lhs, rhs, n, ::xsimd::detail::make_index_sequence<size>());
        }

        /******************
         * bitwise_lshift *
         ******************/

        namespace detail
        {
            template <class A, class T>
            XSIMD_INLINE batch<T, A> bitwise_lshift(batch<T, A> const& /*lhs*/, int /*n*/, ::xsimd::detail::int_sequence<>) noexcept
            {
                assert(false && "bitwise_lshift out of bounds");
                return batch<T, A> {};
            }

            template <class A, class T, int I, int... Is, detail::enable_sized_unsigned_t<T, 1> = 0>
            XSIMD_INLINE batch<T, A> bitwise_lshift(batch<T, A> const& lhs, int n, ::xsimd::detail::int_sequence<I, Is...>) noexcept
            {
                if (n == I)
                {
                    return vshlq_n_u8(lhs, I);
                }
                else
                {
                    return bitwise_lshift(lhs, n, ::xsimd::detail::int_sequence<Is...>());
                }
            }

            template <class A, class T, int I, int... Is, detail::enable_sized_signed_t<T, 1> = 0>
            XSIMD_INLINE batch<T, A> bitwise_lshift(batch<T, A> const& lhs, int n, ::xsimd::detail::int_sequence<I, Is...>) noexcept
            {
                if (n == I)
                {
                    return vshlq_n_s8(lhs, I);
                }
                else
                {
                    return bitwise_lshift(lhs, n, ::xsimd::detail::int_sequence<Is...>());
                }
            }

            template <class A, class T, int I, int... Is, detail::enable_sized_unsigned_t<T, 2> = 0>
            XSIMD_INLINE batch<T, A> bitwise_lshift(batch<T, A> const& lhs, int n, ::xsimd::detail::int_sequence<I, Is...>) noexcept
            {
                if (n == I)
                {
                    return vshlq_n_u16(lhs, I);
                }
                else
                {
                    return bitwise_lshift(lhs, n, ::xsimd::detail::int_sequence<Is...>());
                }
            }

            template <class A, class T, int I, int... Is, detail::enable_sized_signed_t<T, 2> = 0>
            XSIMD_INLINE batch<T, A> bitwise_lshift(batch<T, A> const& lhs, int n, ::xsimd::detail::int_sequence<I, Is...>) noexcept
            {
                if (n == I)
                {
                    return vshlq_n_s16(lhs, I);
                }
                else
                {
                    return bitwise_lshift(lhs, n, ::xsimd::detail::int_sequence<Is...>());
                }
            }

            template <class A, class T, int I, int... Is, detail::enable_sized_unsigned_t<T, 4> = 0>
            XSIMD_INLINE batch<T, A> bitwise_lshift(batch<T, A> const& lhs, int n, ::xsimd::detail::int_sequence<I, Is...>) noexcept
            {
                if (n == I)
                {
                    return vshlq_n_u32(lhs, I);
                }
                else
                {
                    return bitwise_lshift(lhs, n, ::xsimd::detail::int_sequence<Is...>());
                }
            }

            template <class A, class T, int I, int... Is, detail::enable_sized_signed_t<T, 4> = 0>
            XSIMD_INLINE batch<T, A> bitwise_lshift(batch<T, A> const& lhs, int n, ::xsimd::detail::int_sequence<I, Is...>) noexcept
            {
                if (n == I)
                {
                    return vshlq_n_s32(lhs, I);
                }
                else
                {
                    return bitwise_lshift(lhs, n, ::xsimd::detail::int_sequence<Is...>());
                }
            }

            template <class A, class T, int I, int... Is, detail::enable_sized_unsigned_t<T, 8> = 0>
            XSIMD_INLINE batch<T, A> bitwise_lshift(batch<T, A> const& lhs, int n, ::xsimd::detail::int_sequence<I, Is...>) noexcept
            {
                if (n == I)
                {
                    return vshlq_n_u64(lhs, I);
                }
                else
                {
                    return bitwise_lshift(lhs, n, ::xsimd::detail::int_sequence<Is...>());
                }
            }

            template <class A, class T, int I, int... Is, detail::enable_sized_signed_t<T, 8> = 0>
            XSIMD_INLINE batch<T, A> bitwise_lshift(batch<T, A> const& lhs, int n, ::xsimd::detail::int_sequence<I, Is...>) noexcept
            {
                if (n == I)
                {
                    return vshlq_n_s64(lhs, I);
                }
                else
                {
                    return bitwise_lshift(lhs, n, ::xsimd::detail::int_sequence<Is...>());
                }
            }

            template <class A, class T, int... Is>
            XSIMD_INLINE batch<T, A> bitwise_lshift_impl(batch<T, A> const& lhs, int n, ::xsimd::detail::int_sequence<0, Is...>) noexcept
            {
                if (n == 0)
                {
                    return lhs;
                }
                else
                {
                    return bitwise_lshift(lhs, n, ::xsimd::detail::int_sequence<Is...>());
                }
            }
        }

        template <class A, class T>
        XSIMD_INLINE batch<T, A> bitwise_lshift(batch<T, A> const& lhs, int n, requires_arch<neon>) noexcept
        {
            constexpr int size = sizeof(typename batch<T, A>::value_type) * 8;
            assert(0 <= n && n < size && "index in bounds");
            return detail::bitwise_lshift_impl(lhs, n, ::xsimd::detail::make_int_sequence<size>());
        }

        template <class A, class T, detail::enable_sized_unsigned_t<T, 1> = 0>
        XSIMD_INLINE batch<T, A> bitwise_lshift(batch<T, A> const& lhs, batch<as_signed_integer_t<T>, A> const& rhs, requires_arch<neon>) noexcept
        {
            return vshlq_u8(lhs, rhs);
        }

        template <class A, class T, detail::enable_sized_signed_t<T, 1> = 0>
        XSIMD_INLINE batch<T, A> bitwise_lshift(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            return vshlq_s8(lhs, rhs);
        }

        template <class A, class T, detail::enable_sized_unsigned_t<T, 2> = 0>
        XSIMD_INLINE batch<T, A> bitwise_lshift(batch<T, A> const& lhs, batch<as_signed_integer_t<T>, A> const& rhs, requires_arch<neon>) noexcept
        {
            return vshlq_u16(lhs, rhs);
        }

        template <class A, class T, detail::enable_sized_signed_t<T, 2> = 0>
        XSIMD_INLINE batch<T, A> bitwise_lshift(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            return vshlq_s16(lhs, rhs);
        }

        template <class A, class T, detail::enable_sized_unsigned_t<T, 4> = 0>
        XSIMD_INLINE batch<T, A> bitwise_lshift(batch<T, A> const& lhs, batch<as_signed_integer_t<T>, A> const& rhs, requires_arch<neon>) noexcept
        {
            return vshlq_u32(lhs, rhs);
        }

        template <class A, class T, detail::enable_sized_signed_t<T, 4> = 0>
        XSIMD_INLINE batch<T, A> bitwise_lshift(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            return vshlq_s32(lhs, rhs);
        }

        template <class A, class T, detail::enable_sized_unsigned_t<T, 8> = 0>
        XSIMD_INLINE batch<T, A> bitwise_lshift(batch<T, A> const& lhs, batch<as_signed_integer_t<T>, A> const& rhs, requires_arch<neon>) noexcept
        {
            return vshlq_u64(lhs, rhs);
        }

        template <class A, class T, detail::enable_sized_signed_t<T, 8> = 0>
        XSIMD_INLINE batch<T, A> bitwise_lshift(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            return vshlq_s64(lhs, rhs);
        }

        /******************
         * bitwise_rshift *
         ******************/

        namespace detail
        {
            template <class A, class T>
            XSIMD_INLINE batch<T, A> bitwise_rshift(batch<T, A> const& /*lhs*/, int /*n*/, ::xsimd::detail::int_sequence<>) noexcept
            {
                assert(false && "bitwise_rshift out of bounds");
                return batch<T, A> {};
            }

            template <class A, class T, int I, int... Is, detail::enable_sized_unsigned_t<T, 1> = 0>
            XSIMD_INLINE batch<T, A> bitwise_rshift(batch<T, A> const& lhs, int n, ::xsimd::detail::int_sequence<I, Is...>) noexcept
            {
                if (n == I)
                {
                    return vshrq_n_u8(lhs, I);
                }
                else
                {
                    return bitwise_rshift(lhs, n, ::xsimd::detail::int_sequence<Is...>());
                }
            }

            template <class A, class T, int I, int... Is, detail::enable_sized_signed_t<T, 1> = 0>
            XSIMD_INLINE batch<T, A> bitwise_rshift(batch<T, A> const& lhs, int n, ::xsimd::detail::int_sequence<I, Is...>) noexcept
            {
                if (n == I)
                {
                    return vshrq_n_s8(lhs, I);
                }
                else
                {
                    return bitwise_rshift(lhs, n, ::xsimd::detail::int_sequence<Is...>());
                }
            }

            template <class A, class T, int I, int... Is, detail::enable_sized_unsigned_t<T, 2> = 0>
            XSIMD_INLINE batch<T, A> bitwise_rshift(batch<T, A> const& lhs, int n, ::xsimd::detail::int_sequence<I, Is...>) noexcept
            {
                if (n == I)
                {
                    return vshrq_n_u16(lhs, I);
                }
                else
                {
                    return bitwise_rshift(lhs, n, ::xsimd::detail::int_sequence<Is...>());
                }
            }

            template <class A, class T, int I, int... Is, detail::enable_sized_signed_t<T, 2> = 0>
            XSIMD_INLINE batch<T, A> bitwise_rshift(batch<T, A> const& lhs, int n, ::xsimd::detail::int_sequence<I, Is...>) noexcept
            {
                if (n == I)
                {
                    return vshrq_n_s16(lhs, I);
                }
                else
                {
                    return bitwise_rshift(lhs, n, ::xsimd::detail::int_sequence<Is...>());
                }
            }

            template <class A, class T, int I, int... Is, detail::enable_sized_unsigned_t<T, 4> = 0>
            XSIMD_INLINE batch<T, A> bitwise_rshift(batch<T, A> const& lhs, int n, ::xsimd::detail::int_sequence<I, Is...>) noexcept
            {
                if (n == I)
                {
                    return vshrq_n_u32(lhs, I);
                }
                else
                {
                    return bitwise_rshift(lhs, n, ::xsimd::detail::int_sequence<Is...>());
                }
            }

            template <class A, class T, int I, int... Is, detail::enable_sized_signed_t<T, 4> = 0>
            XSIMD_INLINE batch<T, A> bitwise_rshift(batch<T, A> const& lhs, int n, ::xsimd::detail::int_sequence<I, Is...>) noexcept
            {
                if (n == I)
                {
                    return vshrq_n_s32(lhs, I);
                }
                else
                {
                    return bitwise_rshift(lhs, n, ::xsimd::detail::int_sequence<Is...>());
                }
            }

            template <class A, class T, int I, int... Is, detail::enable_sized_unsigned_t<T, 8> = 0>
            XSIMD_INLINE batch<T, A> bitwise_rshift(batch<T, A> const& lhs, int n, ::xsimd::detail::int_sequence<I, Is...>) noexcept
            {
                if (n == I)
                {
                    return vshrq_n_u64(lhs, I);
                }
                else
                {
                    return bitwise_rshift(lhs, n, ::xsimd::detail::int_sequence<Is...>());
                }
            }

            template <class A, class T, int I, int... Is, detail::enable_sized_signed_t<T, 8> = 0>
            XSIMD_INLINE batch<T, A> bitwise_rshift(batch<T, A> const& lhs, int n, ::xsimd::detail::int_sequence<I, Is...>) noexcept
            {
                if (n == I)
                {
                    return vshrq_n_s64(lhs, I);
                }
                else
                {
                    return bitwise_rshift(lhs, n, ::xsimd::detail::int_sequence<Is...>());
                }
            }

            template <class A, class T, int... Is>
            XSIMD_INLINE batch<T, A> bitwise_rshift_impl(batch<T, A> const& lhs, int n, ::xsimd::detail::int_sequence<0, Is...>) noexcept
            {
                if (n == 0)
                {
                    return lhs;
                }
                else
                {
                    return bitwise_rshift(lhs, n, ::xsimd::detail::int_sequence<Is...>());
                }
            }
        }

        template <class A, class T>
        XSIMD_INLINE batch<T, A> bitwise_rshift(batch<T, A> const& lhs, int n, requires_arch<neon>) noexcept
        {
            constexpr int size = sizeof(typename batch<T, A>::value_type) * 8;
            assert(0 <= n && n < size && "index in bounds");
            return detail::bitwise_rshift_impl(lhs, n, ::xsimd::detail::make_int_sequence<size>());
        }

        template <class A, class T, detail::enable_sized_unsigned_t<T, 1> = 0>
        XSIMD_INLINE batch<T, A> bitwise_rshift(batch<T, A> const& lhs, batch<as_signed_integer_t<T>, A> const& rhs, requires_arch<neon>) noexcept
        {
            return vshlq_u8(lhs, vnegq_s8(rhs));
        }

        template <class A, class T, detail::enable_sized_signed_t<T, 1> = 0>
        XSIMD_INLINE batch<T, A> bitwise_rshift(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            return vshlq_s8(lhs, vnegq_s8(rhs));
        }

        template <class A, class T, detail::enable_sized_unsigned_t<T, 2> = 0>
        XSIMD_INLINE batch<T, A> bitwise_rshift(batch<T, A> const& lhs, batch<as_signed_integer_t<T>, A> const& rhs, requires_arch<neon>) noexcept
        {
            return vshlq_u16(lhs, vnegq_s16(rhs));
        }

        template <class A, class T, detail::enable_sized_signed_t<T, 2> = 0>
        XSIMD_INLINE batch<T, A> bitwise_rshift(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            return vshlq_s16(lhs, vnegq_s16(rhs));
        }

        template <class A, class T, detail::enable_sized_unsigned_t<T, 4> = 0>
        XSIMD_INLINE batch<T, A> bitwise_rshift(batch<T, A> const& lhs, batch<as_signed_integer_t<T>, A> const& rhs, requires_arch<neon>) noexcept
        {
            return vshlq_u32(lhs, vnegq_s32(rhs));
        }

        template <class A, class T, detail::enable_sized_signed_t<T, 4> = 0>
        XSIMD_INLINE batch<T, A> bitwise_rshift(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<neon>) noexcept
        {
            return vshlq_s32(lhs, vnegq_s32(rhs));
        }

        // Overloads of bitwise shifts accepting two batches of uint64/int64 are not available with ARMv7

        /*******
         * all *
         *******/

        template <class A, class T, detail::enable_sized_t<T, 8> = 0>
        XSIMD_INLINE bool all(batch_bool<T, A> const& arg, requires_arch<neon>) noexcept
        {
            uint64x1_t tmp = vand_u64(vget_low_u64(arg), vget_high_u64(arg));
            return vget_lane_u64(tmp, 0) == ~0ULL;
        }

        template <class A, class T, detail::enable_sized_t<T, 1> = 0>
        XSIMD_INLINE bool all(batch_bool<T, A> const& arg, requires_arch<neon>) noexcept
        {
            return all(batch_bool<uint64_t, A>(vreinterpretq_u64_u8(arg)), neon {});
        }

        template <class A, class T, detail::enable_sized_t<T, 2> = 0>
        XSIMD_INLINE bool all(batch_bool<T, A> const& arg, requires_arch<neon>) noexcept
        {
            return all(batch_bool<uint64_t, A>(vreinterpretq_u64_u16(arg)), neon {});
        }

        template <class A, class T, detail::enable_sized_t<T, 4> = 0>
        XSIMD_INLINE bool all(batch_bool<T, A> const& arg, requires_arch<neon>) noexcept
        {
            return all(batch_bool<uint64_t, A>(vreinterpretq_u64_u32(arg)), neon {});
        }

        /*******
         * any *
         *******/

        template <class A, class T, detail::enable_sized_t<T, 8> = 0>
        XSIMD_INLINE bool any(batch_bool<T, A> const& arg, requires_arch<neon>) noexcept
        {
            uint32x2_t tmp = vqmovn_u64(arg);
            return vget_lane_u64(vreinterpret_u64_u32(tmp), 0) != 0;
        }

        template <class A, class T, detail::enable_sized_t<T, 1> = 0>
        XSIMD_INLINE bool any(batch_bool<T, A> const& arg, requires_arch<neon>) noexcept
        {
            return any(batch_bool<uint64_t, A>(vreinterpretq_u64_u8(arg)), neon {});
        }

        template <class A, class T, detail::enable_sized_t<T, 2> = 0>
        XSIMD_INLINE bool any(batch_bool<T, A> const& arg, requires_arch<neon>) noexcept
        {
            return any(batch_bool<uint64_t, A>(vreinterpretq_u64_u16(arg)), neon {});
        }

        template <class A, class T, detail::enable_sized_t<T, 4> = 0>
        XSIMD_INLINE bool any(batch_bool<T, A> const& arg, requires_arch<neon>) noexcept
        {
            return any(batch_bool<uint64_t, A>(vreinterpretq_u64_u32(arg)), neon {});
        }

        /****************
         * bitwise_cast *
         ****************/

#define WRAP_CAST(SUFFIX, TYPE)                                                \
    namespace wrap                                                             \
    {                                                                          \
        XSIMD_INLINE TYPE vreinterpretq_##SUFFIX##_u8(uint8x16_t a) noexcept   \
        {                                                                      \
            return ::vreinterpretq_##SUFFIX##_u8(a);                           \
        }                                                                      \
        XSIMD_INLINE TYPE vreinterpretq_##SUFFIX##_s8(int8x16_t a) noexcept    \
        {                                                                      \
            return ::vreinterpretq_##SUFFIX##_s8(a);                           \
        }                                                                      \
        XSIMD_INLINE TYPE vreinterpretq_##SUFFIX##_u16(uint16x8_t a) noexcept  \
        {                                                                      \
            return ::vreinterpretq_##SUFFIX##_u16(a);                          \
        }                                                                      \
        XSIMD_INLINE TYPE vreinterpretq_##SUFFIX##_s16(int16x8_t a) noexcept   \
        {                                                                      \
            return ::vreinterpretq_##SUFFIX##_s16(a);                          \
        }                                                                      \
        XSIMD_INLINE TYPE vreinterpretq_##SUFFIX##_u32(uint32x4_t a) noexcept  \
        {                                                                      \
            return ::vreinterpretq_##SUFFIX##_u32(a);                          \
        }                                                                      \
        XSIMD_INLINE TYPE vreinterpretq_##SUFFIX##_s32(int32x4_t a) noexcept   \
        {                                                                      \
            return ::vreinterpretq_##SUFFIX##_s32(a);                          \
        }                                                                      \
        XSIMD_INLINE TYPE vreinterpretq_##SUFFIX##_u64(uint64x2_t a) noexcept  \
        {                                                                      \
            return ::vreinterpretq_##SUFFIX##_u64(a);                          \
        }                                                                      \
        XSIMD_INLINE TYPE vreinterpretq_##SUFFIX##_s64(int64x2_t a) noexcept   \
        {                                                                      \
            return ::vreinterpretq_##SUFFIX##_s64(a);                          \
        }                                                                      \
        XSIMD_INLINE TYPE vreinterpretq_##SUFFIX##_f32(float32x4_t a) noexcept \
        {                                                                      \
            return ::vreinterpretq_##SUFFIX##_f32(a);                          \
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

        namespace detail
        {
            template <class R, class... T>
            struct bitwise_caster_impl
            {
                using container_type = std::tuple<R (*)(T)...>;
                container_type m_func;

                template <class U>
                R apply(U rhs) const noexcept
                {
                    using func_type = R (*)(U);
                    auto func = xsimd::detail::get<func_type>(m_func);
                    return func(rhs);
                }
            };

            template <class R, class... T>
            XSIMD_INLINE const bitwise_caster_impl<R, T...> make_bitwise_caster_impl(R (*... arg)(T)) noexcept
            {
                return { std::make_tuple(arg...) };
            }

            template <class... T>
            struct type_list
            {
            };

            template <class RTL, class TTL>
            struct bitwise_caster;

            template <class... R, class... T>
            struct bitwise_caster<type_list<R...>, type_list<T...>>
            {
                using container_type = std::tuple<bitwise_caster_impl<R, T...>...>;
                container_type m_caster;

                template <class V, class U>
                V apply(U rhs) const noexcept
                {
                    using caster_type = bitwise_caster_impl<V, T...>;
                    auto caster = xsimd::detail::get<caster_type>(m_caster);
                    return caster.apply(rhs);
                }
            };

            template <class... T>
            using bitwise_caster_t = bitwise_caster<type_list<T...>, type_list<T...>>;

            using neon_bitwise_caster = bitwise_caster_t<uint8x16_t, int8x16_t,
                                                         uint16x8_t, int16x8_t,
                                                         uint32x4_t, int32x4_t,
                                                         uint64x2_t, int64x2_t,
                                                         float32x4_t>;
        }

        template <class A, class T, class R>
        XSIMD_INLINE batch<R, A> bitwise_cast(batch<T, A> const& arg, batch<R, A> const&, requires_arch<neon>) noexcept
        {
            const detail::neon_bitwise_caster caster = {
                std::make_tuple(
                    detail::make_bitwise_caster_impl(wrap::vreinterpretq_u8_u8, wrap::vreinterpretq_u8_s8, wrap::vreinterpretq_u8_u16, wrap::vreinterpretq_u8_s16,
                                                     wrap::vreinterpretq_u8_u32, wrap::vreinterpretq_u8_s32, wrap::vreinterpretq_u8_u64, wrap::vreinterpretq_u8_s64,
                                                     wrap::vreinterpretq_u8_f32),
                    detail::make_bitwise_caster_impl(wrap::vreinterpretq_s8_u8, wrap::vreinterpretq_s8_s8, wrap::vreinterpretq_s8_u16, wrap::vreinterpretq_s8_s16,
                                                     wrap::vreinterpretq_s8_u32, wrap::vreinterpretq_s8_s32, wrap::vreinterpretq_s8_u64, wrap::vreinterpretq_s8_s64,
                                                     wrap::vreinterpretq_s8_f32),
                    detail::make_bitwise_caster_impl(wrap::vreinterpretq_u16_u8, wrap::vreinterpretq_u16_s8, wrap::vreinterpretq_u16_u16, wrap::vreinterpretq_u16_s16,
                                                     wrap::vreinterpretq_u16_u32, wrap::vreinterpretq_u16_s32, wrap::vreinterpretq_u16_u64, wrap::vreinterpretq_u16_s64,
                                                     wrap::vreinterpretq_u16_f32),
                    detail::make_bitwise_caster_impl(wrap::vreinterpretq_s16_u8, wrap::vreinterpretq_s16_s8, wrap::vreinterpretq_s16_u16, wrap::vreinterpretq_s16_s16,
                                                     wrap::vreinterpretq_s16_u32, wrap::vreinterpretq_s16_s32, wrap::vreinterpretq_s16_u64, wrap::vreinterpretq_s16_s64,
                                                     wrap::vreinterpretq_s16_f32),
                    detail::make_bitwise_caster_impl(wrap::vreinterpretq_u32_u8, wrap::vreinterpretq_u32_s8, wrap::vreinterpretq_u32_u16, wrap::vreinterpretq_u32_s16,
                                                     wrap::vreinterpretq_u32_u32, wrap::vreinterpretq_u32_s32, wrap::vreinterpretq_u32_u64, wrap::vreinterpretq_u32_s64,
                                                     wrap::vreinterpretq_u32_f32),
                    detail::make_bitwise_caster_impl(wrap::vreinterpretq_s32_u8, wrap::vreinterpretq_s32_s8, wrap::vreinterpretq_s32_u16, wrap::vreinterpretq_s32_s16,
                                                     wrap::vreinterpretq_s32_u32, wrap::vreinterpretq_s32_s32, wrap::vreinterpretq_s32_u64, wrap::vreinterpretq_s32_s64,
                                                     wrap::vreinterpretq_s32_f32),
                    detail::make_bitwise_caster_impl(wrap::vreinterpretq_u64_u8, wrap::vreinterpretq_u64_s8, wrap::vreinterpretq_u64_u16, wrap::vreinterpretq_u64_s16,
                                                     wrap::vreinterpretq_u64_u32, wrap::vreinterpretq_u64_s32, wrap::vreinterpretq_u64_u64, wrap::vreinterpretq_u64_s64,
                                                     wrap::vreinterpretq_u64_f32),
                    detail::make_bitwise_caster_impl(wrap::vreinterpretq_s64_u8, wrap::vreinterpretq_s64_s8, wrap::vreinterpretq_s64_u16, wrap::vreinterpretq_s64_s16,
                                                     wrap::vreinterpretq_s64_u32, wrap::vreinterpretq_s64_s32, wrap::vreinterpretq_s64_u64, wrap::vreinterpretq_s64_s64,
                                                     wrap::vreinterpretq_s64_f32),
                    detail::make_bitwise_caster_impl(wrap::vreinterpretq_f32_u8, wrap::vreinterpretq_f32_s8, wrap::vreinterpretq_f32_u16, wrap::vreinterpretq_f32_s16,
                                                     wrap::vreinterpretq_f32_u32, wrap::vreinterpretq_f32_s32, wrap::vreinterpretq_f32_u64, wrap::vreinterpretq_f32_s64,
                                                     wrap::vreinterpretq_f32_f32))
            };
            using src_register_type = typename batch<T, A>::register_type;
            using dst_register_type = typename batch<R, A>::register_type;
            return caster.apply<dst_register_type>(src_register_type(arg));
        }

        /*********
         * isnan *
         *********/

        template <class A>
        XSIMD_INLINE batch_bool<float, A> isnan(batch<float, A> const& arg, requires_arch<neon>) noexcept
        {
            return !(arg == arg);
        }

        // slide_left
        namespace detail
        {
            template <size_t N>
            struct slider_left
            {
                template <class A, class T>
                XSIMD_INLINE batch<T, A> operator()(batch<T, A> const& x, requires_arch<neon>) noexcept
                {
                    const auto left = vdupq_n_u8(0);
                    const auto right = bitwise_cast<uint8_t>(x).data;
                    const batch<uint8_t, A> res(vextq_u8(left, right, 16 - N));
                    return bitwise_cast<T>(res);
                }
            };

            template <>
            struct slider_left<0>
            {
                template <class A, class T>
                XSIMD_INLINE batch<T, A> operator()(batch<T, A> const& x, requires_arch<neon>) noexcept
                {
                    return x;
                }
            };
        } // namespace detail

        template <size_t N, class A, class T>
        XSIMD_INLINE batch<T, A> slide_left(batch<T, A> const& x, requires_arch<neon>) noexcept
        {
            return detail::slider_left<N> {}(x, A {});
        }

        // slide_right
        namespace detail
        {
            template <size_t N>
            struct slider_right
            {
                template <class A, class T>
                XSIMD_INLINE batch<T, A> operator()(batch<T, A> const& x, requires_arch<neon>) noexcept
                {
                    const auto left = bitwise_cast<uint8_t>(x).data;
                    const auto right = vdupq_n_u8(0);
                    const batch<uint8_t, A> res(vextq_u8(left, right, N));
                    return bitwise_cast<T>(res);
                }
            };

            template <>
            struct slider_right<16>
            {
                template <class A, class T>
                XSIMD_INLINE batch<T, A> operator()(batch<T, A> const&, requires_arch<neon>) noexcept
                {
                    return batch<T, A> {};
                }
            };
        } // namespace detail

        template <size_t N, class A, class T>
        XSIMD_INLINE batch<T, A> slide_right(batch<T, A> const& x, requires_arch<neon>) noexcept
        {
            return detail::slider_right<N> {}(x, A {});
        }

        /****************
         * rotate_right *
         ****************/
        namespace wrap
        {
            template <size_t N>
            XSIMD_INLINE uint8x16_t rotate_right_u8(uint8x16_t a, uint8x16_t b) noexcept { return vextq_u8(a, b, N); }
            template <size_t N>
            XSIMD_INLINE int8x16_t rotate_right_s8(int8x16_t a, int8x16_t b) noexcept { return vextq_s8(a, b, N); }
            template <size_t N>
            XSIMD_INLINE uint16x8_t rotate_right_u16(uint16x8_t a, uint16x8_t b) noexcept { return vextq_u16(a, b, N); }
            template <size_t N>
            XSIMD_INLINE int16x8_t rotate_right_s16(int16x8_t a, int16x8_t b) noexcept { return vextq_s16(a, b, N); }
            template <size_t N>
            XSIMD_INLINE uint32x4_t rotate_right_u32(uint32x4_t a, uint32x4_t b) noexcept { return vextq_u32(a, b, N); }
            template <size_t N>
            XSIMD_INLINE int32x4_t rotate_right_s32(int32x4_t a, int32x4_t b) noexcept { return vextq_s32(a, b, N); }
            template <size_t N>
            XSIMD_INLINE uint64x2_t rotate_right_u64(uint64x2_t a, uint64x2_t b) noexcept { return vextq_u64(a, b, N); }
            template <size_t N>
            XSIMD_INLINE int64x2_t rotate_right_s64(int64x2_t a, int64x2_t b) noexcept { return vextq_s64(a, b, N); }
            template <size_t N>
            XSIMD_INLINE float32x4_t rotate_right_f32(float32x4_t a, float32x4_t b) noexcept { return vextq_f32(a, b, N); }
        }

        template <size_t N, class A, class T, detail::enable_neon_type_t<T> = 0>
        XSIMD_INLINE batch<T, A> rotate_right(batch<T, A> const& a, requires_arch<neon>) noexcept
        {
            using register_type = typename batch<T, A>::register_type;
            const detail::neon_dispatcher::binary dispatcher = {
                std::make_tuple(wrap::rotate_right_u8<N>, wrap::rotate_right_s8<N>, wrap::rotate_right_u16<N>, wrap::rotate_right_s16<N>,
                                wrap::rotate_right_u32<N>, wrap::rotate_right_s32<N>, wrap::rotate_right_u64<N>, wrap::rotate_right_s64<N>,
                                wrap::rotate_right_f32<N>)
            };
            return dispatcher.apply(register_type(a), register_type(a));
        }
    }

    template <typename T, class A, T... Values>
    struct batch_constant;

    namespace kernel
    {
        /***********
         * swizzle *
         ***********/

        template <class A, class T, class I, I... idx>
        XSIMD_INLINE batch<T, A> swizzle(batch<T, A> const& self,
                                         batch_constant<I, A, idx...>,
                                         requires_arch<neon>) noexcept
        {
            static_assert(batch<T, A>::size == sizeof...(idx), "valid swizzle indices");
            std::array<T, batch<T, A>::size> data;
            self.store_aligned(data.data());
            return set(batch<T, A>(), A(), data[idx]...);
        }
    }
}

#undef WRAP_BINARY_INT_EXCLUDING_64
#undef WRAP_BINARY_INT
#undef WRAP_BINARY_FLOAT
#undef WRAP_UNARY_INT_EXCLUDING_64
#undef WRAP_UNARY_INT
#undef WRAP_UNARY_FLOAT

#endif
