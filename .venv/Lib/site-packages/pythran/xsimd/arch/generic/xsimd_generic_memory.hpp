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

#ifndef XSIMD_GENERIC_MEMORY_HPP
#define XSIMD_GENERIC_MEMORY_HPP

#include <algorithm>
#include <complex>
#include <stdexcept>

#include "../../types/xsimd_batch_constant.hpp"
#include "./xsimd_generic_details.hpp"

namespace xsimd
{
    template <typename T, class A, T... Values>
    struct batch_constant;

    template <typename T, class A, bool... Values>
    struct batch_bool_constant;

    namespace kernel
    {

        using namespace types;

        // compress
        namespace detail
        {
            template <class IT, class A, class I, size_t... Is>
            XSIMD_INLINE batch<IT, A> create_compress_swizzle_mask(I bitmask, ::xsimd::detail::index_sequence<Is...>)
            {
                batch<IT, A> swizzle_mask(IT(0));
                alignas(A::alignment()) IT mask_buffer[batch<IT, A>::size] = { Is... };
                size_t inserted = 0;
                for (size_t i = 0; i < sizeof...(Is); ++i)
                    if ((bitmask >> i) & 1u)
                        std::swap(mask_buffer[inserted++], mask_buffer[i]);
                return batch<IT, A>::load_aligned(&mask_buffer[0]);
            }
        }

        template <typename A, typename T>
        XSIMD_INLINE batch<T, A>
        compress(batch<T, A> const& x, batch_bool<T, A> const& mask,
                 kernel::requires_arch<generic>) noexcept
        {
            using IT = as_unsigned_integer_t<T>;
            constexpr std::size_t size = batch_bool<T, A>::size;
            auto bitmask = mask.mask();
            auto z = select(mask, x, batch<T, A>((T)0));
            auto compress_mask = detail::create_compress_swizzle_mask<IT, A>(bitmask, ::xsimd::detail::make_index_sequence<size>());
            return swizzle(z, compress_mask);
        }

        // expand
        namespace detail
        {
            template <class IT, class A, class I, size_t... Is>
            XSIMD_INLINE batch<IT, A> create_expand_swizzle_mask(I bitmask, ::xsimd::detail::index_sequence<Is...>)
            {
                batch<IT, A> swizzle_mask(IT(0));
                IT j = 0;
                (void)std::initializer_list<bool> { ((swizzle_mask = insert(swizzle_mask, j, index<Is>())), (j += ((bitmask >> Is) & 1u)), true)... };
                return swizzle_mask;
            }
        }

        template <typename A, typename T>
        XSIMD_INLINE batch<T, A>
        expand(batch<T, A> const& x, batch_bool<T, A> const& mask,
               kernel::requires_arch<generic>) noexcept
        {
            constexpr std::size_t size = batch_bool<T, A>::size;
            auto bitmask = mask.mask();
            auto swizzle_mask = detail::create_expand_swizzle_mask<as_unsigned_integer_t<T>, A>(bitmask, ::xsimd::detail::make_index_sequence<size>());
            auto z = swizzle(x, swizzle_mask);
            return select(mask, z, batch<T, A>(T(0)));
        }

        // extract_pair
        template <class A, class T>
        XSIMD_INLINE batch<T, A> extract_pair(batch<T, A> const& self, batch<T, A> const& other, std::size_t i, requires_arch<generic>) noexcept
        {
            constexpr std::size_t size = batch<T, A>::size;
            assert(i < size && "index in bounds");

            alignas(A::alignment()) T self_buffer[size];
            self.store_aligned(self_buffer);

            alignas(A::alignment()) T other_buffer[size];
            other.store_aligned(other_buffer);

            alignas(A::alignment()) T concat_buffer[size];

            for (std::size_t j = 0; j < (size - i); ++j)
            {
                concat_buffer[j] = other_buffer[i + j];
                if (j < i)
                {
                    concat_buffer[size - 1 - j] = self_buffer[i - 1 - j];
                }
            }
            return batch<T, A>::load_aligned(concat_buffer);
        }

        // gather
        namespace detail
        {
            // Not using XSIMD_INLINE here as it makes msvc hand got ever on avx512
            template <size_t N, typename T, typename A, typename U, typename V, typename std::enable_if<N == 0, int>::type = 0>
            inline batch<T, A> gather(U const* src, batch<V, A> const& index,
                                      ::xsimd::index<N> I) noexcept
            {
                return insert(batch<T, A> {}, static_cast<T>(src[index.get(I)]), I);
            }

            template <size_t N, typename T, typename A, typename U, typename V, typename std::enable_if<N != 0, int>::type = 0>
            inline batch<T, A>
            gather(U const* src, batch<V, A> const& index, ::xsimd::index<N> I) noexcept
            {
                static_assert(N <= batch<V, A>::size, "Incorrect value in recursion!");

                const auto test = gather<N - 1, T, A>(src, index, {});
                return insert(test, static_cast<T>(src[index.get(I)]), I);
            }
        } // namespace detail

        template <typename T, typename A, typename V>
        XSIMD_INLINE batch<T, A>
        gather(batch<T, A> const&, T const* src, batch<V, A> const& index,
               kernel::requires_arch<generic>) noexcept
        {
            static_assert(batch<T, A>::size == batch<V, A>::size,
                          "Index and destination sizes must match");

            return detail::gather<batch<V, A>::size - 1, T, A>(src, index, {});
        }

        // Gather with runtime indexes and mismatched strides.
        template <typename T, typename A, typename U, typename V>
        XSIMD_INLINE detail::sizes_mismatch_t<T, U, batch<T, A>>
        gather(batch<T, A> const&, U const* src, batch<V, A> const& index,
               kernel::requires_arch<generic>) noexcept
        {
            static_assert(batch<T, A>::size == batch<V, A>::size,
                          "Index and destination sizes must match");

            return detail::gather<batch<V, A>::size - 1, T, A>(src, index, {});
        }

        // Gather with runtime indexes and matching strides.
        template <typename T, typename A, typename U, typename V>
        XSIMD_INLINE detail::stride_match_t<T, U, batch<T, A>>
        gather(batch<T, A> const&, U const* src, batch<V, A> const& index,
               kernel::requires_arch<generic>) noexcept
        {
            static_assert(batch<T, A>::size == batch<V, A>::size,
                          "Index and destination sizes must match");

            return batch_cast<T>(kernel::gather(batch<U, A> {}, src, index, A {}));
        }

        // insert
        template <class A, class T, size_t I>
        XSIMD_INLINE batch<T, A> insert(batch<T, A> const& self, T val, index<I>, requires_arch<generic>) noexcept
        {
            struct index_mask
            {
                static constexpr bool get(size_t index, size_t /* size*/)
                {
                    return index != I;
                }
            };
            batch<T, A> tmp(val);
            return select(make_batch_bool_constant<T, A, index_mask>(), self, tmp);
        }

        // get
        template <class A, size_t I, class T>
        XSIMD_INLINE T get(batch<T, A> const& self, ::xsimd::index<I>, requires_arch<generic>) noexcept
        {
            alignas(A::alignment()) T buffer[batch<T, A>::size];
            self.store_aligned(&buffer[0]);
            return buffer[I];
        }

        template <class A, size_t I, class T>
        XSIMD_INLINE T get(batch_bool<T, A> const& self, ::xsimd::index<I>, requires_arch<generic>) noexcept
        {
            alignas(A::alignment()) T buffer[batch_bool<T, A>::size];
            self.store_aligned(&buffer[0]);
            return buffer[I];
        }

        template <class A, size_t I, class T>
        XSIMD_INLINE auto get(batch<std::complex<T>, A> const& self, ::xsimd::index<I>, requires_arch<generic>) noexcept -> typename batch<std::complex<T>, A>::value_type
        {
            alignas(A::alignment()) T buffer[batch<std::complex<T>, A>::size];
            self.store_aligned(&buffer[0]);
            return buffer[I];
        }

        template <class A, class T>
        XSIMD_INLINE T get(batch<T, A> const& self, std::size_t i, requires_arch<generic>) noexcept
        {
            alignas(A::alignment()) T buffer[batch<T, A>::size];
            self.store_aligned(&buffer[0]);
            return buffer[i];
        }

        template <class A, class T>
        XSIMD_INLINE T get(batch_bool<T, A> const& self, std::size_t i, requires_arch<generic>) noexcept
        {
            alignas(A::alignment()) bool buffer[batch_bool<T, A>::size];
            self.store_aligned(&buffer[0]);
            return buffer[i];
        }

        template <class A, class T>
        XSIMD_INLINE auto get(batch<std::complex<T>, A> const& self, std::size_t i, requires_arch<generic>) noexcept -> typename batch<std::complex<T>, A>::value_type
        {
            using T2 = typename batch<std::complex<T>, A>::value_type;
            alignas(A::alignment()) T2 buffer[batch<std::complex<T>, A>::size];
            self.store_aligned(&buffer[0]);
            return buffer[i];
        }

        // load_aligned
        namespace detail
        {
            template <class A, class T_in, class T_out>
            XSIMD_INLINE batch<T_out, A> load_aligned(T_in const* mem, convert<T_out>, requires_arch<generic>, with_fast_conversion) noexcept
            {
                using batch_type_in = batch<T_in, A>;
                using batch_type_out = batch<T_out, A>;
                return fast_cast(batch_type_in::load_aligned(mem), batch_type_out(), A {});
            }
            template <class A, class T_in, class T_out>
            XSIMD_INLINE batch<T_out, A> load_aligned(T_in const* mem, convert<T_out>, requires_arch<generic>, with_slow_conversion) noexcept
            {
                static_assert(!std::is_same<T_in, T_out>::value, "there should be a direct load for this type combination");
                using batch_type_out = batch<T_out, A>;
                alignas(A::alignment()) T_out buffer[batch_type_out::size];
                std::copy(mem, mem + batch_type_out::size, std::begin(buffer));
                return batch_type_out::load_aligned(buffer);
            }
        }
        template <class A, class T_in, class T_out>
        XSIMD_INLINE batch<T_out, A> load_aligned(T_in const* mem, convert<T_out> cvt, requires_arch<generic>) noexcept
        {
            return detail::load_aligned<A>(mem, cvt, A {}, detail::conversion_type<A, T_in, T_out> {});
        }

        // load_unaligned
        namespace detail
        {
            template <class A, class T_in, class T_out>
            XSIMD_INLINE batch<T_out, A> load_unaligned(T_in const* mem, convert<T_out>, requires_arch<generic>, with_fast_conversion) noexcept
            {
                using batch_type_in = batch<T_in, A>;
                using batch_type_out = batch<T_out, A>;
                return fast_cast(batch_type_in::load_unaligned(mem), batch_type_out(), A {});
            }

            template <class A, class T_in, class T_out>
            XSIMD_INLINE batch<T_out, A> load_unaligned(T_in const* mem, convert<T_out> cvt, requires_arch<generic>, with_slow_conversion) noexcept
            {
                static_assert(!std::is_same<T_in, T_out>::value, "there should be a direct load for this type combination");
                return load_aligned<A>(mem, cvt, generic {}, with_slow_conversion {});
            }
        }
        template <class A, class T_in, class T_out>
        XSIMD_INLINE batch<T_out, A> load_unaligned(T_in const* mem, convert<T_out> cvt, requires_arch<generic>) noexcept
        {
            return detail::load_unaligned<A>(mem, cvt, generic {}, detail::conversion_type<A, T_in, T_out> {});
        }

        // rotate_left
        template <size_t N, class A, class T>
        XSIMD_INLINE batch<T, A> rotate_left(batch<T, A> const& self, requires_arch<generic>) noexcept
        {
            struct rotate_generator
            {
                static constexpr size_t get(size_t index, size_t size)
                {
                    return (index - N) % size;
                }
            };

            return swizzle(self, make_batch_constant<as_unsigned_integer_t<T>, A, rotate_generator>(), A {});
        }

        template <size_t N, class A, class T>
        XSIMD_INLINE batch<std::complex<T>, A> rotate_left(batch<std::complex<T>, A> const& self, requires_arch<generic>) noexcept
        {
            return { rotate_left<N>(self.real()), rotate_left<N>(self.imag()) };
        }

        // rotate_right
        template <size_t N, class A, class T>
        XSIMD_INLINE batch<T, A> rotate_right(batch<T, A> const& self, requires_arch<generic>) noexcept
        {
            struct rotate_generator
            {
                static constexpr size_t get(size_t index, size_t size)
                {
                    return (index + N) % size;
                }
            };

            return swizzle(self, make_batch_constant<as_unsigned_integer_t<T>, A, rotate_generator>(), A {});
        }

        template <size_t N, class A, class T>
        XSIMD_INLINE batch<std::complex<T>, A> rotate_right(batch<std::complex<T>, A> const& self, requires_arch<generic>) noexcept
        {
            return { rotate_right<N>(self.real()), rotate_right<N>(self.imag()) };
        }

        // Scatter with runtime indexes.
        namespace detail
        {
            template <size_t N, typename T, typename A, typename U, typename V, typename std::enable_if<N == 0, int>::type = 0>
            XSIMD_INLINE void scatter(batch<T, A> const& src, U* dst,
                                      batch<V, A> const& index,
                                      ::xsimd::index<N> I) noexcept
            {
                dst[index.get(I)] = static_cast<U>(src.get(I));
            }

            template <size_t N, typename T, typename A, typename U, typename V, typename std::enable_if<N != 0, int>::type = 0>
            XSIMD_INLINE void
            scatter(batch<T, A> const& src, U* dst, batch<V, A> const& index,
                    ::xsimd::index<N> I) noexcept
            {
                static_assert(N <= batch<V, A>::size, "Incorrect value in recursion!");

                kernel::detail::scatter<N - 1, T, A, U, V>(
                    src, dst, index, {});
                dst[index.get(I)] = static_cast<U>(src.get(I));
            }
        } // namespace detail

        template <typename A, typename T, typename V>
        XSIMD_INLINE void
        scatter(batch<T, A> const& src, T* dst,
                batch<V, A> const& index,
                kernel::requires_arch<generic>) noexcept
        {
            static_assert(batch<T, A>::size == batch<V, A>::size,
                          "Source and index sizes must match");
            kernel::detail::scatter<batch<V, A>::size - 1, T, A, T, V>(
                src, dst, index, {});
        }

        template <typename A, typename T, typename U, typename V>
        XSIMD_INLINE detail::sizes_mismatch_t<T, U, void>
        scatter(batch<T, A> const& src, U* dst,
                batch<V, A> const& index,
                kernel::requires_arch<generic>) noexcept
        {
            static_assert(batch<T, A>::size == batch<V, A>::size,
                          "Source and index sizes must match");
            kernel::detail::scatter<batch<V, A>::size - 1, T, A, U, V>(
                src, dst, index, {});
        }

        template <typename A, typename T, typename U, typename V>
        XSIMD_INLINE detail::stride_match_t<T, U, void>
        scatter(batch<T, A> const& src, U* dst,
                batch<V, A> const& index,
                kernel::requires_arch<generic>) noexcept
        {
            static_assert(batch<T, A>::size == batch<V, A>::size,
                          "Source and index sizes must match");
            const auto tmp = batch_cast<U>(src);
            kernel::scatter<A>(tmp, dst, index, A {});
        }

        // shuffle
        namespace detail
        {
            constexpr bool is_swizzle_fst(size_t)
            {
                return true;
            }
            template <typename ITy, typename... ITys>
            constexpr bool is_swizzle_fst(size_t bsize, ITy index, ITys... indices)
            {
                return index < bsize && is_swizzle_fst(bsize, indices...);
            }
            constexpr bool is_swizzle_snd(size_t)
            {
                return true;
            }
            template <typename ITy, typename... ITys>
            constexpr bool is_swizzle_snd(size_t bsize, ITy index, ITys... indices)
            {
                return index >= bsize && is_swizzle_snd(bsize, indices...);
            }

            constexpr bool is_zip_lo(size_t)
            {
                return true;
            }

            template <typename ITy>
            constexpr bool is_zip_lo(size_t, ITy)
            {
                return false;
            }

            template <typename ITy0, typename ITy1, typename... ITys>
            constexpr bool is_zip_lo(size_t bsize, ITy0 index0, ITy1 index1, ITys... indices)
            {
                return index0 == (bsize - (sizeof...(indices) + 2)) && index1 == (2 * bsize - (sizeof...(indices) + 2)) && is_zip_lo(bsize, indices...);
            }

            constexpr bool is_zip_hi(size_t)
            {
                return true;
            }

            template <typename ITy>
            constexpr bool is_zip_hi(size_t, ITy)
            {
                return false;
            }

            template <typename ITy0, typename ITy1, typename... ITys>
            constexpr bool is_zip_hi(size_t bsize, ITy0 index0, ITy1 index1, ITys... indices)
            {
                return index0 == (bsize / 2 + bsize - (sizeof...(indices) + 2)) && index1 == (bsize / 2 + 2 * bsize - (sizeof...(indices) + 2)) && is_zip_hi(bsize, indices...);
            }

            constexpr bool is_select(size_t)
            {
                return true;
            }

            template <typename ITy, typename... ITys>
            constexpr bool is_select(size_t bsize, ITy index, ITys... indices)
            {
                return (index < bsize ? index : index - bsize) == (bsize - sizeof...(ITys)) && is_select(bsize, indices...);
            }

        }

        template <class A, typename T, typename ITy, ITy... Indices>
        XSIMD_INLINE batch<T, A> shuffle(batch<T, A> const& x, batch<T, A> const& y, batch_constant<ITy, A, Indices...>, requires_arch<generic>) noexcept
        {
            constexpr size_t bsize = sizeof...(Indices);

            // Detect common patterns
            XSIMD_IF_CONSTEXPR(detail::is_swizzle_fst(bsize, Indices...))
            {
                return swizzle(x, batch_constant<ITy, A, ((Indices >= bsize) ? 0 /* never happens */ : Indices)...>());
            }

            XSIMD_IF_CONSTEXPR(detail::is_swizzle_snd(bsize, Indices...))
            {
                return swizzle(y, batch_constant<ITy, A, ((Indices >= bsize) ? (Indices - bsize) : 0 /* never happens */)...>());
            }

            XSIMD_IF_CONSTEXPR(detail::is_zip_lo(bsize, Indices...))
            {
                return zip_lo(x, y);
            }

            XSIMD_IF_CONSTEXPR(detail::is_zip_hi(bsize, Indices...))
            {
                return zip_hi(x, y);
            }

            XSIMD_IF_CONSTEXPR(detail::is_select(bsize, Indices...))
            {
                return select(batch_bool_constant<T, A, (Indices < bsize)...>(), x, y);
            }

#if defined(__has_builtin)
#if __has_builtin(__builtin_shuffle_vector)
#define builtin_shuffle __builtin_shuffle_vector
#endif
#endif

#if defined(builtin_shuffle)
            return builtin_shuffle(x.data, y.data, Indices...);

// FIXME: my experiments show that GCC only correctly optimizes this builtin
// starting at GCC 13, where it already has __builtin_shuffle_vector
//
// #elif __has_builtin(__builtin_shuffle) || GCC >= 6
//            typedef ITy integer_vector_type __attribute__((vector_size(sizeof(batch<ITy, A>))));
//            return __builtin_shuffle(x.data, y.data, integer_vector_type{Indices...});
#else
            // Use a generic_pattern. It is suboptimal but clang optimizes this
            // pretty well.
            batch<T, A> x_lane = swizzle(x, batch_constant<ITy, A, ((Indices >= bsize) ? (Indices - bsize) : Indices)...>());
            batch<T, A> y_lane = swizzle(y, batch_constant<ITy, A, ((Indices >= bsize) ? (Indices - bsize) : Indices)...>());
            batch_bool_constant<T, A, (Indices < bsize)...> select_x_lane;
            return select(select_x_lane, x_lane, y_lane);
#endif
        }

        // store
        template <class T, class A>
        XSIMD_INLINE void store(batch_bool<T, A> const& self, bool* mem, requires_arch<generic>) noexcept
        {
            using batch_type = batch<T, A>;
            constexpr auto size = batch_bool<T, A>::size;
            alignas(A::alignment()) T buffer[size];
            kernel::store_aligned<A>(&buffer[0], batch_type(self), A {});
            for (std::size_t i = 0; i < size; ++i)
                mem[i] = bool(buffer[i]);
        }

        // store_aligned
        template <class A, class T_in, class T_out>
        XSIMD_INLINE void store_aligned(T_out* mem, batch<T_in, A> const& self, requires_arch<generic>) noexcept
        {
            static_assert(!std::is_same<T_in, T_out>::value, "there should be a direct store for this type combination");
            alignas(A::alignment()) T_in buffer[batch<T_in, A>::size];
            store_aligned(&buffer[0], self);
            std::copy(std::begin(buffer), std::end(buffer), mem);
        }

        // store_unaligned
        template <class A, class T_in, class T_out>
        XSIMD_INLINE void store_unaligned(T_out* mem, batch<T_in, A> const& self, requires_arch<generic>) noexcept
        {
            static_assert(!std::is_same<T_in, T_out>::value, "there should be a direct store for this type combination");
            return store_aligned<A>(mem, self, generic {});
        }

        // swizzle
        template <class A, class T, class ITy, ITy... Vs>
        XSIMD_INLINE batch<std::complex<T>, A> swizzle(batch<std::complex<T>, A> const& self, batch_constant<ITy, A, Vs...> mask, requires_arch<generic>) noexcept
        {
            return { swizzle(self.real(), mask), swizzle(self.imag(), mask) };
        }

        template <class A, class T, class ITy>
        XSIMD_INLINE batch<T, A> swizzle(batch<T, A> const& self, batch<ITy, A> mask, requires_arch<generic>) noexcept
        {
            constexpr size_t size = batch<T, A>::size;
            alignas(A::alignment()) T self_buffer[size];
            store_aligned(&self_buffer[0], self);

            alignas(A::alignment()) ITy mask_buffer[size];
            store_aligned(&mask_buffer[0], mask);

            alignas(A::alignment()) T out_buffer[size];
            for (size_t i = 0; i < size; ++i)
                out_buffer[i] = self_buffer[mask_buffer[i]];
            return batch<T, A>::load_aligned(out_buffer);
        }

        template <class A, class T, class ITy>
        XSIMD_INLINE batch<std::complex<T>, A> swizzle(batch<std::complex<T>, A> const& self, batch<ITy, A> mask, requires_arch<generic>) noexcept
        {
            return { swizzle(self.real(), mask), swizzle(self.imag(), mask) };
        }

        // load_complex_aligned
        namespace detail
        {
            template <class A, class T>
            XSIMD_INLINE batch<std::complex<T>, A> load_complex(batch<T, A> const& /*hi*/, batch<T, A> const& /*lo*/, requires_arch<generic>) noexcept
            {
                static_assert(std::is_same<T, void>::value, "load_complex not implemented for the required architecture");
            }

            template <class A, class T>
            XSIMD_INLINE batch<T, A> complex_high(batch<std::complex<T>, A> const& /*src*/, requires_arch<generic>) noexcept
            {
                static_assert(std::is_same<T, void>::value, "complex_high not implemented for the required architecture");
            }

            template <class A, class T>
            XSIMD_INLINE batch<T, A> complex_low(batch<std::complex<T>, A> const& /*src*/, requires_arch<generic>) noexcept
            {
                static_assert(std::is_same<T, void>::value, "complex_low not implemented for the required architecture");
            }
        }

        template <class A, class T_out, class T_in>
        XSIMD_INLINE batch<std::complex<T_out>, A> load_complex_aligned(std::complex<T_in> const* mem, convert<std::complex<T_out>>, requires_arch<generic>) noexcept
        {
            using real_batch = batch<T_out, A>;
            T_in const* buffer = reinterpret_cast<T_in const*>(mem);
            real_batch hi = real_batch::load_aligned(buffer),
                       lo = real_batch::load_aligned(buffer + real_batch::size);
            return detail::load_complex(hi, lo, A {});
        }

        // load_complex_unaligned
        template <class A, class T_out, class T_in>
        XSIMD_INLINE batch<std::complex<T_out>, A> load_complex_unaligned(std::complex<T_in> const* mem, convert<std::complex<T_out>>, requires_arch<generic>) noexcept
        {
            using real_batch = batch<T_out, A>;
            T_in const* buffer = reinterpret_cast<T_in const*>(mem);
            real_batch hi = real_batch::load_unaligned(buffer),
                       lo = real_batch::load_unaligned(buffer + real_batch::size);
            return detail::load_complex(hi, lo, A {});
        }

        // store_complex_aligned
        template <class A, class T_out, class T_in>
        XSIMD_INLINE void store_complex_aligned(std::complex<T_out>* dst, batch<std::complex<T_in>, A> const& src, requires_arch<generic>) noexcept
        {
            using real_batch = batch<T_in, A>;
            real_batch hi = detail::complex_high(src, A {});
            real_batch lo = detail::complex_low(src, A {});
            T_out* buffer = reinterpret_cast<T_out*>(dst);
            lo.store_aligned(buffer);
            hi.store_aligned(buffer + real_batch::size);
        }

        // store_compelx_unaligned
        template <class A, class T_out, class T_in>
        XSIMD_INLINE void store_complex_unaligned(std::complex<T_out>* dst, batch<std::complex<T_in>, A> const& src, requires_arch<generic>) noexcept
        {
            using real_batch = batch<T_in, A>;
            real_batch hi = detail::complex_high(src, A {});
            real_batch lo = detail::complex_low(src, A {});
            T_out* buffer = reinterpret_cast<T_out*>(dst);
            lo.store_unaligned(buffer);
            hi.store_unaligned(buffer + real_batch::size);
        }

    }

}

#endif
