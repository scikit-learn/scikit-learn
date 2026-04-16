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

#ifndef XSIMD_AVX_HPP
#define XSIMD_AVX_HPP

#include <complex>
#include <limits>
#include <type_traits>

#include "../types/xsimd_avx_register.hpp"

namespace xsimd
{

    namespace kernel
    {
        using namespace types;

        // fwd
        template <class A, class T, size_t I>
        XSIMD_INLINE batch<T, A> insert(batch<T, A> const& self, T val, index<I>, requires_arch<generic>) noexcept;

        namespace detail
        {
            XSIMD_INLINE void split_avx(__m256i val, __m128i& low, __m128i& high) noexcept
            {
                low = _mm256_castsi256_si128(val);
                high = _mm256_extractf128_si256(val, 1);
            }
            XSIMD_INLINE void split_avx(__m256 val, __m128& low, __m128& high) noexcept
            {
                low = _mm256_castps256_ps128(val);
                high = _mm256_extractf128_ps(val, 1);
            }
            XSIMD_INLINE void split_avx(__m256d val, __m128d& low, __m128d& high) noexcept
            {
                low = _mm256_castpd256_pd128(val);
                high = _mm256_extractf128_pd(val, 1);
            }
            XSIMD_INLINE __m256i merge_sse(__m128i low, __m128i high) noexcept
            {
                return _mm256_insertf128_si256(_mm256_castsi128_si256(low), high, 1);
            }
            XSIMD_INLINE __m256 merge_sse(__m128 low, __m128 high) noexcept
            {
                return _mm256_insertf128_ps(_mm256_castps128_ps256(low), high, 1);
            }
            XSIMD_INLINE __m256d merge_sse(__m128d low, __m128d high) noexcept
            {
                return _mm256_insertf128_pd(_mm256_castpd128_pd256(low), high, 1);
            }
            template <class F>
            XSIMD_INLINE __m256i fwd_to_sse(F f, __m256i self) noexcept
            {
                __m128i self_low, self_high;
                split_avx(self, self_low, self_high);
                __m128i res_low = f(self_low);
                __m128i res_high = f(self_high);
                return merge_sse(res_low, res_high);
            }
            template <class F>
            XSIMD_INLINE __m256i fwd_to_sse(F f, __m256i self, __m256i other) noexcept
            {
                __m128i self_low, self_high, other_low, other_high;
                split_avx(self, self_low, self_high);
                split_avx(other, other_low, other_high);
                __m128i res_low = f(self_low, other_low);
                __m128i res_high = f(self_high, other_high);
                return merge_sse(res_low, res_high);
            }
            template <class F>
            XSIMD_INLINE __m256i fwd_to_sse(F f, __m256i self, int32_t other) noexcept
            {
                __m128i self_low, self_high;
                split_avx(self, self_low, self_high);
                __m128i res_low = f(self_low, other);
                __m128i res_high = f(self_high, other);
                return merge_sse(res_low, res_high);
            }
        }

        // abs
        template <class A>
        XSIMD_INLINE batch<float, A> abs(batch<float, A> const& self, requires_arch<avx>) noexcept
        {
            __m256 sign_mask = _mm256_set1_ps(-0.f); // -0.f = 1 << 31
            return _mm256_andnot_ps(sign_mask, self);
        }
        template <class A>
        XSIMD_INLINE batch<double, A> abs(batch<double, A> const& self, requires_arch<avx>) noexcept
        {
            __m256d sign_mask = _mm256_set1_pd(-0.f); // -0.f = 1 << 31
            return _mm256_andnot_pd(sign_mask, self);
        }

        // add
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> add(batch<T, A> const& self, batch<T, A> const& other, requires_arch<avx>) noexcept
        {
            return detail::fwd_to_sse([](__m128i s, __m128i o) noexcept
                                      { return add(batch<T, sse4_2>(s), batch<T, sse4_2>(o)); },
                                      self, other);
        }
        template <class A>
        XSIMD_INLINE batch<float, A> add(batch<float, A> const& self, batch<float, A> const& other, requires_arch<avx>) noexcept
        {
            return _mm256_add_ps(self, other);
        }
        template <class A>
        XSIMD_INLINE batch<double, A> add(batch<double, A> const& self, batch<double, A> const& other, requires_arch<avx>) noexcept
        {
            return _mm256_add_pd(self, other);
        }

        // all
        template <class A>
        XSIMD_INLINE bool all(batch_bool<float, A> const& self, requires_arch<avx>) noexcept
        {
            return _mm256_testc_ps(self, batch_bool<float, A>(true)) != 0;
        }
        template <class A>
        XSIMD_INLINE bool all(batch_bool<double, A> const& self, requires_arch<avx>) noexcept
        {
            return _mm256_testc_pd(self, batch_bool<double, A>(true)) != 0;
        }
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE bool all(batch_bool<T, A> const& self, requires_arch<avx>) noexcept
        {
            return _mm256_testc_si256(self, batch_bool<T, A>(true)) != 0;
        }

        // any
        template <class A>
        XSIMD_INLINE bool any(batch_bool<float, A> const& self, requires_arch<avx>) noexcept
        {
            return !_mm256_testz_ps(self, self);
        }
        template <class A>
        XSIMD_INLINE bool any(batch_bool<double, A> const& self, requires_arch<avx>) noexcept
        {
            return !_mm256_testz_pd(self, self);
        }
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE bool any(batch_bool<T, A> const& self, requires_arch<avx>) noexcept
        {
            return !_mm256_testz_si256(self, self);
        }

        // batch_bool_cast
        template <class A, class T_out, class T_in>
        XSIMD_INLINE batch_bool<T_out, A> batch_bool_cast(batch_bool<T_in, A> const& self, batch_bool<T_out, A> const&, requires_arch<avx>) noexcept
        {
            return { bitwise_cast<T_out>(batch<T_in, A>(self.data)).data };
        }

        // bitwise_and
        template <class A>
        XSIMD_INLINE batch<float, A> bitwise_and(batch<float, A> const& self, batch<float, A> const& other, requires_arch<avx>) noexcept
        {
            return _mm256_and_ps(self, other);
        }
        template <class A>
        XSIMD_INLINE batch<double, A> bitwise_and(batch<double, A> const& self, batch<double, A> const& other, requires_arch<avx>) noexcept
        {
            return _mm256_and_pd(self, other);
        }

        template <class A>
        XSIMD_INLINE batch_bool<float, A> bitwise_and(batch_bool<float, A> const& self, batch_bool<float, A> const& other, requires_arch<avx>) noexcept
        {
            return _mm256_and_ps(self, other);
        }
        template <class A>
        XSIMD_INLINE batch_bool<double, A> bitwise_and(batch_bool<double, A> const& self, batch_bool<double, A> const& other, requires_arch<avx>) noexcept
        {
            return _mm256_and_pd(self, other);
        }

        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> bitwise_and(batch<T, A> const& self, batch<T, A> const& other, requires_arch<avx>) noexcept
        {
            return detail::fwd_to_sse([](__m128i s, __m128i o) noexcept
                                      { return bitwise_and(batch<T, sse4_2>(s), batch<T, sse4_2>(o)); },
                                      self, other);
        }
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch_bool<T, A> bitwise_and(batch_bool<T, A> const& self, batch_bool<T, A> const& other, requires_arch<avx>) noexcept
        {
            return detail::fwd_to_sse([](__m128i s, __m128i o) noexcept
                                      { return bitwise_and(batch<T, sse4_2>(s), batch<T, sse4_2>(o)); },
                                      self, other);
        }

        // bitwise_andnot
        template <class A>
        XSIMD_INLINE batch<float, A> bitwise_andnot(batch<float, A> const& self, batch<float, A> const& other, requires_arch<avx>) noexcept
        {
            return _mm256_andnot_ps(other, self);
        }
        template <class A>
        XSIMD_INLINE batch<double, A> bitwise_andnot(batch<double, A> const& self, batch<double, A> const& other, requires_arch<avx>) noexcept
        {
            return _mm256_andnot_pd(other, self);
        }

        template <class A>
        XSIMD_INLINE batch_bool<float, A> bitwise_andnot(batch_bool<float, A> const& self, batch_bool<float, A> const& other, requires_arch<avx>) noexcept
        {
            return _mm256_andnot_ps(other, self);
        }
        template <class A>
        XSIMD_INLINE batch_bool<double, A> bitwise_andnot(batch_bool<double, A> const& self, batch_bool<double, A> const& other, requires_arch<avx>) noexcept
        {
            return _mm256_andnot_pd(other, self);
        }

        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> bitwise_andnot(batch<T, A> const& self, batch<T, A> const& other, requires_arch<avx>) noexcept
        {
            return detail::fwd_to_sse([](__m128i s, __m128i o) noexcept
                                      { return bitwise_andnot(batch<T, sse4_2>(s), batch<T, sse4_2>(o)); },
                                      self, other);
        }
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch_bool<T, A> bitwise_andnot(batch_bool<T, A> const& self, batch_bool<T, A> const& other, requires_arch<avx>) noexcept
        {
            return detail::fwd_to_sse([](__m128i s, __m128i o) noexcept
                                      { return bitwise_andnot(batch<T, sse4_2>(s), batch<T, sse4_2>(o)); },
                                      self, other);
        }

        // bitwise_lshift
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> bitwise_lshift(batch<T, A> const& self, int32_t other, requires_arch<avx>) noexcept
        {
            return detail::fwd_to_sse([](__m128i s, int32_t o) noexcept
                                      { return bitwise_lshift(batch<T, sse4_2>(s), o, sse4_2 {}); },
                                      self, other);
        }

        // bitwise_not
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> bitwise_not(batch<T, A> const& self, requires_arch<avx>) noexcept
        {
            return detail::fwd_to_sse([](__m128i s) noexcept
                                      { return bitwise_not(batch<T, sse4_2>(s), sse4_2 {}); },
                                      self);
        }
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch_bool<T, A> bitwise_not(batch_bool<T, A> const& self, requires_arch<avx>) noexcept
        {
            return detail::fwd_to_sse([](__m128i s) noexcept
                                      { return bitwise_not(batch_bool<T, sse4_2>(s), sse4_2 {}); },
                                      self);
        }

        // bitwise_or
        template <class A>
        XSIMD_INLINE batch<float, A> bitwise_or(batch<float, A> const& self, batch<float, A> const& other, requires_arch<avx>) noexcept
        {
            return _mm256_or_ps(self, other);
        }
        template <class A>
        XSIMD_INLINE batch<double, A> bitwise_or(batch<double, A> const& self, batch<double, A> const& other, requires_arch<avx>) noexcept
        {
            return _mm256_or_pd(self, other);
        }
        template <class A>
        XSIMD_INLINE batch_bool<float, A> bitwise_or(batch_bool<float, A> const& self, batch_bool<float, A> const& other, requires_arch<avx>) noexcept
        {
            return _mm256_or_ps(self, other);
        }
        template <class A>
        XSIMD_INLINE batch_bool<double, A> bitwise_or(batch_bool<double, A> const& self, batch_bool<double, A> const& other, requires_arch<avx>) noexcept
        {
            return _mm256_or_pd(self, other);
        }
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> bitwise_or(batch<T, A> const& self, batch<T, A> const& other, requires_arch<avx>) noexcept
        {
            return detail::fwd_to_sse([](__m128i s, __m128i o) noexcept
                                      { return bitwise_or(batch<T, sse4_2>(s), batch<T, sse4_2>(o)); },
                                      self, other);
        }
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch_bool<T, A> bitwise_or(batch_bool<T, A> const& self, batch_bool<T, A> const& other, requires_arch<avx>) noexcept
        {
            return detail::fwd_to_sse([](__m128i s, __m128i o) noexcept
                                      { return bitwise_or(batch_bool<T, sse4_2>(s), batch_bool<T, sse4_2>(o)); },
                                      self, other);
        }

        // bitwise_rshift
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> bitwise_rshift(batch<T, A> const& self, int32_t other, requires_arch<avx>) noexcept
        {
            return detail::fwd_to_sse([](__m128i s, int32_t o) noexcept
                                      { return bitwise_rshift(batch<T, sse4_2>(s), o, sse4_2 {}); },
                                      self, other);
        }

        // bitwise_xor
        template <class A>
        XSIMD_INLINE batch<float, A> bitwise_xor(batch<float, A> const& self, batch<float, A> const& other, requires_arch<avx>) noexcept
        {
            return _mm256_xor_ps(self, other);
        }
        template <class A>
        XSIMD_INLINE batch<double, A> bitwise_xor(batch<double, A> const& self, batch<double, A> const& other, requires_arch<avx>) noexcept
        {
            return _mm256_xor_pd(self, other);
        }
        template <class A>
        XSIMD_INLINE batch_bool<float, A> bitwise_xor(batch_bool<float, A> const& self, batch_bool<float, A> const& other, requires_arch<avx>) noexcept
        {
            return _mm256_xor_ps(self, other);
        }
        template <class A>
        XSIMD_INLINE batch_bool<double, A> bitwise_xor(batch_bool<double, A> const& self, batch_bool<double, A> const& other, requires_arch<avx>) noexcept
        {
            return _mm256_xor_pd(self, other);
        }
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> bitwise_xor(batch<T, A> const& self, batch<T, A> const& other, requires_arch<avx>) noexcept
        {
            return detail::fwd_to_sse([](__m128i s, __m128i o) noexcept
                                      { return bitwise_xor(batch<T, sse4_2>(s), batch<T, sse4_2>(o), sse4_2 {}); },
                                      self, other);
        }
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> bitwise_xor(batch_bool<T, A> const& self, batch_bool<T, A> const& other, requires_arch<avx>) noexcept
        {
            return detail::fwd_to_sse([](__m128i s, __m128i o) noexcept
                                      { return bitwise_xor(batch_bool<T, sse4_2>(s), batch_bool<T, sse4_2>(o), sse4_2 {}); },
                                      self, other);
        }

        // bitwise_cast
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<float, A> bitwise_cast(batch<T, A> const& self, batch<float, A> const&, requires_arch<avx>) noexcept
        {
            return _mm256_castsi256_ps(self);
        }
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<double, A> bitwise_cast(batch<T, A> const& self, batch<double, A> const&, requires_arch<avx>) noexcept
        {
            return _mm256_castsi256_pd(self);
        }
        template <class A, class T, class Tp, class = typename std::enable_if<std::is_integral<std::common_type_t<T, Tp>>::value, void>::type>
        XSIMD_INLINE batch<Tp, A> bitwise_cast(batch<T, A> const& self, batch<Tp, A> const&, requires_arch<avx>) noexcept
        {
            return batch<Tp, A>(self.data);
        }
        template <class A>
        XSIMD_INLINE batch<double, A> bitwise_cast(batch<float, A> const& self, batch<double, A> const&, requires_arch<avx>) noexcept
        {
            return _mm256_castps_pd(self);
        }
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> bitwise_cast(batch<float, A> const& self, batch<T, A> const&, requires_arch<avx>) noexcept
        {
            return _mm256_castps_si256(self);
        }
        template <class A>
        XSIMD_INLINE batch<float, A> bitwise_cast(batch<double, A> const& self, batch<float, A> const&, requires_arch<avx>) noexcept
        {
            return _mm256_castpd_ps(self);
        }
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> bitwise_cast(batch<double, A> const& self, batch<T, A> const&, requires_arch<avx>) noexcept
        {
            return _mm256_castpd_si256(self);
        }

        // bitwise_not
        template <class A>
        XSIMD_INLINE batch<float, A> bitwise_not(batch<float, A> const& self, requires_arch<avx>) noexcept
        {
            return _mm256_xor_ps(self, _mm256_castsi256_ps(_mm256_set1_epi32(-1)));
        }
        template <class A>
        XSIMD_INLINE batch<double, A> bitwise_not(batch<double, A> const& self, requires_arch<avx>) noexcept
        {
            return _mm256_xor_pd(self, _mm256_castsi256_pd(_mm256_set1_epi32(-1)));
        }
        template <class A>
        XSIMD_INLINE batch_bool<float, A> bitwise_not(batch_bool<float, A> const& self, requires_arch<avx>) noexcept
        {
            return _mm256_xor_ps(self, _mm256_castsi256_ps(_mm256_set1_epi32(-1)));
        }
        template <class A>
        XSIMD_INLINE batch_bool<double, A> bitwise_not(batch_bool<double, A> const& self, requires_arch<avx>) noexcept
        {
            return _mm256_xor_pd(self, _mm256_castsi256_pd(_mm256_set1_epi32(-1)));
        }

        // broadcast
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> broadcast(T val, requires_arch<avx>) noexcept
        {
            XSIMD_IF_CONSTEXPR(sizeof(T) == 1)
            {
                return _mm256_set1_epi8(val);
            }
            else XSIMD_IF_CONSTEXPR(sizeof(T) == 2)
            {
                return _mm256_set1_epi16(val);
            }
            else XSIMD_IF_CONSTEXPR(sizeof(T) == 4)
            {
                return _mm256_set1_epi32(val);
            }
            else XSIMD_IF_CONSTEXPR(sizeof(T) == 8)
            {
                return _mm256_set1_epi64x(val);
            }
            else
            {
                assert(false && "unsupported");
                return {};
            }
        }
        template <class A>
        XSIMD_INLINE batch<float, A> broadcast(float val, requires_arch<avx>) noexcept
        {
            return _mm256_set1_ps(val);
        }
        template <class A>
        XSIMD_INLINE batch<double, A> broadcast(double val, requires_arch<avx>) noexcept
        {
            return _mm256_set1_pd(val);
        }

        // ceil
        template <class A>
        XSIMD_INLINE batch<float, A> ceil(batch<float, A> const& self, requires_arch<avx>) noexcept
        {
            return _mm256_ceil_ps(self);
        }
        template <class A>
        XSIMD_INLINE batch<double, A> ceil(batch<double, A> const& self, requires_arch<avx>) noexcept
        {
            return _mm256_ceil_pd(self);
        }

        namespace detail
        {
            // On clang, _mm256_extractf128_ps is built upon build_shufflevector
            // which require index parameter to be a constant
            template <int index, class B>
            XSIMD_INLINE B get_half_complex_f(const B& real, const B& imag) noexcept
            {
                __m128 tmp0 = _mm256_extractf128_ps(real, index);
                __m128 tmp1 = _mm256_extractf128_ps(imag, index);
                __m128 tmp2 = _mm_unpackhi_ps(tmp0, tmp1);
                tmp0 = _mm_unpacklo_ps(tmp0, tmp1);
                __m256 res = real;
                res = _mm256_insertf128_ps(res, tmp0, 0);
                res = _mm256_insertf128_ps(res, tmp2, 1);
                return res;
            }
            template <int index, class B>
            XSIMD_INLINE B get_half_complex_d(const B& real, const B& imag) noexcept
            {
                __m128d tmp0 = _mm256_extractf128_pd(real, index);
                __m128d tmp1 = _mm256_extractf128_pd(imag, index);
                __m128d tmp2 = _mm_unpackhi_pd(tmp0, tmp1);
                tmp0 = _mm_unpacklo_pd(tmp0, tmp1);
                __m256d res = real;
                res = _mm256_insertf128_pd(res, tmp0, 0);
                res = _mm256_insertf128_pd(res, tmp2, 1);
                return res;
            }

            // complex_low
            template <class A>
            XSIMD_INLINE batch<float, A> complex_low(batch<std::complex<float>, A> const& self, requires_arch<avx>) noexcept
            {
                return get_half_complex_f<0>(self.real(), self.imag());
            }
            template <class A>
            XSIMD_INLINE batch<double, A> complex_low(batch<std::complex<double>, A> const& self, requires_arch<avx>) noexcept
            {
                return get_half_complex_d<0>(self.real(), self.imag());
            }

            // complex_high
            template <class A>
            XSIMD_INLINE batch<float, A> complex_high(batch<std::complex<float>, A> const& self, requires_arch<avx>) noexcept
            {
                return get_half_complex_f<1>(self.real(), self.imag());
            }
            template <class A>
            XSIMD_INLINE batch<double, A> complex_high(batch<std::complex<double>, A> const& self, requires_arch<avx>) noexcept
            {
                return get_half_complex_d<1>(self.real(), self.imag());
            }
        }

        // fast_cast
        namespace detail
        {
            template <class A>
            XSIMD_INLINE batch<float, A> fast_cast(batch<int32_t, A> const& self, batch<float, A> const&, requires_arch<avx>) noexcept
            {
                return _mm256_cvtepi32_ps(self);
            }

            template <class A>
            XSIMD_INLINE batch<int32_t, A> fast_cast(batch<float, A> const& self, batch<int32_t, A> const&, requires_arch<avx>) noexcept
            {
                return _mm256_cvttps_epi32(self);
            }
        }

        // decr_if
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> decr_if(batch<T, A> const& self, batch_bool<T, A> const& mask, requires_arch<avx>) noexcept
        {
            return self + batch<T, A>(mask.data);
        }

        // div
        template <class A>
        XSIMD_INLINE batch<float, A> div(batch<float, A> const& self, batch<float, A> const& other, requires_arch<avx>) noexcept
        {
            return _mm256_div_ps(self, other);
        }
        template <class A>
        XSIMD_INLINE batch<double, A> div(batch<double, A> const& self, batch<double, A> const& other, requires_arch<avx>) noexcept
        {
            return _mm256_div_pd(self, other);
        }

        // eq
        template <class A>
        XSIMD_INLINE batch_bool<float, A> eq(batch<float, A> const& self, batch<float, A> const& other, requires_arch<avx>) noexcept
        {
            return _mm256_cmp_ps(self, other, _CMP_EQ_OQ);
        }
        template <class A>
        XSIMD_INLINE batch_bool<double, A> eq(batch<double, A> const& self, batch<double, A> const& other, requires_arch<avx>) noexcept
        {
            return _mm256_cmp_pd(self, other, _CMP_EQ_OQ);
        }
        template <class A>
        XSIMD_INLINE batch_bool<float, A> eq(batch_bool<float, A> const& self, batch_bool<float, A> const& other, requires_arch<avx>) noexcept
        {
            return ~(self != other);
        }
        template <class A>
        XSIMD_INLINE batch_bool<double, A> eq(batch_bool<double, A> const& self, batch_bool<double, A> const& other, requires_arch<avx>) noexcept
        {
            return ~(self != other);
        }
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch_bool<T, A> eq(batch<T, A> const& self, batch<T, A> const& other, requires_arch<avx>) noexcept
        {
            return detail::fwd_to_sse([](__m128i s, __m128i o) noexcept
                                      { return eq(batch<T, sse4_2>(s), batch<T, sse4_2>(o), sse4_2 {}); },
                                      self, other);
        }

        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch_bool<T, A> eq(batch_bool<T, A> const& self, batch_bool<T, A> const& other, requires_arch<avx>) noexcept
        {
            return ~(self != other);
        }

        // floor
        template <class A>
        XSIMD_INLINE batch<float, A> floor(batch<float, A> const& self, requires_arch<avx>) noexcept
        {
            return _mm256_floor_ps(self);
        }
        template <class A>
        XSIMD_INLINE batch<double, A> floor(batch<double, A> const& self, requires_arch<avx>) noexcept
        {
            return _mm256_floor_pd(self);
        }

        // from_mask
        template <class A>
        XSIMD_INLINE batch_bool<float, A> from_mask(batch_bool<float, A> const&, uint64_t mask, requires_arch<avx>) noexcept
        {
            alignas(A::alignment()) static const uint64_t lut32[] = {
                0x0000000000000000ul,
                0x00000000FFFFFFFFul,
                0xFFFFFFFF00000000ul,
                0xFFFFFFFFFFFFFFFFul,
            };
            assert(!(mask & ~0xFFul) && "inbound mask");
            return _mm256_castsi256_ps(_mm256_setr_epi64x(lut32[mask & 0x3], lut32[(mask >> 2) & 0x3], lut32[(mask >> 4) & 0x3], lut32[mask >> 6]));
        }
        template <class A>
        XSIMD_INLINE batch_bool<double, A> from_mask(batch_bool<double, A> const&, uint64_t mask, requires_arch<avx>) noexcept
        {
            alignas(A::alignment()) static const uint64_t lut64[][4] = {
                { 0x0000000000000000ul, 0x0000000000000000ul, 0x0000000000000000ul, 0x0000000000000000ul },
                { 0xFFFFFFFFFFFFFFFFul, 0x0000000000000000ul, 0x0000000000000000ul, 0x0000000000000000ul },
                { 0x0000000000000000ul, 0xFFFFFFFFFFFFFFFFul, 0x0000000000000000ul, 0x0000000000000000ul },
                { 0xFFFFFFFFFFFFFFFFul, 0xFFFFFFFFFFFFFFFFul, 0x0000000000000000ul, 0x0000000000000000ul },
                { 0x0000000000000000ul, 0x0000000000000000ul, 0xFFFFFFFFFFFFFFFFul, 0x0000000000000000ul },
                { 0xFFFFFFFFFFFFFFFFul, 0x0000000000000000ul, 0xFFFFFFFFFFFFFFFFul, 0x0000000000000000ul },
                { 0x0000000000000000ul, 0xFFFFFFFFFFFFFFFFul, 0xFFFFFFFFFFFFFFFFul, 0x0000000000000000ul },
                { 0xFFFFFFFFFFFFFFFFul, 0xFFFFFFFFFFFFFFFFul, 0xFFFFFFFFFFFFFFFFul, 0x0000000000000000ul },
                { 0x0000000000000000ul, 0x0000000000000000ul, 0x0000000000000000ul, 0xFFFFFFFFFFFFFFFFul },
                { 0xFFFFFFFFFFFFFFFFul, 0x0000000000000000ul, 0x0000000000000000ul, 0xFFFFFFFFFFFFFFFFul },
                { 0x0000000000000000ul, 0xFFFFFFFFFFFFFFFFul, 0x0000000000000000ul, 0xFFFFFFFFFFFFFFFFul },
                { 0xFFFFFFFFFFFFFFFFul, 0xFFFFFFFFFFFFFFFFul, 0x0000000000000000ul, 0xFFFFFFFFFFFFFFFFul },
                { 0x0000000000000000ul, 0x0000000000000000ul, 0xFFFFFFFFFFFFFFFFul, 0xFFFFFFFFFFFFFFFFul },
                { 0xFFFFFFFFFFFFFFFFul, 0x0000000000000000ul, 0xFFFFFFFFFFFFFFFFul, 0xFFFFFFFFFFFFFFFFul },
                { 0x0000000000000000ul, 0xFFFFFFFFFFFFFFFFul, 0xFFFFFFFFFFFFFFFFul, 0xFFFFFFFFFFFFFFFFul },
                { 0xFFFFFFFFFFFFFFFFul, 0xFFFFFFFFFFFFFFFFul, 0xFFFFFFFFFFFFFFFFul, 0xFFFFFFFFFFFFFFFFul },
            };
            assert(!(mask & ~0xFul) && "inbound mask");
            return _mm256_castsi256_pd(_mm256_load_si256((const __m256i*)lut64[mask]));
        }
        template <class T, class A, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch_bool<T, A> from_mask(batch_bool<T, A> const&, uint64_t mask, requires_arch<avx>) noexcept
        {
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
            alignas(A::alignment()) static const uint64_t lut64[] = {
                0x0000000000000000ul,
                0x000000000000FFFFul,
                0x00000000FFFF0000ul,
                0x00000000FFFFFFFFul,
                0x0000FFFF00000000ul,
                0x0000FFFF0000FFFFul,
                0x0000FFFFFFFF0000ul,
                0x0000FFFFFFFFFFFFul,
                0xFFFF000000000000ul,
                0xFFFF00000000FFFFul,
                0xFFFF0000FFFF0000ul,
                0xFFFF0000FFFFFFFFul,
                0xFFFFFFFF00000000ul,
                0xFFFFFFFF0000FFFFul,
                0xFFFFFFFFFFFF0000ul,
                0xFFFFFFFFFFFFFFFFul,
            };
            XSIMD_IF_CONSTEXPR(sizeof(T) == 1)
            {
                assert(!(mask & ~0xFFFFFFFFul) && "inbound mask");
                return _mm256_setr_epi32(lut32[mask & 0xF], lut32[(mask >> 4) & 0xF],
                                         lut32[(mask >> 8) & 0xF], lut32[(mask >> 12) & 0xF],
                                         lut32[(mask >> 16) & 0xF], lut32[(mask >> 20) & 0xF],
                                         lut32[(mask >> 24) & 0xF], lut32[mask >> 28]);
            }
            else XSIMD_IF_CONSTEXPR(sizeof(T) == 2)
            {
                assert(!(mask & ~0xFFFFul) && "inbound mask");
                return _mm256_setr_epi64x(lut64[mask & 0xF], lut64[(mask >> 4) & 0xF], lut64[(mask >> 8) & 0xF], lut64[(mask >> 12) & 0xF]);
            }
            else XSIMD_IF_CONSTEXPR(sizeof(T) == 4)
            {
                return _mm256_castps_si256(from_mask(batch_bool<float, A> {}, mask, avx {}));
            }
            else XSIMD_IF_CONSTEXPR(sizeof(T) == 8)
            {
                return _mm256_castpd_si256(from_mask(batch_bool<double, A> {}, mask, avx {}));
            }
        }

        // haddp
        template <class A>
        XSIMD_INLINE batch<float, A> haddp(batch<float, A> const* row, requires_arch<avx>) noexcept
        {
            // row = (a,b,c,d,e,f,g,h)
            // tmp0 = (a0+a1, a2+a3, b0+b1, b2+b3, a4+a5, a6+a7, b4+b5, b6+b7)
            __m256 tmp0 = _mm256_hadd_ps(row[0], row[1]);
            // tmp1 = (c0+c1, c2+c3, d1+d2, d2+d3, c4+c5, c6+c7, d4+d5, d6+d7)
            __m256 tmp1 = _mm256_hadd_ps(row[2], row[3]);
            // tmp1 = (a0+a1+a2+a3, b0+b1+b2+b3, c0+c1+c2+c3, d0+d1+d2+d3,
            // a4+a5+a6+a7, b4+b5+b6+b7, c4+c5+c6+c7, d4+d5+d6+d7)
            tmp1 = _mm256_hadd_ps(tmp0, tmp1);
            // tmp0 = (e0+e1, e2+e3, f0+f1, f2+f3, e4+e5, e6+e7, f4+f5, f6+f7)
            tmp0 = _mm256_hadd_ps(row[4], row[5]);
            // tmp2 = (g0+g1, g2+g3, h0+h1, h2+h3, g4+g5, g6+g7, h4+h5, h6+h7)
            __m256 tmp2 = _mm256_hadd_ps(row[6], row[7]);
            // tmp2 = (e0+e1+e2+e3, f0+f1+f2+f3, g0+g1+g2+g3, h0+h1+h2+h3,
            // e4+e5+e6+e7, f4+f5+f6+f7, g4+g5+g6+g7, h4+h5+h6+h7)
            tmp2 = _mm256_hadd_ps(tmp0, tmp2);
            // tmp0 = (a0+a1+a2+a3, b0+b1+b2+b3, c0+c1+c2+c3, d0+d1+d2+d3,
            // e4+e5+e6+e7, f4+f5+f6+f7, g4+g5+g6+g7, h4+h5+h6+h7)
            tmp0 = _mm256_blend_ps(tmp1, tmp2, 0b11110000);
            // tmp1 = (a4+a5+a6+a7, b4+b5+b6+b7, c4+c5+c6+c7, d4+d5+d6+d7,
            // e0+e1+e2+e3, f0+f1+f2+f3, g0+g1+g2+g3, h0+h1+h2+h3)
            tmp1 = _mm256_permute2f128_ps(tmp1, tmp2, 0x21);
            return _mm256_add_ps(tmp0, tmp1);
        }
        template <class A>
        XSIMD_INLINE batch<double, A> haddp(batch<double, A> const* row, requires_arch<avx>) noexcept
        {
            // row = (a,b,c,d)
            // tmp0 = (a0+a1, b0+b1, a2+a3, b2+b3)
            __m256d tmp0 = _mm256_hadd_pd(row[0], row[1]);
            // tmp1 = (c0+c1, d0+d1, c2+c3, d2+d3)
            __m256d tmp1 = _mm256_hadd_pd(row[2], row[3]);
            // tmp2 = (a0+a1, b0+b1, c2+c3, d2+d3)
            __m256d tmp2 = _mm256_blend_pd(tmp0, tmp1, 0b1100);
            // tmp1 = (a2+a3, b2+b3, c2+c3, d2+d3)
            tmp1 = _mm256_permute2f128_pd(tmp0, tmp1, 0x21);
            return _mm256_add_pd(tmp1, tmp2);
        }

        // incr_if
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> incr_if(batch<T, A> const& self, batch_bool<T, A> const& mask, requires_arch<avx>) noexcept
        {
            return self - batch<T, A>(mask.data);
        }

        // insert
        template <class A, class T, size_t I, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> insert(batch<T, A> const& self, T val, index<I> pos, requires_arch<avx>) noexcept
        {
#if !defined(_MSC_VER) || _MSC_VER > 1900
            XSIMD_IF_CONSTEXPR(sizeof(T) == 1)
            {
                return _mm256_insert_epi8(self, val, I);
            }
            else XSIMD_IF_CONSTEXPR(sizeof(T) == 2)
            {
                return _mm256_insert_epi16(self, val, I);
            }
            else XSIMD_IF_CONSTEXPR(sizeof(T) == 4)
            {
                return _mm256_insert_epi32(self, val, I);
            }
            else
            {
                return insert(self, val, pos, generic {});
            }
#endif
            return insert(self, val, pos, generic {});
        }

        // isnan
        template <class A>
        XSIMD_INLINE batch_bool<float, A> isnan(batch<float, A> const& self, requires_arch<avx>) noexcept
        {
            return _mm256_cmp_ps(self, self, _CMP_UNORD_Q);
        }
        template <class A>
        XSIMD_INLINE batch_bool<double, A> isnan(batch<double, A> const& self, requires_arch<avx>) noexcept
        {
            return _mm256_cmp_pd(self, self, _CMP_UNORD_Q);
        }

        // le
        template <class A>
        XSIMD_INLINE batch_bool<float, A> le(batch<float, A> const& self, batch<float, A> const& other, requires_arch<avx>) noexcept
        {
            return _mm256_cmp_ps(self, other, _CMP_LE_OQ);
        }
        template <class A>
        XSIMD_INLINE batch_bool<double, A> le(batch<double, A> const& self, batch<double, A> const& other, requires_arch<avx>) noexcept
        {
            return _mm256_cmp_pd(self, other, _CMP_LE_OQ);
        }

        // load_aligned
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> load_aligned(T const* mem, convert<T>, requires_arch<avx>) noexcept
        {
            return _mm256_load_si256((__m256i const*)mem);
        }
        template <class A>
        XSIMD_INLINE batch<float, A> load_aligned(float const* mem, convert<float>, requires_arch<avx>) noexcept
        {
            return _mm256_load_ps(mem);
        }
        template <class A>
        XSIMD_INLINE batch<double, A> load_aligned(double const* mem, convert<double>, requires_arch<avx>) noexcept
        {
            return _mm256_load_pd(mem);
        }

        namespace detail
        {
            // load_complex
            template <class A>
            XSIMD_INLINE batch<std::complex<float>, A> load_complex(batch<float, A> const& hi, batch<float, A> const& lo, requires_arch<avx>) noexcept
            {
                using batch_type = batch<float, A>;
                __m128 tmp0 = _mm256_extractf128_ps(hi, 0);
                __m128 tmp1 = _mm256_extractf128_ps(hi, 1);
                __m128 tmp_real = _mm_shuffle_ps(tmp0, tmp1, _MM_SHUFFLE(2, 0, 2, 0));
                __m128 tmp_imag = _mm_shuffle_ps(tmp0, tmp1, _MM_SHUFFLE(3, 1, 3, 1));
                batch_type real = _mm256_castps128_ps256(tmp_real);
                batch_type imag = _mm256_castps128_ps256(tmp_imag);

                tmp0 = _mm256_extractf128_ps(lo, 0);
                tmp1 = _mm256_extractf128_ps(lo, 1);
                tmp_real = _mm_shuffle_ps(tmp0, tmp1, _MM_SHUFFLE(2, 0, 2, 0));
                tmp_imag = _mm_shuffle_ps(tmp0, tmp1, _MM_SHUFFLE(3, 1, 3, 1));
                real = _mm256_insertf128_ps(real, tmp_real, 1);
                imag = _mm256_insertf128_ps(imag, tmp_imag, 1);
                return { real, imag };
            }
            template <class A>
            XSIMD_INLINE batch<std::complex<double>, A> load_complex(batch<double, A> const& hi, batch<double, A> const& lo, requires_arch<avx>) noexcept
            {
                using batch_type = batch<double, A>;
                __m128d tmp0 = _mm256_extractf128_pd(hi, 0);
                __m128d tmp1 = _mm256_extractf128_pd(hi, 1);
                batch_type real = _mm256_castpd128_pd256(_mm_unpacklo_pd(tmp0, tmp1));
                batch_type imag = _mm256_castpd128_pd256(_mm_unpackhi_pd(tmp0, tmp1));

                tmp0 = _mm256_extractf128_pd(lo, 0);
                tmp1 = _mm256_extractf128_pd(lo, 1);
                __m256d re_tmp1 = _mm256_insertf128_pd(real, _mm_unpacklo_pd(tmp0, tmp1), 1);
                __m256d im_tmp1 = _mm256_insertf128_pd(imag, _mm_unpackhi_pd(tmp0, tmp1), 1);
                real = _mm256_blend_pd(real, re_tmp1, 12);
                imag = _mm256_blend_pd(imag, im_tmp1, 12);
                return { real, imag };
            }
        }

        // load_unaligned
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> load_unaligned(T const* mem, convert<T>, requires_arch<avx>) noexcept
        {
            return _mm256_loadu_si256((__m256i const*)mem);
        }
        template <class A>
        XSIMD_INLINE batch<float, A> load_unaligned(float const* mem, convert<float>, requires_arch<avx>) noexcept
        {
            return _mm256_loadu_ps(mem);
        }
        template <class A>
        XSIMD_INLINE batch<double, A> load_unaligned(double const* mem, convert<double>, requires_arch<avx>) noexcept
        {
            return _mm256_loadu_pd(mem);
        }

        // lt
        template <class A>
        XSIMD_INLINE batch_bool<float, A> lt(batch<float, A> const& self, batch<float, A> const& other, requires_arch<avx>) noexcept
        {
            return _mm256_cmp_ps(self, other, _CMP_LT_OQ);
        }
        template <class A>
        XSIMD_INLINE batch_bool<double, A> lt(batch<double, A> const& self, batch<double, A> const& other, requires_arch<avx>) noexcept
        {
            return _mm256_cmp_pd(self, other, _CMP_LT_OQ);
        }

        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch_bool<T, A> lt(batch<T, A> const& self, batch<T, A> const& other, requires_arch<avx>) noexcept
        {
            return detail::fwd_to_sse([](__m128i s, __m128i o) noexcept
                                      { return lt(batch<T, sse4_2>(s), batch<T, sse4_2>(o)); },
                                      self, other);
        }

        // mask
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE uint64_t mask(batch_bool<T, A> const& self, requires_arch<avx>) noexcept
        {
            XSIMD_IF_CONSTEXPR(sizeof(T) == 1 || sizeof(T) == 2)
            {
                __m128i self_low, self_high;
                detail::split_avx(self, self_low, self_high);
                return mask(batch_bool<T, sse4_2>(self_low), sse4_2 {}) | (mask(batch_bool<T, sse4_2>(self_high), sse4_2 {}) << (128 / (8 * sizeof(T))));
            }
            else XSIMD_IF_CONSTEXPR(sizeof(T) == 4)
            {
                return _mm256_movemask_ps(_mm256_castsi256_ps(self));
            }
            else XSIMD_IF_CONSTEXPR(sizeof(T) == 8)
            {
                return _mm256_movemask_pd(_mm256_castsi256_pd(self));
            }
            else
            {
                assert(false && "unsupported arch/op combination");
                return {};
            }
        }
        template <class A>
        XSIMD_INLINE uint64_t mask(batch_bool<float, A> const& self, requires_arch<avx>) noexcept
        {
            return _mm256_movemask_ps(self);
        }

        template <class A>
        XSIMD_INLINE uint64_t mask(batch_bool<double, A> const& self, requires_arch<avx>) noexcept
        {
            return _mm256_movemask_pd(self);
        }

        // max
        template <class A>
        XSIMD_INLINE batch<float, A> max(batch<float, A> const& self, batch<float, A> const& other, requires_arch<avx>) noexcept
        {
            return _mm256_max_ps(self, other);
        }
        template <class A>
        XSIMD_INLINE batch<double, A> max(batch<double, A> const& self, batch<double, A> const& other, requires_arch<avx>) noexcept
        {
            return _mm256_max_pd(self, other);
        }
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> max(batch<T, A> const& self, batch<T, A> const& other, requires_arch<avx>) noexcept
        {
            return select(self > other, self, other);
        }

        // min
        template <class A>
        XSIMD_INLINE batch<float, A> min(batch<float, A> const& self, batch<float, A> const& other, requires_arch<avx>) noexcept
        {
            return _mm256_min_ps(self, other);
        }
        template <class A>
        XSIMD_INLINE batch<double, A> min(batch<double, A> const& self, batch<double, A> const& other, requires_arch<avx>) noexcept
        {
            return _mm256_min_pd(self, other);
        }
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> min(batch<T, A> const& self, batch<T, A> const& other, requires_arch<avx>) noexcept
        {
            return select(self <= other, self, other);
        }

        // mul
        template <class A>
        XSIMD_INLINE batch<float, A> mul(batch<float, A> const& self, batch<float, A> const& other, requires_arch<avx>) noexcept
        {
            return _mm256_mul_ps(self, other);
        }
        template <class A>
        XSIMD_INLINE batch<double, A> mul(batch<double, A> const& self, batch<double, A> const& other, requires_arch<avx>) noexcept
        {
            return _mm256_mul_pd(self, other);
        }

        // nearbyint
        template <class A>
        XSIMD_INLINE batch<float, A> nearbyint(batch<float, A> const& self, requires_arch<avx>) noexcept
        {
            return _mm256_round_ps(self, _MM_FROUND_TO_NEAREST_INT);
        }
        template <class A>
        XSIMD_INLINE batch<double, A> nearbyint(batch<double, A> const& self, requires_arch<avx>) noexcept
        {
            return _mm256_round_pd(self, _MM_FROUND_TO_NEAREST_INT);
        }

        // nearbyint_as_int
        template <class A>
        XSIMD_INLINE batch<int32_t, A> nearbyint_as_int(batch<float, A> const& self,
                                                        requires_arch<avx>) noexcept
        {
            return _mm256_cvtps_epi32(self);
        }

        // neg
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> neg(batch<T, A> const& self, requires_arch<avx>) noexcept
        {
            return 0 - self;
        }
        template <class A>
        batch<float, A> neg(batch<float, A> const& self, requires_arch<avx>)
        {
            return _mm256_xor_ps(self, _mm256_castsi256_ps(_mm256_set1_epi32(0x80000000)));
        }
        template <class A>
        XSIMD_INLINE batch<double, A> neg(batch<double, A> const& self, requires_arch<avx>) noexcept
        {
            return _mm256_xor_pd(self, _mm256_castsi256_pd(_mm256_set1_epi64x(0x8000000000000000)));
        }

        // neq
        template <class A>
        XSIMD_INLINE batch_bool<float, A> neq(batch<float, A> const& self, batch<float, A> const& other, requires_arch<avx>) noexcept
        {
            return _mm256_cmp_ps(self, other, _CMP_NEQ_UQ);
        }
        template <class A>
        XSIMD_INLINE batch_bool<double, A> neq(batch<double, A> const& self, batch<double, A> const& other, requires_arch<avx>) noexcept
        {
            return _mm256_cmp_pd(self, other, _CMP_NEQ_UQ);
        }
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch_bool<T, A> neq(batch<T, A> const& self, batch<T, A> const& other, requires_arch<avx>) noexcept
        {
            return ~(self == other);
        }

        template <class A>
        XSIMD_INLINE batch_bool<float, A> neq(batch_bool<float, A> const& self, batch_bool<float, A> const& other, requires_arch<avx>) noexcept
        {
            return _mm256_xor_ps(self, other);
        }
        template <class A>
        XSIMD_INLINE batch_bool<double, A> neq(batch_bool<double, A> const& self, batch_bool<double, A> const& other, requires_arch<avx>) noexcept
        {
            return _mm256_xor_pd(self, other);
        }
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch_bool<T, A> neq(batch_bool<T, A> const& self, batch_bool<T, A> const& other, requires_arch<avx>) noexcept
        {
            return _mm256_castps_si256(_mm256_xor_ps(_mm256_castsi256_ps(self.data), _mm256_castsi256_ps(other.data)));
        }

        // reciprocal
        template <class A>
        XSIMD_INLINE batch<float, A> reciprocal(batch<float, A> const& self,
                                                kernel::requires_arch<avx>) noexcept
        {
            return _mm256_rcp_ps(self);
        }

        // reduce_add
        template <class A>
        XSIMD_INLINE float reduce_add(batch<float, A> const& rhs, requires_arch<avx>) noexcept
        {
            // Warning about _mm256_hadd_ps:
            // _mm256_hadd_ps(a,b) gives
            // (a0+a1,a2+a3,b0+b1,b2+b3,a4+a5,a6+a7,b4+b5,b6+b7). Hence we can't
            // rely on a naive use of this method
            // rhs = (x0, x1, x2, x3, x4, x5, x6, x7)
            // tmp = (x4, x5, x6, x7, x0, x1, x2, x3)
            __m256 tmp = _mm256_permute2f128_ps(rhs, rhs, 1);
            // tmp = (x4+x0, x5+x1, x6+x2, x7+x3, x0+x4, x1+x5, x2+x6, x3+x7)
            tmp = _mm256_add_ps(rhs, tmp);
            // tmp = (x4+x0+x5+x1, x6+x2+x7+x3, -, -, -, -, -, -)
            tmp = _mm256_hadd_ps(tmp, tmp);
            // tmp = (x4+x0+x5+x1+x6+x2+x7+x3, -, -, -, -, -, -, -)
            tmp = _mm256_hadd_ps(tmp, tmp);
            return _mm_cvtss_f32(_mm256_extractf128_ps(tmp, 0));
        }
        template <class A>
        XSIMD_INLINE double reduce_add(batch<double, A> const& rhs, requires_arch<avx>) noexcept
        {
            // rhs = (x0, x1, x2, x3)
            // tmp = (x2, x3, x0, x1)
            __m256d tmp = _mm256_permute2f128_pd(rhs, rhs, 1);
            // tmp = (x2+x0, x3+x1, -, -)
            tmp = _mm256_add_pd(rhs, tmp);
            // tmp = (x2+x0+x3+x1, -, -, -)
            tmp = _mm256_hadd_pd(tmp, tmp);
            return _mm_cvtsd_f64(_mm256_extractf128_pd(tmp, 0));
        }
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE T reduce_add(batch<T, A> const& self, requires_arch<avx>) noexcept
        {
            __m128i low, high;
            detail::split_avx(self, low, high);
            batch<T, sse4_2> blow(low), bhigh(high);
            return reduce_add(blow) + reduce_add(bhigh);
        }

        // reduce_max
        template <class A, class T, class _ = typename std::enable_if<(sizeof(T) <= 2), void>::type>
        XSIMD_INLINE T reduce_max(batch<T, A> const& self, requires_arch<avx>) noexcept
        {
            constexpr auto mask = detail::shuffle(1, 0);
            batch<T, A> step = _mm256_permute2f128_si256(self, self, mask);
            batch<T, A> acc = max(self, step);
            __m128i low = _mm256_castsi256_si128(acc);
            return reduce_max(batch<T, sse4_2>(low));
        }

        // reduce_min
        template <class A, class T, class _ = typename std::enable_if<(sizeof(T) <= 2), void>::type>
        XSIMD_INLINE T reduce_min(batch<T, A> const& self, requires_arch<avx>) noexcept
        {
            constexpr auto mask = detail::shuffle(1, 0);
            batch<T, A> step = _mm256_permute2f128_si256(self, self, mask);
            batch<T, A> acc = min(self, step);
            __m128i low = _mm256_castsi256_si128(acc);
            return reduce_min(batch<T, sse4_2>(low));
        }

        // rsqrt
        template <class A>
        XSIMD_INLINE batch<float, A> rsqrt(batch<float, A> const& val, requires_arch<avx>) noexcept
        {
            return _mm256_rsqrt_ps(val);
        }
        template <class A>
        XSIMD_INLINE batch<double, A> rsqrt(batch<double, A> const& val, requires_arch<avx>) noexcept
        {
            return _mm256_cvtps_pd(_mm_rsqrt_ps(_mm256_cvtpd_ps(val)));
        }

        // sadd
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> sadd(batch<T, A> const& self, batch<T, A> const& other, requires_arch<avx>) noexcept
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

        // select
        template <class A>
        XSIMD_INLINE batch<float, A> select(batch_bool<float, A> const& cond, batch<float, A> const& true_br, batch<float, A> const& false_br, requires_arch<avx>) noexcept
        {
            return _mm256_blendv_ps(false_br, true_br, cond);
        }
        template <class A>
        XSIMD_INLINE batch<double, A> select(batch_bool<double, A> const& cond, batch<double, A> const& true_br, batch<double, A> const& false_br, requires_arch<avx>) noexcept
        {
            return _mm256_blendv_pd(false_br, true_br, cond);
        }
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> select(batch_bool<T, A> const& cond, batch<T, A> const& true_br, batch<T, A> const& false_br, requires_arch<avx>) noexcept
        {
            __m128i cond_low, cond_hi;
            detail::split_avx(cond, cond_low, cond_hi);

            __m128i true_low, true_hi;
            detail::split_avx(true_br, true_low, true_hi);

            __m128i false_low, false_hi;
            detail::split_avx(false_br, false_low, false_hi);

            __m128i res_low = select(batch_bool<T, sse4_2>(cond_low), batch<T, sse4_2>(true_low), batch<T, sse4_2>(false_low), sse4_2 {});
            __m128i res_hi = select(batch_bool<T, sse4_2>(cond_hi), batch<T, sse4_2>(true_hi), batch<T, sse4_2>(false_hi), sse4_2 {});
            return detail::merge_sse(res_low, res_hi);
        }
        template <class A, class T, bool... Values, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> select(batch_bool_constant<T, A, Values...> const&, batch<T, A> const& true_br, batch<T, A> const& false_br, requires_arch<avx>) noexcept
        {
            return select(batch_bool<T, A> { Values... }, true_br, false_br, avx2 {});
        }

        template <class A, bool... Values>
        XSIMD_INLINE batch<float, A> select(batch_bool_constant<float, A, Values...> const&, batch<float, A> const& true_br, batch<float, A> const& false_br, requires_arch<avx>) noexcept
        {
            constexpr auto mask = batch_bool_constant<float, A, Values...>::mask();
            return _mm256_blend_ps(false_br, true_br, mask);
        }

        template <class A, bool... Values>
        XSIMD_INLINE batch<double, A> select(batch_bool_constant<double, A, Values...> const&, batch<double, A> const& true_br, batch<double, A> const& false_br, requires_arch<avx>) noexcept
        {
            constexpr auto mask = batch_bool_constant<double, A, Values...>::mask();
            return _mm256_blend_pd(false_br, true_br, mask);
        }

        // set
        template <class A, class... Values>
        XSIMD_INLINE batch<float, A> set(batch<float, A> const&, requires_arch<avx>, Values... values) noexcept
        {
            static_assert(sizeof...(Values) == batch<float, A>::size, "consistent init");
            return _mm256_setr_ps(values...);
        }

        template <class A, class... Values>
        XSIMD_INLINE batch<double, A> set(batch<double, A> const&, requires_arch<avx>, Values... values) noexcept
        {
            static_assert(sizeof...(Values) == batch<double, A>::size, "consistent init");
            return _mm256_setr_pd(values...);
        }
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> set(batch<T, A> const&, requires_arch<avx>, T v0, T v1, T v2, T v3) noexcept
        {
            return _mm256_set_epi64x(v3, v2, v1, v0);
        }
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> set(batch<T, A> const&, requires_arch<avx>, T v0, T v1, T v2, T v3, T v4, T v5, T v6, T v7) noexcept
        {
            return _mm256_setr_epi32(v0, v1, v2, v3, v4, v5, v6, v7);
        }
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> set(batch<T, A> const&, requires_arch<avx>, T v0, T v1, T v2, T v3, T v4, T v5, T v6, T v7, T v8, T v9, T v10, T v11, T v12, T v13, T v14, T v15) noexcept
        {
            return _mm256_setr_epi16(v0, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, v15);
        }
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> set(batch<T, A> const&, requires_arch<avx>, T v0, T v1, T v2, T v3, T v4, T v5, T v6, T v7, T v8, T v9, T v10, T v11, T v12, T v13, T v14, T v15,
                                     T v16, T v17, T v18, T v19, T v20, T v21, T v22, T v23, T v24, T v25, T v26, T v27, T v28, T v29, T v30, T v31) noexcept
        {
            return _mm256_setr_epi8(v0, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, v15, v16, v17, v18, v19, v20, v21, v22, v23, v24, v25, v26, v27, v28, v29, v30, v31);
        }

        template <class A, class T, class... Values, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch_bool<T, A> set(batch_bool<T, A> const&, requires_arch<avx>, Values... values) noexcept
        {
            return set(batch<T, A>(), A {}, static_cast<T>(values ? -1LL : 0LL)...).data;
        }

        template <class A, class... Values>
        XSIMD_INLINE batch_bool<float, A> set(batch_bool<float, A> const&, requires_arch<avx>, Values... values) noexcept
        {
            static_assert(sizeof...(Values) == batch_bool<float, A>::size, "consistent init");
            return _mm256_castsi256_ps(set(batch<int32_t, A>(), A {}, static_cast<int32_t>(values ? -1LL : 0LL)...).data);
        }

        template <class A, class... Values>
        XSIMD_INLINE batch_bool<double, A> set(batch_bool<double, A> const&, requires_arch<avx>, Values... values) noexcept
        {
            static_assert(sizeof...(Values) == batch_bool<double, A>::size, "consistent init");
            return _mm256_castsi256_pd(set(batch<int64_t, A>(), A {}, static_cast<int64_t>(values ? -1LL : 0LL)...).data);
        }

        // shuffle
        template <class A, class ITy, ITy I0, ITy I1, ITy I2, ITy I3, ITy I4, ITy I5, ITy I6, ITy I7>
        XSIMD_INLINE batch<float, A> shuffle(batch<float, A> const& x, batch<float, A> const& y, batch_constant<ITy, A, I0, I1, I2, I3, I4, I5, I6, I7> mask, requires_arch<avx>) noexcept
        {
            constexpr uint32_t smask = detail::mod_shuffle(I0, I1, I2, I3);
            // shuffle within lane
            if (I4 == (I0 + 4) && I5 == (I1 + 4) && I6 == (I2 + 4) && I7 == (I3 + 4) && I0 < 4 && I1 < 4 && I2 >= 8 && I2 < 12 && I3 >= 8 && I3 < 12)
                return _mm256_shuffle_ps(x, y, smask);

            // shuffle within opposite lane
            if (I4 == (I0 + 4) && I5 == (I1 + 4) && I6 == (I2 + 4) && I7 == (I3 + 4) && I2 < 4 && I3 < 4 && I0 >= 8 && I0 < 12 && I1 >= 8 && I1 < 12)
                return _mm256_shuffle_ps(y, x, smask);

            return shuffle(x, y, mask, generic {});
        }

        template <class A, class ITy, ITy I0, ITy I1, ITy I2, ITy I3>
        XSIMD_INLINE batch<double, A> shuffle(batch<double, A> const& x, batch<double, A> const& y, batch_constant<ITy, A, I0, I1, I2, I3> mask, requires_arch<avx>) noexcept
        {
            constexpr uint32_t smask = (I0 & 0x1) | ((I1 & 0x1) << 1) | ((I2 & 0x1) << 2) | ((I3 & 0x1) << 3);
            // shuffle within lane
            if (I0 < 2 && I1 >= 4 && I1 < 6 && I2 >= 2 && I2 < 4 && I3 >= 6)
                return _mm256_shuffle_pd(x, y, smask);

            // shuffle within opposite lane
            if (I1 < 2 && I0 >= 4 && I0 < 6 && I3 >= 2 && I3 < 4 && I2 >= 6)
                return _mm256_shuffle_pd(y, x, smask);

            return shuffle(x, y, mask, generic {});
        }

        // slide_left
        template <size_t N, class A, class T>
        XSIMD_INLINE batch<T, A> slide_left(batch<T, A> const& x, requires_arch<avx>) noexcept
        {
            constexpr unsigned BitCount = N * 8;
            if (BitCount == 0)
            {
                return x;
            }
            if (BitCount >= 256)
            {
                return batch<T, A>(T(0));
            }
            if (BitCount > 128)
            {
                constexpr unsigned M = (BitCount - 128) / 8;
                __m128i low = _mm256_castsi256_si128(x);
                auto y = _mm_slli_si128(low, M);
                __m256i zero = _mm256_setzero_si256();
                return _mm256_insertf128_si256(zero, y, 1);
            }
            if (BitCount == 128)
            {
                __m128i low = _mm256_castsi256_si128(x);
                __m256i zero = _mm256_setzero_si256();
                return _mm256_insertf128_si256(zero, low, 1);
            }
            // shifting by [0, 128[ bits
            constexpr unsigned M = BitCount / 8;

            __m128i low = _mm256_castsi256_si128(x);
            auto ylow = _mm_slli_si128(low, M);
            auto zlow = _mm_srli_si128(low, 16 - M);

            __m128i high = _mm256_extractf128_si256(x, 1);
            auto yhigh = _mm_slli_si128(high, M);

            __m256i res = _mm256_castsi128_si256(ylow);
            return _mm256_insertf128_si256(res, _mm_or_si128(yhigh, zlow), 1);
        }

        // slide_right
        template <size_t N, class A, class T>
        XSIMD_INLINE batch<T, A> slide_right(batch<T, A> const& x, requires_arch<avx>) noexcept
        {
            constexpr unsigned BitCount = N * 8;
            if (BitCount == 0)
            {
                return x;
            }
            if (BitCount >= 256)
            {
                return batch<T, A>(T(0));
            }
            if (BitCount > 128)
            {
                constexpr unsigned M = (BitCount - 128) / 8;
                __m128i high = _mm256_extractf128_si256(x, 1);
                __m128i y = _mm_srli_si128(high, M);
                __m256i zero = _mm256_setzero_si256();
                return _mm256_insertf128_si256(zero, y, 0);
            }
            if (BitCount == 128)
            {
                __m128i high = _mm256_extractf128_si256(x, 1);
                return _mm256_castsi128_si256(high);
            }
            // shifting by [0, 128[ bits
            constexpr unsigned M = BitCount / 8;

            __m128i low = _mm256_castsi256_si128(x);
            auto ylow = _mm_srli_si128(low, M);

            __m128i high = _mm256_extractf128_si256(x, 1);
            auto yhigh = _mm_srli_si128(high, M);
            auto zhigh = _mm_slli_si128(high, 16 - M);

            __m256i res = _mm256_castsi128_si256(_mm_or_si128(ylow, zhigh));
            return _mm256_insertf128_si256(res, yhigh, 1);
        }

        // sqrt
        template <class A>
        XSIMD_INLINE batch<float, A> sqrt(batch<float, A> const& val, requires_arch<avx>) noexcept
        {
            return _mm256_sqrt_ps(val);
        }
        template <class A>
        XSIMD_INLINE batch<double, A> sqrt(batch<double, A> const& val, requires_arch<avx>) noexcept
        {
            return _mm256_sqrt_pd(val);
        }

        // ssub
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> ssub(batch<T, A> const& self, batch<T, A> const& other, requires_arch<avx>) noexcept
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

        // store_aligned
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE void store_aligned(T* mem, batch<T, A> const& self, requires_arch<avx>) noexcept
        {
            return _mm256_store_si256((__m256i*)mem, self);
        }
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE void store_aligned(T* mem, batch_bool<T, A> const& self, requires_arch<avx>) noexcept
        {
            return _mm256_store_si256((__m256i*)mem, self);
        }
        template <class A>
        XSIMD_INLINE void store_aligned(float* mem, batch<float, A> const& self, requires_arch<avx>) noexcept
        {
            return _mm256_store_ps(mem, self);
        }
        template <class A>
        XSIMD_INLINE void store_aligned(double* mem, batch<double, A> const& self, requires_arch<avx>) noexcept
        {
            return _mm256_store_pd(mem, self);
        }

        // store_unaligned
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE void store_unaligned(T* mem, batch<T, A> const& self, requires_arch<avx>) noexcept
        {
            return _mm256_storeu_si256((__m256i*)mem, self);
        }
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE void store_unaligned(T* mem, batch_bool<T, A> const& self, requires_arch<avx>) noexcept
        {
            return _mm256_storeu_si256((__m256i*)mem, self);
        }
        template <class A>
        XSIMD_INLINE void store_unaligned(float* mem, batch<float, A> const& self, requires_arch<avx>) noexcept
        {
            return _mm256_storeu_ps(mem, self);
        }
        template <class A>
        XSIMD_INLINE void store_unaligned(double* mem, batch<double, A> const& self, requires_arch<avx>) noexcept
        {
            return _mm256_storeu_pd(mem, self);
        }

        // sub
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> sub(batch<T, A> const& self, batch<T, A> const& other, requires_arch<avx>) noexcept
        {
            return detail::fwd_to_sse([](__m128i s, __m128i o) noexcept
                                      { return sub(batch<T, sse4_2>(s), batch<T, sse4_2>(o)); },
                                      self, other);
        }
        template <class A>
        XSIMD_INLINE batch<float, A> sub(batch<float, A> const& self, batch<float, A> const& other, requires_arch<avx>) noexcept
        {
            return _mm256_sub_ps(self, other);
        }
        template <class A>
        XSIMD_INLINE batch<double, A> sub(batch<double, A> const& self, batch<double, A> const& other, requires_arch<avx>) noexcept
        {
            return _mm256_sub_pd(self, other);
        }

        // swizzle (dynamic mask)
        template <class A>
        XSIMD_INLINE batch<float, A> swizzle(batch<float, A> const& self, batch<uint32_t, A> mask, requires_arch<avx>) noexcept
        {
            // duplicate low and high part of input
            __m256 hi = _mm256_castps128_ps256(_mm256_extractf128_ps(self, 1));
            __m256 hi_hi = _mm256_insertf128_ps(self, _mm256_castps256_ps128(hi), 0);

            __m256 low = _mm256_castps128_ps256(_mm256_castps256_ps128(self));
            __m256 low_low = _mm256_insertf128_ps(self, _mm256_castps256_ps128(low), 1);

            // normalize mask
            batch<uint32_t, A> half_mask = mask % 4;

            // permute within each lane
            __m256 r0 = _mm256_permutevar_ps(low_low, half_mask);
            __m256 r1 = _mm256_permutevar_ps(hi_hi, half_mask);

            // mask to choose the right lane
            batch_bool<uint32_t, A> blend_mask = mask >= 4;

            // blend the two permutes
            return _mm256_blendv_ps(r0, r1, batch_bool_cast<float>(blend_mask));
        }

        template <class A>
        XSIMD_INLINE batch<double, A> swizzle(batch<double, A> const& self, batch<uint64_t, A> mask, requires_arch<avx>) noexcept
        {
            // duplicate low and high part of input
            __m256d hi = _mm256_castpd128_pd256(_mm256_extractf128_pd(self, 1));
            __m256d hi_hi = _mm256_insertf128_pd(self, _mm256_castpd256_pd128(hi), 0);

            __m256d low = _mm256_castpd128_pd256(_mm256_castpd256_pd128(self));
            __m256d low_low = _mm256_insertf128_pd(self, _mm256_castpd256_pd128(low), 1);

            // normalize mask
            batch<uint64_t, A> half_mask = -(mask & 1);

            // permute within each lane
            __m256d r0 = _mm256_permutevar_pd(low_low, half_mask);
            __m256d r1 = _mm256_permutevar_pd(hi_hi, half_mask);

            // mask to choose the right lane
            batch_bool<uint64_t, A> blend_mask = mask >= 2;

            // blend the two permutes
            return _mm256_blendv_pd(r0, r1, batch_bool_cast<double>(blend_mask));
        }

        template <class A, typename T, detail::enable_sized_integral_t<T, 4> = 0>
        XSIMD_INLINE batch<T, A> swizzle(batch<T, A> const& self, batch<uint32_t, A> const& mask, requires_arch<avx>) noexcept
        {
            return bitwise_cast<T>(
                swizzle(bitwise_cast<float>(self), mask));
        }

        template <class A, typename T, detail::enable_sized_integral_t<T, 8> = 0>
        XSIMD_INLINE batch<T, A>
        swizzle(batch<T, A> const& self, batch<uint64_t, A> const& mask, requires_arch<avx>) noexcept
        {
            return bitwise_cast<T>(
                swizzle(bitwise_cast<double>(self), mask));
        }

        // swizzle (constant mask)
        template <class A, uint32_t V0, uint32_t V1, uint32_t V2, uint32_t V3, uint32_t V4, uint32_t V5, uint32_t V6, uint32_t V7>
        XSIMD_INLINE batch<float, A> swizzle(batch<float, A> const& self, batch_constant<uint32_t, A, V0, V1, V2, V3, V4, V5, V6, V7>, requires_arch<avx>) noexcept
        {
            // duplicate low and high part of input
            __m256 hi = _mm256_castps128_ps256(_mm256_extractf128_ps(self, 1));
            __m256 hi_hi = _mm256_insertf128_ps(self, _mm256_castps256_ps128(hi), 0);

            __m256 low = _mm256_castps128_ps256(_mm256_castps256_ps128(self));
            __m256 low_low = _mm256_insertf128_ps(self, _mm256_castps256_ps128(low), 1);

            // normalize mask
            batch_constant<uint32_t, A, (V0 % 4), (V1 % 4), (V2 % 4), (V3 % 4), (V4 % 4), (V5 % 4), (V6 % 4), (V7 % 4)> half_mask;

            // permute within each lane
            __m256 r0 = _mm256_permutevar_ps(low_low, half_mask.as_batch());
            __m256 r1 = _mm256_permutevar_ps(hi_hi, half_mask.as_batch());

            // mask to choose the right lane
            batch_bool_constant<uint32_t, A, (V0 >= 4), (V1 >= 4), (V2 >= 4), (V3 >= 4), (V4 >= 4), (V5 >= 4), (V6 >= 4), (V7 >= 4)> blend_mask;

            // blend the two permutes
            constexpr auto mask = blend_mask.mask();
            return _mm256_blend_ps(r0, r1, mask);
        }

        template <class A, uint64_t V0, uint64_t V1, uint64_t V2, uint64_t V3>
        XSIMD_INLINE batch<double, A> swizzle(batch<double, A> const& self, batch_constant<uint64_t, A, V0, V1, V2, V3>, requires_arch<avx>) noexcept
        {
            // duplicate low and high part of input
            __m256d hi = _mm256_castpd128_pd256(_mm256_extractf128_pd(self, 1));
            __m256d hi_hi = _mm256_insertf128_pd(self, _mm256_castpd256_pd128(hi), 0);

            __m256d low = _mm256_castpd128_pd256(_mm256_castpd256_pd128(self));
            __m256d low_low = _mm256_insertf128_pd(self, _mm256_castpd256_pd128(low), 1);

            // normalize mask
            batch_constant<uint64_t, A, (V0 % 2) * -1, (V1 % 2) * -1, (V2 % 2) * -1, (V3 % 2) * -1> half_mask;

            // permute within each lane
            __m256d r0 = _mm256_permutevar_pd(low_low, half_mask.as_batch());
            __m256d r1 = _mm256_permutevar_pd(hi_hi, half_mask.as_batch());

            // mask to choose the right lane
            batch_bool_constant<uint64_t, A, (V0 >= 2), (V1 >= 2), (V2 >= 2), (V3 >= 2)> blend_mask;

            // blend the two permutes
            constexpr auto mask = blend_mask.mask();
            return _mm256_blend_pd(r0, r1, mask);
        }
        template <class A,
                  typename T,
                  uint32_t V0,
                  uint32_t V1,
                  uint32_t V2,
                  uint32_t V3,
                  uint32_t V4,
                  uint32_t V5,
                  uint32_t V6,
                  uint32_t V7,
                  detail::enable_sized_integral_t<T, 4> = 0>
        XSIMD_INLINE batch<T, A> swizzle(batch<T, A> const& self,
                                         batch_constant<uint32_t, A,
                                                        V0,
                                                        V1,
                                                        V2,
                                                        V3,
                                                        V4,
                                                        V5,
                                                        V6,
                                                        V7> const& mask,
                                         requires_arch<avx>) noexcept
        {
            return bitwise_cast<T>(
                swizzle(bitwise_cast<float>(self), mask));
        }

        template <class A,
                  typename T,
                  uint64_t V0,
                  uint64_t V1,
                  uint64_t V2,
                  uint64_t V3,
                  detail::enable_sized_integral_t<T, 8> = 0>
        XSIMD_INLINE batch<T, A>
        swizzle(batch<T, A> const& self,
                batch_constant<uint64_t, A, V0, V1, V2, V3> const& mask,
                requires_arch<avx>) noexcept
        {
            return bitwise_cast<T>(
                swizzle(bitwise_cast<double>(self), mask));
        }

        // trunc
        template <class A>
        XSIMD_INLINE batch<float, A> trunc(batch<float, A> const& self, requires_arch<avx>) noexcept
        {
            return _mm256_round_ps(self, _MM_FROUND_TO_ZERO);
        }
        template <class A>
        XSIMD_INLINE batch<double, A> trunc(batch<double, A> const& self, requires_arch<avx>) noexcept
        {
            return _mm256_round_pd(self, _MM_FROUND_TO_ZERO);
        }

        // zip_hi
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> zip_hi(batch<T, A> const& self, batch<T, A> const& other, requires_arch<avx>) noexcept
        {
            XSIMD_IF_CONSTEXPR(sizeof(T) == 1 || sizeof(T) == 2)
            {
                // extract high word
                __m128i self_hi = _mm256_extractf128_si256(self, 1);
                __m128i other_hi = _mm256_extractf128_si256(other, 1);

                // interleave
                __m128i res_lo, res_hi;
                XSIMD_IF_CONSTEXPR(sizeof(T) == 1)
                {
                    res_lo = _mm_unpacklo_epi8(self_hi, other_hi);
                    res_hi = _mm_unpackhi_epi8(self_hi, other_hi);
                }
                else
                {
                    res_lo = _mm_unpacklo_epi16(self_hi, other_hi);
                    res_hi = _mm_unpackhi_epi16(self_hi, other_hi);
                }

                // fuse
                return _mm256_castps_si256(
                    _mm256_insertf128_ps(
                        _mm256_castsi256_ps(_mm256_castsi128_si256(res_lo)),
                        _mm_castsi128_ps(res_hi),
                        1));
            }
            else XSIMD_IF_CONSTEXPR(sizeof(T) == 4)
            {
                auto lo = _mm256_unpacklo_ps(_mm256_castsi256_ps(self), _mm256_castsi256_ps(other));
                auto hi = _mm256_unpackhi_ps(_mm256_castsi256_ps(self), _mm256_castsi256_ps(other));
                return _mm256_castps_si256(_mm256_permute2f128_ps(lo, hi, 0x31));
            }
            else XSIMD_IF_CONSTEXPR(sizeof(T) == 8)
            {
                auto lo = _mm256_unpacklo_pd(_mm256_castsi256_pd(self), _mm256_castsi256_pd(other));
                auto hi = _mm256_unpackhi_pd(_mm256_castsi256_pd(self), _mm256_castsi256_pd(other));
                return _mm256_castpd_si256(_mm256_permute2f128_pd(lo, hi, 0x31));
            }
            else
            {
                assert(false && "unsupported arch/op combination");
                return {};
            }
        }
        template <class A>
        XSIMD_INLINE batch<float, A> zip_hi(batch<float, A> const& self, batch<float, A> const& other, requires_arch<avx>) noexcept
        {
            auto lo = _mm256_unpacklo_ps(self, other);
            auto hi = _mm256_unpackhi_ps(self, other);
            return _mm256_permute2f128_ps(lo, hi, 0x31);
        }
        template <class A>
        XSIMD_INLINE batch<double, A> zip_hi(batch<double, A> const& self, batch<double, A> const& other, requires_arch<avx>) noexcept
        {
            auto lo = _mm256_unpacklo_pd(self, other);
            auto hi = _mm256_unpackhi_pd(self, other);
            return _mm256_permute2f128_pd(lo, hi, 0x31);
        }

        // zip_lo
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> zip_lo(batch<T, A> const& self, batch<T, A> const& other, requires_arch<avx>) noexcept
        {
            XSIMD_IF_CONSTEXPR(sizeof(T) == 1 || sizeof(T) == 2)
            {
                // extract low word
                __m128i self_lo = _mm256_extractf128_si256(self, 0);
                __m128i other_lo = _mm256_extractf128_si256(other, 0);

                // interleave
                __m128i res_lo, res_hi;
                XSIMD_IF_CONSTEXPR(sizeof(T) == 1)
                {
                    res_lo = _mm_unpacklo_epi8(self_lo, other_lo);
                    res_hi = _mm_unpackhi_epi8(self_lo, other_lo);
                }
                else
                {
                    res_lo = _mm_unpacklo_epi16(self_lo, other_lo);
                    res_hi = _mm_unpackhi_epi16(self_lo, other_lo);
                }

                // fuse
                return _mm256_castps_si256(
                    _mm256_insertf128_ps(
                        _mm256_castsi256_ps(_mm256_castsi128_si256(res_lo)),
                        _mm_castsi128_ps(res_hi),
                        1));
            }
            else XSIMD_IF_CONSTEXPR(sizeof(T) == 4)
            {
                auto lo = _mm256_unpacklo_ps(_mm256_castsi256_ps(self), _mm256_castsi256_ps(other));
                auto hi = _mm256_unpackhi_ps(_mm256_castsi256_ps(self), _mm256_castsi256_ps(other));
                return _mm256_castps_si256(_mm256_insertf128_ps(lo, _mm256_castps256_ps128(hi), 1));
            }
            else XSIMD_IF_CONSTEXPR(sizeof(T) == 8)
            {
                auto lo = _mm256_unpacklo_pd(_mm256_castsi256_pd(self), _mm256_castsi256_pd(other));
                auto hi = _mm256_unpackhi_pd(_mm256_castsi256_pd(self), _mm256_castsi256_pd(other));
                return _mm256_castpd_si256(_mm256_insertf128_pd(lo, _mm256_castpd256_pd128(hi), 1));
            }
            else
            {
                assert(false && "unsupported arch/op combination");
                return {};
            }
        }

        template <class A>
        XSIMD_INLINE batch<float, A> zip_lo(batch<float, A> const& self, batch<float, A> const& other, requires_arch<avx>) noexcept
        {
            auto lo = _mm256_unpacklo_ps(self, other);
            auto hi = _mm256_unpackhi_ps(self, other);
            return _mm256_insertf128_ps(lo, _mm256_castps256_ps128(hi), 1);
        }
        template <class A>
        XSIMD_INLINE batch<double, A> zip_lo(batch<double, A> const& self, batch<double, A> const& other, requires_arch<avx>) noexcept
        {
            auto lo = _mm256_unpacklo_pd(self, other);
            auto hi = _mm256_unpackhi_pd(self, other);
            return _mm256_insertf128_pd(lo, _mm256_castpd256_pd128(hi), 1);
        }
    }
}

#endif
