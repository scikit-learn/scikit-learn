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

#ifndef XSIMD_API_HPP
#define XSIMD_API_HPP

#include <complex>
#include <cstddef>
#include <limits>
#include <ostream>

#include "../arch/xsimd_isa.hpp"
#include "../types/xsimd_batch.hpp"
#include "../types/xsimd_traits.hpp"

namespace xsimd
{
    /**
     * high level free functions
     *
     * @defgroup batch_arithmetic Arithmetic operators
     * @defgroup batch_constant Constant batches
     * @defgroup batch_data_transfer Memory operators
     * @defgroup batch_math Basic math operators
     * @defgroup batch_math_extra Extra math operators
     * @defgroup batch_fp Floating point manipulation
     * @defgroup batch_rounding Rounding operators
     * @defgroup batch_conversion Conversion operators
     * @defgroup batch_complex_op Complex operators
     * @defgroup batch_logical Logical operators
     * @defgroup batch_bitwise Bitwise operators
     * @defgroup batch_reducers Reducers
     * @defgroup batch_miscellaneous Miscellaneous
     * @defgroup batch_trigo Trigonometry
     *
     * @defgroup batch_bool_logical Boolean logical operators
     * @defgroup batch_bool_reducers Boolean reducers
     */

    /**
     * @ingroup batch_math
     *
     * Computes the absolute values of each scalar in the batch \c x.
     * @param x batch of integer or floating point values.
     * @return the absolute values of \c x.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> abs(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::abs<A>(x, A {});
    }

    /**
     * @ingroup batch_complex
     *
     * Computes the absolute values of each complex in the batch \c z.
     * @param z batch of complex values.
     * @return the absolute values of \c z.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> abs(batch<std::complex<T>, A> const& z) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::abs<A>(z, A {});
    }

    /**
     * @ingroup batch_arithmetic
     *
     * Computes the sum of the batches \c x and \c y.
     * @param x batch or scalar involved in the addition.
     * @param y batch or scalar involved in the addition.
     * @return the sum of \c x and \c y
     */
    template <class T, class A>
    XSIMD_INLINE auto add(batch<T, A> const& x, batch<T, A> const& y) noexcept -> decltype(x + y)
    {
        detail::static_check_supported_config<T, A>();
        return x + y;
    }

    /**
     * @ingroup batch_trigo
     *
     * Computes the arc cosine of the batch \c x.
     * @param x batch of floating point values.
     * @return the arc cosine of \c x.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> acos(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::acos<A>(x, A {});
    }

    /**
     * @ingroup batch_trigo
     *
     * Computes the inverse hyperbolic cosine of the batch \c x.
     * @param x batch of floating point values.
     * @return the inverse hyperbolic cosine of \c x.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> acosh(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::acosh<A>(x, A {});
    }

    /**
     * @ingroup batch_complex
     *
     * Computes the argument of the batch \c z.
     * @param z batch of complex or real values.
     * @return the argument of \c z.
     */
    template <class T, class A>
    XSIMD_INLINE real_batch_type_t<batch<T, A>> arg(batch<T, A> const& z) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::arg<A>(z, A {});
    }

    /**
     * @ingroup batch_trigo
     *
     * Computes the arc sine of the batch \c x.
     * @param x batch of floating point values.
     * @return the arc sine of \c x.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> asin(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::asin<A>(x, A {});
    }

    /**
     * @ingroup batch_trigo
     *
     * Computes the inverse hyperbolic sine of the batch \c x.
     * @param x batch of floating point values.
     * @return the inverse hyperbolic sine of \c x.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> asinh(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::asinh<A>(x, A {});
    }

    /**
     * @ingroup batch_trigo
     *
     * Computes the arc tangent of the batch \c x.
     * @param x batch of floating point values.
     * @return the arc tangent of \c x.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> atan(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::atan<A>(x, A {});
    }

    /**
     * @ingroup batch_trigo
     *
     * Computes the arc tangent of the batch \c x/y, using the signs of the
     * arguments to determine the correct quadrant.
     * @param x batch of floating point values.
     * @param y batch of floating point values.
     * @return the arc tangent of \c x/y.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> atan2(batch<T, A> const& x, batch<T, A> const& y) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::atan2<A>(x, y, A {});
    }

    /**
     * @ingroup batch_trigo
     *
     * Computes the inverse hyperbolic tangent of the batch \c x.
     * @param x batch of floating point values.
     * @return the inverse hyperbolic tangent of \c x.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> atanh(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::atanh<A>(x, A {});
    }

    /**
     * @ingroup batch_math
     *
     * Computes the average of batches \c x and \c y
     * @param x batch of T
     * @param y batch of T
     * @return the average of elements between \c x and \c y.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> avg(batch<T, A> const& x, batch<T, A> const& y) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::avg<A>(x, y, A {});
    }

    /**
     * @ingroup batch_math
     *
     * Computes the rounded average of batches \c x and \c y
     * @param x batch of T
     * @param y batch of T
     * @return the rounded average of elements between \c x and \c y.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> avgr(batch<T, A> const& x, batch<T, A> const& y) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::avgr<A>(x, y, A {});
    }

    /**
     * @ingroup batch_conversion
     *
     * Perform a static_cast from \c T_in to \c T_out on \c \c x.
     * @param x batch_bool of \c T_in
     * @return \c x cast to \c T_out
     */
    template <class T_out, class T_in, class A>
    XSIMD_INLINE batch_bool<T_out, A> batch_bool_cast(batch_bool<T_in, A> const& x) noexcept
    {
        detail::static_check_supported_config<T_out, A>();
        detail::static_check_supported_config<T_in, A>();
        static_assert(batch_bool<T_out, A>::size == batch_bool<T_in, A>::size, "Casting between incompatibles batch_bool types.");
        return kernel::batch_bool_cast<A>(x, batch_bool<T_out, A> {}, A {});
    }

    /**
     * @ingroup batch_conversion
     *
     * Perform a static_cast from \c T_in to \c T_out on \c \c x.
     * @param x batch of \c T_in
     * @return \c x cast to \c T_out
     */
    template <class T_out, class T_in, class A>
    XSIMD_INLINE batch<T_out, A> batch_cast(batch<T_in, A> const& x) noexcept
    {
        detail::static_check_supported_config<T_out, A>();
        detail::static_check_supported_config<T_in, A>();
        return kernel::batch_cast<A>(x, batch<T_out, A> {}, A {});
    }

    /**
     * @ingroup batch_miscellaneous
     *
     * Computes the bit of sign of \c x
     * @param x batch of scalar
     * @return bit of sign of \c x
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> bitofsign(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::bitofsign<A>(x, A {});
    }

    /**
     * @ingroup batch_bitwise
     *
     * Computes the bitwise and of the batches \c x and \c y.
     * @param x batch involved in the operation.
     * @param y batch involved in the operation.
     * @return the result of the bitwise and.
     */
    template <class T, class A>
    XSIMD_INLINE auto bitwise_and(batch<T, A> const& x, batch<T, A> const& y) noexcept -> decltype(x & y)
    {
        detail::static_check_supported_config<T, A>();
        return x & y;
    }

    /**
     * @ingroup batch_bitwise
     *
     * Computes the bitwise and of the batches \c x and \c y.
     * @param x batch involved in the operation.
     * @param y batch involved in the operation.
     * @return the result of the bitwise and.
     */
    template <class T, class A>
    XSIMD_INLINE auto bitwise_and(batch_bool<T, A> const& x, batch_bool<T, A> const& y) noexcept -> decltype(x & y)
    {
        detail::static_check_supported_config<T, A>();
        return x & y;
    }

    /**
     * @ingroup batch_bitwise
     *
     * Computes the bitwise and not of batches \c x and \c y.
     * @param x batch involved in the operation.
     * @param y batch involved in the operation.
     * @return the result of the bitwise and not.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> bitwise_andnot(batch<T, A> const& x, batch<T, A> const& y) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::bitwise_andnot<A>(x, y, A {});
    }

    /**
     * @ingroup batch_bool_logical
     *
     * Computes the bitwise and not of batches \c x and \c y.
     * @param x batch involved in the operation.
     * @param y batch involved in the operation.
     * @return the result of the bitwise and not.
     */
    template <class T, class A>
    XSIMD_INLINE batch_bool<T, A> bitwise_andnot(batch_bool<T, A> const& x, batch_bool<T, A> const& y) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::bitwise_andnot<A>(x, y, A {});
    }

    /**
     * @ingroup batch_conversion
     *
     * Perform a reinterpret_cast from \c T_in to \c T_out on \c x.
     * @param x batch of \c T_in
     * @return \c x reinterpreted as \c T_out
     */
    template <class T_out, class T_in, class A>
    XSIMD_INLINE batch<T_out, A> bitwise_cast(batch<T_in, A> const& x) noexcept
    {
        detail::static_check_supported_config<T_in, A>();
        detail::static_check_supported_config<T_out, A>();
        return kernel::bitwise_cast<A>(x, batch<T_out, A> {}, A {});
    }

    /**
     * @ingroup batch_bitwise
     *
     * Perform a bitwise shift to the left
     * @param x batch of \c T_in
     * @param shift scalar amount to shift
     * @return shifted \c x.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> bitwise_lshift(batch<T, A> const& x, int shift) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::bitwise_lshift<A>(x, shift, A {});
    }
    template <class T, class A>
    XSIMD_INLINE batch<T, A> bitwise_lshift(batch<T, A> const& x, batch<T, A> const& shift) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::bitwise_lshift<A>(x, shift, A {});
    }

    /**
     * @ingroup batch_bitwise
     *
     * Computes the bitwise not of batch \c x.
     * @param x batch involved in the operation.
     * @return the result of the bitwise not.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> bitwise_not(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::bitwise_not<A>(x, A {});
    }

    /**
     * @ingroup batch_bitwise
     *
     * Computes the bitwise not of batch \c x.
     * @param x batch involved in the operation.
     * @return the result of the bitwise not.
     */
    template <class T, class A>
    XSIMD_INLINE batch_bool<T, A> bitwise_not(batch_bool<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::bitwise_not<A>(x, A {});
    }

    /**
     * @ingroup batch_bitwise
     *
     * Computes the bitwise or of the batches \c x and \c y.
     * @param x scalar or batch of scalars
     * @param y scalar or batch of scalars
     * @return the result of the bitwise or.
     */
    template <class T, class A>
    XSIMD_INLINE auto bitwise_or(batch<T, A> const& x, batch<T, A> const& y) noexcept -> decltype(x | y)
    {
        detail::static_check_supported_config<T, A>();
        return x | y;
    }

    /**
     * @ingroup batch_bitwise
     *
     * Computes the bitwise or of the batches \c x and \c y.
     * @param x scalar or batch of scalars
     * @param y scalar or batch of scalars
     * @return the result of the bitwise or.
     */
    template <class T, class A>
    XSIMD_INLINE auto bitwise_or(batch_bool<T, A> const& x, batch_bool<T, A> const& y) noexcept -> decltype(x | y)
    {
        detail::static_check_supported_config<T, A>();
        return x | y;
    }

    /**
     * @ingroup batch_bitwise
     *
     * Perform a bitwise shift to the right
     * @param x batch of \c T_in
     * @param shift scalar amount to shift
     * @return shifted \c x.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> bitwise_rshift(batch<T, A> const& x, int shift) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::bitwise_rshift<A>(x, shift, A {});
    }
    template <class T, class A>
    XSIMD_INLINE batch<T, A> bitwise_rshift(batch<T, A> const& x, batch<T, A> const& shift) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::bitwise_rshift<A>(x, shift, A {});
    }

    /**
     * @ingroup batch_bitwise
     *
     * Computes the bitwise xor of the batches \c x and \c y.
     * @param x scalar or batch of scalars
     * @param y scalar or batch of scalars
     * @return the result of the bitwise xor.
     */
    template <class T, class A>
    XSIMD_INLINE auto bitwise_xor(batch<T, A> const& x, batch<T, A> const& y) noexcept -> decltype(x ^ y)
    {
        detail::static_check_supported_config<T, A>();
        return x ^ y;
    }

    /**
     * @ingroup batch_bitwise
     *
     * Computes the bitwise xor of the batches \c x and \c y.
     * @param x scalar or batch of scalars
     * @param y scalar or batch of scalars
     * @return the result of the bitwise xor.
     */
    template <class T, class A>
    XSIMD_INLINE auto bitwise_xor(batch_bool<T, A> const& x, batch_bool<T, A> const& y) noexcept -> decltype(x ^ y)
    {
        detail::static_check_supported_config<T, A>();
        return x ^ y;
    }

    /**
     * @ingroup batch_data_transfer
     *
     * Creates a batch from the single value \c v.
     * @param v the value used to initialize the batch
     * @return a new batch instance
     */
    template <class T, class A = default_arch>
    XSIMD_INLINE batch<T, A> broadcast(T v) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return batch<T, A>::broadcast(v);
    }

    /**
     * @ingroup batch_data_transfer
     *
     * Creates a batch from the single value \c v and
     * the specified batch value type \c To.
     * @param v the value used to initialize the batch
     * @return a new batch instance
     */
    template <class To, class A = default_arch, class From>
    XSIMD_INLINE simd_return_type<From, To, A> broadcast_as(From v) noexcept
    {
        detail::static_check_supported_config<From, A>();
        using batch_value_type = typename simd_return_type<From, To, A>::value_type;
        using value_type = typename std::conditional<std::is_same<From, bool>::value,
                                                     bool,
                                                     batch_value_type>::type;
        return simd_return_type<From, To, A>(value_type(v));
    }

    /**
     * @ingroup batch_math
     *
     * Computes the cubic root of the batch \c x.
     * @param x batch of floating point values.
     * @return the cubic root of \c x.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> cbrt(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::cbrt<A>(x, A {});
    }

    /**
     * @ingroup batch_rounding
     *
     * Computes the batch of smallest integer values not less than
     * scalars in \c x.
     * @param x batch of floating point values.
     * @return the batch of smallest integer values not less than \c x.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> ceil(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::ceil<A>(x, A {});
    }

    /**
     * @ingroup batch_math
     *
     * Clips the values of the batch \c x between those of the batches \c lo and \c hi.
     * @param x batch of scalar values.
     * @param lo batch of scalar values.
     * @param hi batch of scalar values.
     * @return the result of the clipping.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> clip(batch<T, A> const& x, batch<T, A> const& lo, batch<T, A> const& hi) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::clip(x, lo, hi, A {});
    }

    /**
     * @ingroup batch_data_transfer
     *
     * Pick elements from \c x selected by \c mask, and append them to the
     * resulting vector, zeroing the remaining slots
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> compress(batch<T, A> const& x, batch_bool<T, A> const& mask) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::compress<A>(x, mask, A {});
    }

    /**
     * @ingroup batch_complex
     *
     * Computes the conjugate of the batch \c z.
     * @param z batch of complex values.
     * @return the argument of \c z.
     */
    template <class A, class T>
    XSIMD_INLINE complex_batch_type_t<batch<T, A>> conj(batch<T, A> const& z) noexcept
    {
        return kernel::conj(z, A {});
    }

    /**
     * @ingroup batch_miscellaneous
     *
     * Computes a value whose  absolute  value  matches
     *        that of \c x, but whose sign bit matches that of \c y.
     * @param x batch of scalars
     * @param y batch of scalars
     * @return batch whose absolute  value  matches that of \c x, but whose sign bit
     * matches that of \c y.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> copysign(batch<T, A> const& x, batch<T, A> const& y) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::copysign<A>(x, y, A {});
    }

    /**
     * @ingroup batch_trigo
     *
     * Computes the cosine of the batch \c x.
     * @param x batch of floating point values.
     * @return the cosine of \c x.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> cos(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::cos<A>(x, A {});
    }

    /**
     * @ingroup batch_trigo
     *
     * computes the hyperbolic cosine of the batch \c x.
     * @param x batch of floating point values.
     * @return the hyperbolic cosine of \c x.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> cosh(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::cosh<A>(x, A {});
    }

    /**
     * @ingroup batch_arithmetic
     *
     * Subtract 1 to batch \c x.
     * @param x batch involved in the decrement.
     * @return the subtraction of \c x and 1.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> decr(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::decr<A>(x, A {});
    }

    /**
     * @ingroup batch_arithmetic
     *
     * Subtract 1 to batch \c x for each element where \c mask is true.
     * @param x batch involved in the increment.
     * @param mask whether to perform the increment or not. Can be a \c
     *             batch_bool or a \c batch_bool_constant.
     * @return the subtraction of \c x and 1 when \c mask is true.
     */
    template <class T, class A, class Mask>
    XSIMD_INLINE batch<T, A> decr_if(batch<T, A> const& x, Mask const& mask) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::decr_if<A>(x, mask, A {});
    }

    /**
     * @ingroup batch_arithmetic
     *
     * Computes the division of the batch \c x by the batch \c y.
     * @param x scalar or batch of scalars
     * @param y scalar or batch of scalars
     * @return the result of the division.
     */
    template <class T, class A>
    XSIMD_INLINE auto div(batch<T, A> const& x, batch<T, A> const& y) noexcept -> decltype(x / y)
    {
        detail::static_check_supported_config<T, A>();
        return x / y;
    }

    /**
     * @ingroup batch_logical
     *
     * Element-wise equality comparison of batches \c x and \c y.
     * @param x batch of scalars
     * @param y batch of scalars
     * @return a boolean batch.
     */
    template <class T, class A>
    XSIMD_INLINE auto eq(batch<T, A> const& x, batch<T, A> const& y) noexcept -> decltype(x == y)
    {
        detail::static_check_supported_config<T, A>();
        return x == y;
    }

    /**
     * @ingroup batch_logical
     *
     * Element-wise equality comparison of batches of boolean values \c x and \c y.
     * @param x batch of booleans involved in the comparison.
     * @param y batch of booleans involved in the comparison.
     * @return a boolean batch.
     */
    template <class T, class A>
    XSIMD_INLINE auto eq(batch_bool<T, A> const& x, batch_bool<T, A> const& y) noexcept -> decltype(x == y)
    {
        detail::static_check_supported_config<T, A>();
        return x == y;
    }

    /**
     * @ingroup batch_math
     *
     * Computes the natural exponential of the batch \c x.
     * @param x batch of floating point values.
     * @return the natural exponential of \c x.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> exp(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::exp<A>(x, A {});
    }

    /**
     * @ingroup batch_math
     *
     * Computes the base 10 exponential of the batch \c x.
     * @param x batch of floating point values.
     * @return the base 10 exponential of \c x.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> exp10(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::exp10<A>(x, A {});
    }

    /**
     * @ingroup batch_math
     *
     * Computes the base 2 exponential of the batch \c x.
     * @param x batch of floating point values.
     * @return the base 2 exponential of \c x.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> exp2(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::exp2<A>(x, A {});
    }

    /**
     * @ingroup batch_data_transfer
     *
     * Load contiguous elements from \c x and place them in slots selected by \c
     * mask, zeroing the other slots
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> expand(batch<T, A> const& x, batch_bool<T, A> const& mask) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::expand<A>(x, mask, A {});
    }

    /**
     * @ingroup batch_math
     *
     * Computes the natural exponential of the batch \c x, minus one.
     * @param x batch of floating point values.
     * @return the natural exponential of \c x, minus one.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> expm1(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::expm1<A>(x, A {});
    }

    /**
     * @ingroup batch_math_extra
     *
     * Computes the error function of the batch \c x.
     * @param x batch of floating point values.
     * @return the error function of \c x.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> erf(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::erf<A>(x, A {});
    }

    /**
     * @ingroup batch_math_extra
     *
     * Computes the complementary error function of the batch \c x.
     * @param x batch of floating point values.
     * @return the error function of \c x.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> erfc(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::erfc<A>(x, A {});
    }

    /**
     * Extract vector from pair of vectors
     * extracts the lowest vector elements from the second source \c x
     * and the highest vector elements from the first source \c y
     * Concatenates the results into th Return value.
     * @param x batch of integer or floating point values.
     * @param y batch of integer or floating point values.
     * @param i integer specifying the lowest vector element to extract from the first source register
     * @return.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> extract_pair(batch<T, A> const& x, batch<T, A> const& y, std::size_t i) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::extract_pair<A>(x, y, i, A {});
    }

    /**
     * @ingroup batch_math
     *
     * Computes the absolute values of each scalar in the batch \c x.
     * @param x batch floating point values.
     * @return the absolute values of \c x.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> fabs(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::abs<A>(x, A {});
    }

    /**
     * @ingroup batch_math
     *
     * Computes the positive difference between \c x and \c y, that is,
     * <tt>max(0, x-y)</tt>.
     * @param x batch of floating point values.
     * @param y batch of floating point values.
     * @return the positive difference.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> fdim(batch<T, A> const& x, batch<T, A> const& y) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::fdim<A>(x, y, A {});
    }

    /**
     * @ingroup batch_rounding
     *
     * Computes the batch of largest integer values not greater than
     * scalars in \c x.
     * @param x batch of floating point values.
     * @return the batch of largest integer values not greater than \c x.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> floor(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::floor<A>(x, A {});
    }

    /**
     * @ingroup batch_arithmetic
     *
     * Computes <tt>(x*y) + z</tt> in a single instruction when possible.
     * @param x a batch of integer or floating point values.
     * @param y a batch of integer or floating point values.
     * @param z a batch of integer or floating point values.
     * @return the result of the fused multiply-add operation.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> fma(batch<T, A> const& x, batch<T, A> const& y, batch<T, A> const& z) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::fma<A>(x, y, z, A {});
    }

    /**
     * @ingroup batch_math
     *
     * Computes the larger values of the batches \c x and \c y.
     * @param x a batch of integer or floating point values.
     * @param y a batch of integer or floating point values.
     * @return a batch of the larger values.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> fmax(batch<T, A> const& x, batch<T, A> const& y) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::max<A>(x, y, A {});
    }

    /**
     * @ingroup batch_math
     *
     * Computes the smaller values of the batches \c x and \c y.
     * @param x a batch of integer or floating point values.
     * @param y a batch of integer or floating point values.
     * @return a batch of the smaller values.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> fmin(batch<T, A> const& x, batch<T, A> const& y) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::min<A>(x, y, A {});
    }

    /**
     * @ingroup batch_math
     *
     * Computes the modulo of the batch \c x by the batch \c y.
     * @param x batch involved in the modulo.
     * @param y batch involved in the modulo.
     * @return the result of the modulo.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> fmod(batch<T, A> const& x, batch<T, A> const& y) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::fmod<A>(x, y, A {});
    }

    /**
     * @ingroup batch_arithmetic
     *
     * Computes <tt>(x*y) - z</tt> in a single instruction when possible.
     * @param x a batch of integer or floating point values.
     * @param y a batch of integer or floating point values.
     * @param z a batch of integer or floating point values.
     * @return the result of the fused multiply-sub operation.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> fms(batch<T, A> const& x, batch<T, A> const& y, batch<T, A> const& z) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::fms<A>(x, y, z, A {});
    }

    /**
     * @ingroup batch_arithmetic
     *
     * Computes <tt>-(x*y) + z</tt> in a single instruction when possible.
     * @param x a batch of integer or floating point values.
     * @param y a batch of integer or floating point values.
     * @param z a batch of integer or floating point values.
     * @return the result of the fused negated multiply-add operation.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> fnma(batch<T, A> const& x, batch<T, A> const& y, batch<T, A> const& z) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::fnma<A>(x, y, z, A {});
    }

    /**
     * @ingroup batch_arithmetic
     *
     * Computes <tt>-(x*y) - z</tt> in a single instruction when possible.
     * @param x a batch of integer or floating point values.
     * @param y a batch of integer or floating point values.
     * @param z a batch of integer or floating point values.
     * @return the result of the fused negated multiply-sub operation.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> fnms(batch<T, A> const& x, batch<T, A> const& y, batch<T, A> const& z) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::fnms<A>(x, y, z, A {});
    }

    /**
     * @ingroup batch_fp
     *
     * Split split the number x into a normalized fraction and an exponent which is stored in exp
     * @param x a batch of integer or floating point values.
     * @param y a batch of integer or floating point values.
     * @return the normalized fraction of x
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> frexp(const batch<T, A>& x, batch<as_integer_t<T>, A>& y) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::frexp<A>(x, y, A {});
    }

    /**
     * @ingroup batch_logical
     *
     * Element-wise greater or equal comparison of batches \c x and \c y.
     * @tparam X the actual type of batch.
     * @param x batch involved in the comparison.
     * @param y batch involved in the comparison.
     * @return a boolean batch.
     */
    template <class T, class A>
    XSIMD_INLINE batch_bool<T, A> ge(batch<T, A> const& x, batch<T, A> const& y) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return x >= y;
    }

    /**
     * @ingroup batch_logical
     *
     * Element-wise greater than comparison of batches \c x and \c y.
     * @tparam X the actual type of batch.
     * @param x batch involved in the comparison.
     * @param y batch involved in the comparison.
     * @return a boolean batch.
     */
    template <class T, class A>
    XSIMD_INLINE batch_bool<T, A> gt(batch<T, A> const& x, batch<T, A> const& y) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return x > y;
    }

    /**
     * @ingroup batch_reducers
     *
     * Parallel horizontal addition: adds the scalars of each batch
     * in the array pointed by \c row and store them in a returned
     * batch.
     * @param row an array of \c N batches
     * @return the result of the reduction.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> haddp(batch<T, A> const* row) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::haddp<A>(row, A {});
    }

    /**
     * @ingroup batch_math
     *
     * Computes the square root of the sum of the squares of the batches
     * \c x, and \c y.
     * @param x batch of floating point values.
     * @param y batch of floating point values.
     * @return the square root of the sum of the squares of \c x and \c y.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> hypot(batch<T, A> const& x, batch<T, A> const& y) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::hypot<A>(x, y, A {});
    }

    /**
     * @ingroup batch_complex
     *
     * Computes the imaginary part of the batch \c x.
     * @param x batch of complex or real values.
     * @return the argument of \c x.
     */
    template <class T, class A>
    XSIMD_INLINE real_batch_type_t<batch<T, A>> imag(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::imag<A>(x, A {});
    }

    /**
     * @ingroup batch_arithmetic
     *
     * Add 1 to batch \c x.
     * @param x batch involved in the increment.
     * @return the sum of \c x and 1.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> incr(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::incr<A>(x, A {});
    }

    /**
     * @ingroup batch_arithmetic
     *
     * Add 1 to batch \c x for each element where \c mask is true.
     * @param x batch involved in the increment.
     * @param mask whether to perform the increment or not. Can be a \c
     *             batch_bool or a \c batch_bool_constant.
     * @return the sum of \c x and 1 when \c mask is true.
     */
    template <class T, class A, class Mask>
    XSIMD_INLINE batch<T, A> incr_if(batch<T, A> const& x, Mask const& mask) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::incr_if<A>(x, mask, A {});
    }

    /**
     * @ingroup batch_constant
     *
     * Return a batch of scalars representing positive infinity
     * @return a batch of positive infinity
     */
    template <class B>
    XSIMD_INLINE B infinity()
    {
        using T = typename B::value_type;
        using A = typename B::arch_type;
        detail::static_check_supported_config<T, A>();
        return B(std::numeric_limits<T>::infinity());
    }

    /**
     * @ingroup batch_data_transfer
     *
     * Create a new batch equivalent to \c x but with element \c val set at position \c pos
     * @param x batch
     * @param val value to set
     * @param pos index of the updated slot
     * @return copy of \c x with position \c pos set to \c val
     */
    template <class T, class A, size_t I>
    XSIMD_INLINE batch<T, A> insert(batch<T, A> const& x, T val, index<I> pos) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::insert<A>(x, val, pos, A {});
    }

    /**
     * @ingroup batch_logical
     *
     * Determines if the scalars in the given batch \c x represent an even integer value
     * @param x batch of floating point values.
     * @return a batch of booleans.
     */
    template <class T, class A>
    XSIMD_INLINE batch_bool<T, A> is_even(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::is_even<A>(x, A {});
    }

    /**
     * @ingroup batch_logical
     *
     * Determines if the floating-point scalars in the given batch \c x represent integer value
     * @param x batch of floating point values.
     * @return a batch of booleans.
     */
    template <class T, class A>
    XSIMD_INLINE batch_bool<T, A> is_flint(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::is_flint<A>(x, A {});
    }

    /**
     * @ingroup batch_logical
     *
     * Determines if the scalars in the given batch \c x represent an odd integer value
     * @param x batch of floating point values.
     * @return a batch of booleans.
     */
    template <class T, class A>
    XSIMD_INLINE batch_bool<T, A> is_odd(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::is_odd<A>(x, A {});
    }

    /**
     * @ingroup batch_logical
     *
     * Determines if the scalars in the given batch \c x are inf values.
     * @param x batch of floating point values.
     * @return a batch of booleans.
     */
    template <class T, class A>
    XSIMD_INLINE typename batch<T, A>::batch_bool_type isinf(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::isinf<A>(x, A {});
    }

    /**
     * @ingroup batch_logical
     *
     * Determines if the scalars in the given batch \c x are finite values.
     * @param x batch of floating point values.
     * @return a batch of booleans.
     */
    template <class T, class A>
    XSIMD_INLINE typename batch<T, A>::batch_bool_type isfinite(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::isfinite<A>(x, A {});
    }

    /**
     * @ingroup batch_logical
     *
     * Determines if the scalars in the given batch \c x are NaN values.
     * @param x batch of floating point values.
     * @return a batch of booleans.
     */
    template <class T, class A>
    XSIMD_INLINE typename batch<T, A>::batch_bool_type isnan(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::isnan<A>(x, A {});
    }

    /**
     * @ingroup batch_math_extra
     *
     * Computes the multiplication of the floating point number \c x by 2 raised to the power \c y.
     * @param x batch of floating point values.
     * @param y batch of integer values.
     * @return a batch of floating point values.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> ldexp(const batch<T, A>& x, const batch<as_integer_t<T>, A>& y) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::ldexp<A>(x, y, A {});
    }

    /**
     * @ingroup batch_logical
     *
     * Element-wise lesser or equal to comparison of batches \c x and \c y.
     * @param x batch involved in the comparison.
     * @param y batch involved in the comparison.
     * @return a boolean batch.
     */
    template <class T, class A>
    XSIMD_INLINE batch_bool<T, A> le(batch<T, A> const& x, batch<T, A> const& y) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return x <= y;
    }

    /**
     * @ingroup batch_math_extra
     *
     * Computes the natural logarithm of the gamma function of the batch \c x.
     * @param x batch of floating point values.
     * @return the natural logarithm of the gamma function of \c x.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> lgamma(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::lgamma<A>(x, A {});
    }

    /**
     * @ingroup batch_data_transfer
     *
     * Creates a batch from the buffer \c ptr and the specifed
     * batch value type \c To. The memory needs to be aligned.
     * @param ptr the memory buffer to read
     * @return a new batch instance
     */
    template <class To, class A = default_arch, class From>
    XSIMD_INLINE simd_return_type<From, To, A> load_as(From const* ptr, aligned_mode) noexcept
    {
        using batch_value_type = typename simd_return_type<From, To, A>::value_type;
        detail::static_check_supported_config<From, A>();
        detail::static_check_supported_config<To, A>();
        return kernel::load_aligned<A>(ptr, kernel::convert<batch_value_type> {}, A {});
    }

    template <class To, class A = default_arch>
    XSIMD_INLINE simd_return_type<bool, To, A> load_as(bool const* ptr, aligned_mode) noexcept
    {
        detail::static_check_supported_config<To, A>();
        return simd_return_type<bool, To, A>::load_aligned(ptr);
    }

    template <class To, class A = default_arch, class From>
    XSIMD_INLINE simd_return_type<std::complex<From>, To, A> load_as(std::complex<From> const* ptr, aligned_mode) noexcept
    {
        detail::static_check_supported_config<To, A>();
        using batch_value_type = typename simd_return_type<std::complex<From>, To, A>::value_type;
        return kernel::load_complex_aligned<A>(ptr, kernel::convert<batch_value_type> {}, A {});
    }

#ifdef XSIMD_ENABLE_XTL_COMPLEX
    template <class To, class A = default_arch, class From, bool i3ec>
    XSIMD_INLINE simd_return_type<xtl::xcomplex<From, From, i3ec>, To, A> load_as(xtl::xcomplex<From, From, i3ec> const* ptr, aligned_mode) noexcept
    {
        detail::static_check_supported_config<To, A>();
        detail::static_check_supported_config<From, A>();
        return load_as<To>(reinterpret_cast<std::complex<From> const*>(ptr), aligned_mode());
    }
#endif

    /**
     * @ingroup batch_data_transfer
     *
     * Creates a batch from the buffer \c ptr and the specifed
     * batch value type \c To. The memory does not need to be aligned.
     * @param ptr the memory buffer to read
     * @return a new batch instance
     */
    template <class To, class A = default_arch, class From>
    XSIMD_INLINE simd_return_type<From, To, A> load_as(From const* ptr, unaligned_mode) noexcept
    {
        using batch_value_type = typename simd_return_type<From, To, A>::value_type;
        detail::static_check_supported_config<To, A>();
        detail::static_check_supported_config<From, A>();
        return kernel::load_unaligned<A>(ptr, kernel::convert<batch_value_type> {}, A {});
    }

    template <class To, class A = default_arch>
    XSIMD_INLINE simd_return_type<bool, To, A> load_as(bool const* ptr, unaligned_mode) noexcept
    {
        return simd_return_type<bool, To, A>::load_unaligned(ptr);
    }

    template <class To, class A = default_arch, class From>
    XSIMD_INLINE simd_return_type<std::complex<From>, To, A> load_as(std::complex<From> const* ptr, unaligned_mode) noexcept
    {
        detail::static_check_supported_config<To, A>();
        detail::static_check_supported_config<From, A>();
        using batch_value_type = typename simd_return_type<std::complex<From>, To, A>::value_type;
        return kernel::load_complex_unaligned<A>(ptr, kernel::convert<batch_value_type> {}, A {});
    }

#ifdef XSIMD_ENABLE_XTL_COMPLEX
    template <class To, class A = default_arch, class From, bool i3ec>
    XSIMD_INLINE simd_return_type<xtl::xcomplex<From, From, i3ec>, To, A> load_as(xtl::xcomplex<From, From, i3ec> const* ptr, unaligned_mode) noexcept
    {
        detail::static_check_supported_config<To, A>();
        detail::static_check_supported_config<From, A>();
        return load_as<To>(reinterpret_cast<std::complex<From> const*>(ptr), unaligned_mode());
    }
#endif

    /**
     * @ingroup batch_data_transfer
     *
     * Creates a batch from the buffer \c ptr. The
     * memory needs to be aligned.
     * @param ptr the memory buffer to read
     * @return a new batch instance
     */
    template <class A = default_arch, class From>
    XSIMD_INLINE batch<From, A> load(From const* ptr, aligned_mode = {}) noexcept
    {
        detail::static_check_supported_config<From, A>();
        return load_as<From, A>(ptr, aligned_mode {});
    }

    /**
     * @ingroup batch_data_transfer
     *
     * Creates a batch from the buffer \c ptr. The
     * memory does not need to be aligned.
     * @param ptr the memory buffer to read
     * @return a new batch instance
     */
    template <class A = default_arch, class From>
    XSIMD_INLINE batch<From, A> load(From const* ptr, unaligned_mode) noexcept
    {
        detail::static_check_supported_config<From, A>();
        return load_as<From, A>(ptr, unaligned_mode {});
    }

    /**
     * @ingroup batch_data_transfer
     *
     * Creates a batch from the buffer \c ptr. The
     * memory needs to be aligned.
     * @param ptr the memory buffer to read
     * @return a new batch instance
     */
    template <class A = default_arch, class From>
    XSIMD_INLINE batch<From, A> load_aligned(From const* ptr) noexcept
    {
        detail::static_check_supported_config<From, A>();
        return load_as<From, A>(ptr, aligned_mode {});
    }

    /**
     * @ingroup batch_data_transfer
     *
     * Creates a batch from the buffer \c ptr. The
     * memory does not need to be aligned.
     * @param ptr the memory buffer to read
     * @return a new batch instance
     */
    template <class A = default_arch, class From>
    XSIMD_INLINE batch<From, A> load_unaligned(From const* ptr) noexcept
    {
        detail::static_check_supported_config<From, A>();
        return load_as<From, A>(ptr, unaligned_mode {});
    }

    /**
     * @ingroup batch_math
     *
     * Computes the natural logarithm of the batch \c x.
     * @param x batch of floating point values.
     * @return the natural logarithm of \c x.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> log(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::log<A>(x, A {});
    }

    /**
     * @ingroup batch_math
     * Computes the base 2 logarithm of the batch \c x.
     * @param x batch of floating point values.
     * @return the base 2 logarithm of \c x.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> log2(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::log2<A>(x, A {});
    }

    /**
     * @ingroup batch_math
     * Computes the base 10 logarithm of the batch \c x.
     * @param x batch of floating point values.
     * @return the base 10 logarithm of \c x.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> log10(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::log10<A>(x, A {});
    }

    /**
     * @ingroup batch_math
     * Computes the natural logarithm of one plus the batch \c x.
     * @param x batch of floating point values.
     * @return the natural logarithm of one plus \c x.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> log1p(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::log1p<A>(x, A {});
    }

    /**
     * @ingroup batch_logical
     *
     * Element-wise lesser than comparison of batches \c x and \c y.
     * @param x batch involved in the comparison.
     * @param y batch involved in the comparison.
     * @return a boolean batch.
     */
    template <class T, class A>
    XSIMD_INLINE batch_bool<T, A> lt(batch<T, A> const& x, batch<T, A> const& y) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return x < y;
    }

    /**
     * @ingroup batch_math
     *
     * Computes the larger values of the batches \c x and \c y.
     * @param x a batch of integer or floating point values.
     * @param y a batch of integer or floating point values.
     * @return a batch of the larger values.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> max(batch<T, A> const& x, batch<T, A> const& y) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::max<A>(x, y, A {});
    }

    /**
     * @ingroup batch_math
     *
     * Computes the smaller values of the batches \c x and \c y.
     * @param x a batch of integer or floating point values.
     * @param y a batch of integer or floating point values.
     * @return a batch of the smaller values.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> min(batch<T, A> const& x, batch<T, A> const& y) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::min<A>(x, y, A {});
    }

    /**
     * @ingroup batch_constant
     *
     * Return a batch of scalars representing positive infinity
     * @return a batch of positive infinity
     */
    template <class B>
    XSIMD_INLINE B minusinfinity() noexcept
    {
        using T = typename B::value_type;
        using A = typename B::arch_type;
        detail::static_check_supported_config<T, A>();
        return B(-std::numeric_limits<T>::infinity());
    }

    /**
     * @ingroup batch_arithmetic
     *
     * Computes the integer modulo of the batch \c x by the batch \c y.
     * @param x batch involved in the modulo.
     * @param y batch involved in the modulo.
     * @return the result of the modulo.
     */
    template <class T, class A>
    XSIMD_INLINE auto mod(batch<T, A> const& x, batch<T, A> const& y) noexcept -> decltype(x % y)
    {
        detail::static_check_supported_config<T, A>();
        return x % y;
    }

    /**
     * @ingroup batch_arithmetic
     *
     * Computes the product of the batches \c x and \c y.
     * @tparam X the actual type of batch.
     * @param x batch involved in the product.
     * @param y batch involved in the product.
     * @return the result of the product.
     */
    template <class T, class A>
    XSIMD_INLINE auto mul(batch<T, A> const& x, batch<T, A> const& y) noexcept -> decltype(x * y)
    {
        detail::static_check_supported_config<T, A>();
        return x * y;
    }

    /**
     * @ingroup batch_rounding
     *
     * Rounds the scalars in \c x to integer values (in floating point format), using
     * the current rounding mode.
     * @param x batch of floating point values.
     * @return the batch of nearest integer values.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> nearbyint(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::nearbyint<A>(x, A {});
    }

    /**
     * @ingroup batch_rounding
     *
     * Rounds the scalars in \c x to integer values (in integer format) using
     * the current rounding mode.
     * @param x batch of floating point values.
     * @return the batch of nearest integer values.
     *
     * @warning For very large values the conversion to int silently overflows.
     */
    template <class T, class A>
    XSIMD_INLINE batch<as_integer_t<T>, A>
    nearbyint_as_int(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::nearbyint_as_int(x, A {});
    }

    /**
     * @ingroup batch_logical
     *
     * Element-wise inequality comparison of batches \c x and \c y.
     * @param x batch involved in the comparison.
     * @param y batch involved in the comparison.
     * @return a boolean batch.
     */
    template <class T, class A>
    XSIMD_INLINE auto neq(batch<T, A> const& x, batch<T, A> const& y) noexcept -> decltype(x != y)
    {
        detail::static_check_supported_config<T, A>();
        return x != y;
    }

    /**
     * @ingroup batch_logical
     *
     * Element-wise inequality comparison of batches of boolean values \c x and \c y.
     * @param x batch of booleans involved in the comparison.
     * @param y batch of booleans involved in the comparison.
     * @return a boolean batch.
     */
    template <class T, class A>
    XSIMD_INLINE auto neq(batch_bool<T, A> const& x, batch_bool<T, A> const& y) noexcept -> decltype(x != y)
    {
        detail::static_check_supported_config<T, A>();
        return x != y;
    }

    /**
     * @ingroup batch_arithmetic
     *
     * Computes the opposite of the batch \c x.
     * @param x batch involved in the operation.
     * @return the opposite of \c x.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> neg(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return -x;
    }

    /**
     * @ingroup batch_math_extra
     *
     * Computes  the next representable  floating-point
     *        value  following  x  in the direction of y
     * @param x batch of floating point values.
     * @param y batch of floating point values.
     * @return \c x raised to the power \c y.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> nextafter(batch<T, A> const& x, batch<T, A> const& y) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::nextafter<A>(x, y, A {});
    }

    /**
     * @ingroup batch_complex
     *
     * Computes the norm of the batch \c x.
     * @param x batch of complex or real values.
     * @return the norm of \c x.
     */
    template <class T, class A>
    XSIMD_INLINE real_batch_type_t<batch<T, A>> norm(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::norm(x, A {});
    }

    /**
     * @ingroup batch_math
     *
     * Returns a complex batch with magnitude \c r and phase angle \c theta.
     * @param r The magnitude of the desired complex result.
     * @param theta The phase angle of the desired complex result.
     * @return \c r exp(i * \c theta).
     */
    template <class T, class A>
    XSIMD_INLINE complex_batch_type_t<batch<T, A>> polar(batch<T, A> const& r, batch<T, A> const& theta = batch<T, A> {}) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::polar<A>(r, theta, A {});
    }

    /**
     * @ingroup batch_arithmetic
     *
     * No-op on \c x.
     * @param x batch involved in the operation.
     * @return \c x.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> pos(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return +x;
    }

    /**
     * @ingroup batch_math
     *
     * Computes the value of the batch \c x raised to the power
     * \c y.
     * @param x batch of floating point values.
     * @param y batch of floating point values.
     * @return \c x raised to the power \c y.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> pow(batch<T, A> const& x, batch<T, A> const& y) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::pow<A>(x, y, A {});
    }

    /**
     * @ingroup batch_math
     *
     * Computes the value of the batch \c x raised to the power
     * \c y.
     * @param x batch of integral values.
     * @param y batch of integral values.
     * @return \c x raised to the power \c y.
     */
    template <class T, class ITy, class A, class = typename std::enable_if<std::is_integral<ITy>::value, void>::type>
    XSIMD_INLINE batch<T, A> pow(batch<T, A> const& x, ITy y) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::ipow<A>(x, y, A {});
    }

    /**
     * @ingroup batch_complex
     *
     * Computes the projection of the batch \c z.
     * @param z batch of complex or real values.
     * @return the projection of \c z.
     */
    template <class T, class A>
    XSIMD_INLINE complex_batch_type_t<batch<T, A>> proj(batch<T, A> const& z) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::proj(z, A {});
    }

    /**
     * @ingroup batch_complex
     *
     * Computes the real part of the batch \c z.
     * @param z batch of complex or real values.
     * @return the argument of \c z.
     */
    template <class T, class A>
    XSIMD_INLINE real_batch_type_t<batch<T, A>> real(batch<T, A> const& z) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::real<A>(z, A {});
    }

    /**
     * @ingroup batch_arithmetic
     *
     * Computes the approximate reciprocal of the batch \c x.
     * The maximum relative error for this approximation is
     * less than 1.5*2^-12.
     * @param x batch of floating point numbers.
     * @return the reciprocal.
     */
    template <class T, class A, class = typename std::enable_if<std::is_floating_point<T>::value, void>::type>
    XSIMD_INLINE batch<T, A> reciprocal(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::reciprocal(x, A {});
    }

    /**
     * @ingroup batch_reducers
     *
     * Generic reducer using only batch operations
     * @param f reducing function, accepting `batch ()(batch, batch)`
     * @param x batch involved in the reduction
     * @return the result of the reduction, as a scalar.
     */
    template <class T, class A, class F>
    XSIMD_INLINE T reduce(F&& f, batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::detail::reduce(std::forward<F>(f), x, std::integral_constant<unsigned, batch<T, A>::size>());
    }

    /**
     * @ingroup batch_reducers
     *
     * Adds all the scalars of the batch \c x.
     * @param x batch involved in the reduction
     * @return the result of the reduction.
     */
    template <class T, class A>
    XSIMD_INLINE T reduce_add(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::reduce_add<A>(x, A {});
    }

    /**
     * @ingroup batch_reducers
     *
     * Max of all the scalars of the batch \c x.
     * @param x batch involved in the reduction
     * @return the result of the reduction.
     */
    template <class T, class A>
    XSIMD_INLINE T reduce_max(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::reduce_max<A>(x, A {});
    }

    /**
     * @ingroup batch_reducers
     *
     * Min of all the scalars of the batch \c x.
     * @param x batch involved in the reduction
     * @return the result of the reduction.
     */
    template <class T, class A>
    XSIMD_INLINE T reduce_min(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::reduce_min<A>(x, A {});
    }

    /**
     * @ingroup batch_math
     *
     * Computes the remainder of dividing \c x by \c y
     * @param x batch of scalar values
     * @param y batch of scalar values
     * @return the result of the addition.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> remainder(batch<T, A> const& x, batch<T, A> const& y) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::remainder<A>(x, y, A {});
    }

    /**
     * @ingroup batch_rounding
     *
     * Rounds the scalars in \c x to integer values (in floating point format), using
     * the current rounding mode.
     * @param x batch of floating point values.
     * @return the batch of rounded values.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> rint(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return nearbyint(x);
    }

    /**
     * @ingroup rotate_left
     *
     * Slide the whole batch to the left by \c n bytes, and reintroduce the
     * slided out elements from the right. This is different from
     * \c rol that rotates each batch element to the left.
     *
     * @tparam N Amount of bytes to rotated to the left.
     * @param x batch of integer values.
     * @return rotated batch.
     */
    template <size_t N, class T, class A>
    XSIMD_INLINE batch<T, A> rotate_left(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::rotate_left<N, A>(x, A {});
    }

    /**
     * @ingroup rotate_right
     *
     * Slide the whole batch to the right by \c n bytes, and reintroduce the
     * slided out elements from the left. This is different from
     * \c rol that rotates each batch element to the left.
     *
     * @tparam N Amount of bytes to rotate to the right.
     * @param x batch of integer values.
     * @return rotated batch.
     */
    template <size_t N, class T, class A>
    XSIMD_INLINE batch<T, A> rotate_right(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::rotate_right<N, A>(x, A {});
    }

    /**
     * @ingroup batch_bitwise
     *
     * Perform a bitwise shift to the left, reintroducing the shifted out bits
     * to the right
     * @param x batch to rotate
     * @param shift scalar amount to shift
     * @return rotated \c x.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> rotl(batch<T, A> const& x, int shift) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::rotl<A>(x, shift, A {});
    }
    template <class T, class A>
    XSIMD_INLINE batch<T, A> rotl(batch<T, A> const& x, batch<T, A> const& shift) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::rotl<A>(x, shift, A {});
    }

    /**
     * @ingroup batch_bitwise
     *
     * Perform a bitwise shift to the right, reintroducing the shifted out bits
     * to the left.
     * @param x batch to rotate
     * @param shift scalar amount to shift
     * @return rotated \c x.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> rotr(batch<T, A> const& x, int shift) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::rotr<A>(x, shift, A {});
    }
    template <class T, class A>
    XSIMD_INLINE batch<T, A> rotr(batch<T, A> const& x, batch<T, A> const& shift) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::rotr<A>(x, shift, A {});
    }

    /**
     * @ingroup batch_rounding
     *
     * Computes the batch of nearest integer values to scalars in \c x (in
     * floating point format), rounding halfway cases away from zero, regardless
     * of the current rounding mode.
     * @param x batch of flaoting point values.
     * @return the batch of nearest integer values.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> round(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::round<A>(x, A {});
    }

    /**
     * @ingroup batch_math
     *
     * Computes an estimate of the inverse square root of the batch \c x.
     *
     * @warning Unlike most xsimd function, this does not return the same result as the
     * equivalent scalar operation, trading accuracy for speed.
     *
     * @param x batch of floating point values.
     * @return the inverse square root of \c x.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> rsqrt(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::rsqrt<A>(x, A {});
    }

    /**
     * @ingroup batch_arithmetic
     *
     * Computes the saturate sum of the batch \c x and the batch \c y.

     * @tparam X the actual type of batch.
     * @param x batch involved in the saturated addition.
     * @param y batch involved in the saturated addition.
     * @return the result of the saturated addition.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> sadd(batch<T, A> const& x, batch<T, A> const& y) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::sadd<A>(x, y, A {});
    }

    /**
     * @ingroup batch_miscellaneous
     *
     * Ternary operator for batches: selects values from the batches \c true_br or \c false_br
     * depending on the boolean values in the constant batch \c cond. Equivalent to
     * \code{.cpp}
     * for(std::size_t i = 0; i < N; ++i)
     *     res[i] = cond[i] ? true_br[i] : false_br[i];
     * \endcode
     * @param cond batch condition.
     * @param true_br batch values for truthy condition.
     * @param false_br batch value for falsy condition.
     * @return the result of the selection.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> select(batch_bool<T, A> const& cond, batch<T, A> const& true_br, batch<T, A> const& false_br) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::select<A>(cond, true_br, false_br, A {});
    }

    /**
     * @ingroup batch_miscellaneous
     *
     * Ternary operator for batches: selects values from the batches \c true_br or \c false_br
     * depending on the boolean values in the constant batch \c cond. Equivalent to
     * \code{.cpp}
     * for(std::size_t i = 0; i < N; ++i)
     *     res[i] = cond[i] ? true_br[i] : false_br[i];
     * \endcode
     * @param cond batch condition.
     * @param true_br batch values for truthy condition.
     * @param false_br batch value for falsy condition.
     * @return the result of the selection.
     */
    template <class T, class A>
    XSIMD_INLINE batch<std::complex<T>, A> select(batch_bool<T, A> const& cond, batch<std::complex<T>, A> const& true_br, batch<std::complex<T>, A> const& false_br) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::select<A>(cond, true_br, false_br, A {});
    }

    /**
     * @ingroup batch_miscellaneous
     *
     * Ternary operator for batches: selects values from the batches \c true_br or \c false_br
     * depending on the boolean values in the constant batch \c cond. Equivalent to
     * \code{.cpp}
     * for(std::size_t i = 0; i < N; ++i)
     *     res[i] = cond[i] ? true_br[i] : false_br[i];
     * \endcode
     * @param cond constant batch condition.
     * @param true_br batch values for truthy condition.
     * @param false_br batch value for falsy condition.
     * @return the result of the selection.
     */
    template <class T, class A, bool... Values>
    XSIMD_INLINE batch<T, A> select(batch_bool_constant<T, A, Values...> const& cond, batch<T, A> const& true_br, batch<T, A> const& false_br) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::select<A>(cond, true_br, false_br, A {});
    }

    /**
     * @ingroup batch_data_transfer
     *
     * Combine elements from \c x and \c y according to selector \c mask
     * @param x batch
     * @param y batch
     * @param mask constant batch mask of integer elements of the same size as
     * element of \c x and \c y. Each element of the mask index the vector that
     * would be formed by the concatenation of \c x and \c y. For instance
     * \code{.cpp}
     * batch_constant<uint32_t, sse2, 0, 4, 3, 7>
     * \endcode
     * Picks \c x[0], \c y[0], \c x[3], \c y[3]
     *
     * @return combined batch
     */
    template <class T, class A, class Vt, Vt... Values>
    XSIMD_INLINE typename std::enable_if<std::is_arithmetic<T>::value, batch<T, A>>::type
    shuffle(batch<T, A> const& x, batch<T, A> const& y, batch_constant<Vt, A, Values...> mask) noexcept
    {
        static_assert(sizeof(T) == sizeof(Vt), "consistent mask");
        detail::static_check_supported_config<T, A>();
        return kernel::shuffle<A>(x, y, mask, A {});
    }

    /**
     * @ingroup batch_miscellaneous
     *
     * Computes the sign of \c x
     * @param x batch
     * @return -1 for each negative element, -1 or +1 for each null element and +1 for each element
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> sign(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::sign<A>(x, A {});
    }

    /**
     * @ingroup batch_miscellaneous
     *
     * Computes the sign of \c x, assuming x doesn't have any zero
     * @param x batch
     * @return -1 for each negative element, -1 or +1 for each null element and +1 for each element
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> signnz(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::signnz<A>(x, A {});
    }

    /**
     * @ingroup batch_trigo
     *
     * Computes the sine of the batch \c x.
     * @param x batch of floating point values.
     * @return the sine of \c x.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> sin(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::sin<A>(x, A {});
    }

    /**
     * @ingroup batch_trigo
     *
     * Computes the sine and the cosine of the batch \c x. This method is faster
     * than calling sine and cosine independently.
     * @param x batch of floating point values.
     * @return a pair containing the sine then the cosine of  batch \c x
     */
    template <class T, class A>
    XSIMD_INLINE std::pair<batch<T, A>, batch<T, A>> sincos(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::sincos<A>(x, A {});
    }

    /**
     * @ingroup batch_trigo
     *
     * Computes the hyperbolic sine of the batch \c x.
     * @param x batch of floating point values.
     * @return the hyperbolic sine of \c x.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> sinh(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::sinh<A>(x, A {});
    }

    /**
     * @ingroup batch_data_transfer
     *
     * Slide the whole batch to the left by \c n bytes. This is different from
     * \c bitwise_lshift that shifts each batch element to the left.
     *
     * @tparam N Amount of bytes to slide to the left.
     * @param x batch of integer values.
     * @return slided batch.
     */
    template <size_t N, class T, class A>
    XSIMD_INLINE batch<T, A> slide_left(batch<T, A> const& x) noexcept
    {
        static_assert(std::is_integral<T>::value, "can only slide batch of integers");
        detail::static_check_supported_config<T, A>();
        return kernel::slide_left<N, A>(x, A {});
    }

    /**
     * @ingroup batch_data_transfer
     *
     * Slide the whole batch to the right by \c N bytes. This is different from
     * \c bitwise_rshift that shifts each batch element to the right.
     *
     * @tparam N Amount of bytes to slide to the right.
     * @param x batch of integer values.
     * @return slided batch.
     */
    template <size_t N, class T, class A>
    XSIMD_INLINE batch<T, A> slide_right(batch<T, A> const& x) noexcept
    {
        static_assert(std::is_integral<T>::value, "can only slide batch of integers");
        detail::static_check_supported_config<T, A>();
        return kernel::slide_right<N, A>(x, A {});
    }

    /**
     * @ingroup batch_math
     *
     * Computes the square root of the batch \c x.
     * @param x batch of floating point values.
     * @return the square root of \c x.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> sqrt(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::sqrt<A>(x, A {});
    }

    /**
     * @ingroup batch_arithmetic
     *
     * Computes the saturate difference of the batch \c x and the batch \c y.
     * @tparam X the actual type of batch.
     * @param x batch involved in the saturated difference.
     * @param y batch involved in the saturated difference.
     * @return the result of the saturated difference.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> ssub(batch<T, A> const& x, batch<T, A> const& y) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::ssub<A>(x, y, A {});
    }

    /**
     * @ingroup batch_data_transfer
     *
     * Copy content of batch \c src to the buffer \c dst. The
     * memory needs to be aligned.
     * @param dst the memory buffer to write to
     * @param src the batch to copy
     */
    template <class To, class A = default_arch, class From>
    XSIMD_INLINE void store_as(To* dst, batch<From, A> const& src, aligned_mode) noexcept
    {
        detail::static_check_supported_config<From, A>();
        kernel::store_aligned<A>(dst, src, A {});
    }

    template <class A = default_arch, class From>
    XSIMD_INLINE void store_as(bool* dst, batch_bool<From, A> const& src, aligned_mode) noexcept
    {
        detail::static_check_supported_config<From, A>();
        kernel::store<A>(src, dst, A {});
    }

    template <class To, class A = default_arch, class From>
    XSIMD_INLINE void store_as(std::complex<To>* dst, batch<std::complex<From>, A> const& src, aligned_mode) noexcept
    {
        detail::static_check_supported_config<std::complex<From>, A>();
        kernel::store_complex_aligned<A>(dst, src, A {});
    }

#ifdef XSIMD_ENABLE_XTL_COMPLEX
    template <class To, class A = default_arch, class From, bool i3ec>
    XSIMD_INLINE void store_as(xtl::xcomplex<To, To, i3ec>* dst, batch<std::complex<From>, A> const& src, aligned_mode) noexcept
    {
        store_as(reinterpret_cast<std::complex<To>*>(dst), src, aligned_mode());
    }
#endif

    /**
     * @ingroup batch_data_transfer
     *
     * Copy content of batch \c src to the buffer \c dst. The
     * memory does not need to be aligned.
     * @param dst the memory buffer to write to
     * @param src the batch to copy
     */
    template <class To, class A = default_arch, class From>
    XSIMD_INLINE void store_as(To* dst, batch<From, A> const& src, unaligned_mode) noexcept
    {
        detail::static_check_supported_config<From, A>();
        kernel::store_unaligned<A>(dst, src, A {});
    }

    template <class A = default_arch, class From>
    XSIMD_INLINE void store_as(bool* dst, batch_bool<From, A> const& src, unaligned_mode) noexcept
    {
        detail::static_check_supported_config<From, A>();
        kernel::store<A>(src, dst, A {});
    }

    template <class To, class A = default_arch, class From>
    XSIMD_INLINE void store_as(std::complex<To>* dst, batch<std::complex<From>, A> const& src, unaligned_mode) noexcept
    {
        detail::static_check_supported_config<std::complex<From>, A>();
        kernel::store_complex_unaligned<A>(dst, src, A {});
    }

#ifdef XSIMD_ENABLE_XTL_COMPLEX
    template <class To, class A = default_arch, class From, bool i3ec>
    XSIMD_INLINE void store_as(xtl::xcomplex<To, To, i3ec>* dst, batch<std::complex<From>, A> const& src, unaligned_mode) noexcept
    {
        detail::static_check_supported_config<std::complex<From>, A>();
        store_as(reinterpret_cast<std::complex<To>*>(dst), src, unaligned_mode());
    }
#endif

    /**
     * @ingroup batch_data_transfer
     *
     * Copy content of batch \c val to the buffer \c mem. The
     * memory does not need to be aligned.
     * @param mem the memory buffer to write to
     * @param val the batch to copy from
     */
    template <class A, class T>
    XSIMD_INLINE void store(T* mem, batch<T, A> const& val, aligned_mode = {}) noexcept
    {
        store_as<T, A>(mem, val, aligned_mode {});
    }

    /**
     * @ingroup batch_data_transfer
     *
     * Copy content of batch \c val to the buffer \c mem. The
     * memory does not need to be aligned.
     * @param mem the memory buffer to write to
     * @param val the batch to copy from
     */
    template <class A, class T>
    XSIMD_INLINE void store(T* mem, batch<T, A> const& val, unaligned_mode) noexcept
    {
        store_as<T, A>(mem, val, unaligned_mode {});
    }

    /**
     * @ingroup batch_data_transfer
     *
     * Copy content of batch \c val to the buffer \c mem. The
     * memory needs to be aligned.
     * @param mem the memory buffer to write to
     * @param val the batch to copy from
     */
    template <class A, class T>
    XSIMD_INLINE void store_aligned(T* mem, batch<T, A> const& val) noexcept
    {
        store_as<T, A>(mem, val, aligned_mode {});
    }

    /**
     * @ingroup batch_data_transfer
     *
     * Copy content of batch \c val to the buffer \c mem. The
     * memory does not need to be aligned.
     * @param mem the memory buffer to write to
     * @param val the batch to copy
     */
    template <class A, class T>
    XSIMD_INLINE void store_unaligned(T* mem, batch<T, A> const& val) noexcept
    {
        store_as<T, A>(mem, val, unaligned_mode {});
    }

    /**
     * @ingroup batch_arithmetic
     *
     * Computes the difference between \c x and \c y
     * @tparam X the actual type of batch.
     * @param x scalar or batch of scalars
     * @param y scalar or batch of scalars
     * @return the difference between \c x and \c y
     */
    template <class T, class A>
    XSIMD_INLINE auto sub(batch<T, A> const& x, batch<T, A> const& y) noexcept -> decltype(x - y)
    {
        detail::static_check_supported_config<T, A>();
        return x - y;
    }

    /**
     * @ingroup batch_data_transfer
     *
     * Rearrange elements from \c x according to constant mask \c mask
     * @param x batch
     * @param mask constant batch mask of integer elements of the same size as
     * element of \c x
     * @return swizzled batch
     */
    template <class T, class A, class Vt, Vt... Values>
    XSIMD_INLINE typename std::enable_if<std::is_arithmetic<T>::value, batch<T, A>>::type
    swizzle(batch<T, A> const& x, batch_constant<Vt, A, Values...> mask) noexcept
    {
        static_assert(sizeof(T) == sizeof(Vt), "consistent mask");
        detail::static_check_supported_config<T, A>();
        return kernel::swizzle<A>(x, mask, A {});
    }
    template <class T, class A, class Vt, Vt... Values>
    XSIMD_INLINE batch<std::complex<T>, A> swizzle(batch<std::complex<T>, A> const& x, batch_constant<Vt, A, Values...> mask) noexcept
    {
        static_assert(sizeof(T) == sizeof(Vt), "consistent mask");
        detail::static_check_supported_config<T, A>();
        return kernel::swizzle<A>(x, mask, A {});
    }

    /**
     * @ingroup batch_data_transfer
     *
     * Rearrange elements from \c x according to mask \c mask
     * @param x batch
     * @param mask batch mask of integer elements of the same size as
     * element of \c x
     * @return swizzled batch
     */
    template <class T, class A, class Vt>
    XSIMD_INLINE typename std::enable_if<std::is_arithmetic<T>::value, batch<T, A>>::type
    swizzle(batch<T, A> const& x, batch<Vt, A> mask) noexcept
    {
        static_assert(sizeof(T) == sizeof(Vt), "consistent mask");
        detail::static_check_supported_config<T, A>();
        return kernel::swizzle<A>(x, mask, A {});
    }

    template <class T, class A, class Vt>
    XSIMD_INLINE batch<std::complex<T>, A> swizzle(batch<std::complex<T>, A> const& x, batch<Vt, A> mask) noexcept
    {
        static_assert(sizeof(T) == sizeof(Vt), "consistent mask");
        detail::static_check_supported_config<T, A>();
        return kernel::swizzle<A>(x, mask, A {});
    }

    /**
     * @ingroup batch_trigo
     *
     * Computes the tangent of the batch \c x.
     * @param x batch of floating point values.
     * @return the tangent of \c x.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> tan(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::tan<A>(x, A {});
    }

    /**
     * @ingroup batch_trigo
     *
     * Computes the hyperbolic tangent of the batch \c x.
     * @param x batch of floating point values.
     * @return the hyperbolic tangent of \c x.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> tanh(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::tanh<A>(x, A {});
    }

    /**
     * @ingroup batch_math_extra
     *
     * Computes the gamma function of the batch \c x.
     * @param x batch of floating point values.
     * @return the gamma function of \c x.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> tgamma(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::tgamma<A>(x, A {});
    }

    /**
     * @ingroup batch_conversion
     *
     * Perform a conversion from \c i to a value of an floating point type of the same size as \c T.
     * This is equivalent to \c batch_cast<as_float_t<T>>(i)
     * @param i batch of integers.
     * @return \c i converted to a value of an floating point type of the same size as \c T
     */
    template <class T, class A>
    XSIMD_INLINE batch<as_float_t<T>, A> to_float(batch<T, A> const& i) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return batch_cast<as_float_t<T>>(i);
    }

    /**
     * @ingroup batch_conversion
     *
     * Perform a conversion from \c x to a value of an integer type of the same size as \c T
     * This is equivalent to \c batch_cast<as_integer_t<T>>(x)
     * @param x batch.
     * @return \c x converted to a value of an integer type of the same size as \c T
     */
    template <class T, class A>
    XSIMD_INLINE batch<as_integer_t<T>, A> to_int(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return batch_cast<as_integer_t<T>>(x);
    }

    /**
     * @ingroup batch_rounding
     *
     * Computes the batch of nearest integer values not greater in magnitude
     * than scalars in \c x.
     * @param x batch of floating point values.
     * @return the batch of nearest integer values not greater in magnitude than \c x.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> trunc(batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::trunc<A>(x, A {});
    }

    /**
     * @ingroup batch_data_transfer
     *
     * Unpack and interleave data from the HIGH half of batches \c x and \c y.
     * Store the results in the Return value.
     * @param x a batch of integer or floating point or double precision values.
     * @param y a batch of integer or floating point or double precision values.
     * @return a batch of the high part of shuffled values.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> zip_hi(batch<T, A> const& x, batch<T, A> const& y) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::zip_hi<A>(x, y, A {});
    }

    /**
     * @ingroup batch_data_transfer
     *
     * Unpack and interleave data from the LOW half of batches \c x and \c y.
     * Store the results in the Return value.
     * @param x a batch of integer or floating point or double precision values.
     * @param y a batch of integer or floating point or double precision values.
     * @return a batch of the low part of shuffled values.
     */
    template <class T, class A>
    XSIMD_INLINE batch<T, A> zip_lo(batch<T, A> const& x, batch<T, A> const& y) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::zip_lo<A>(x, y, A {});
    }

    /**
     * @ingroup batch_conversion
     *
     * Cast a \c batch_bool of \c T into a \c batch of the same type using the
     * following rule: if an element of \c self is true, it maps to -1 in the
     * returned integral batch, otherwise it maps to 0.
     *
     * @param self batch_bool of \c T
     * @return \c self cast to a \c batch of \c T
     */
    template <class T, class A, typename std::enable_if<std::is_integral<T>::value, int>::type = 3>
    XSIMD_INLINE batch<T, A> bitwise_cast(batch_bool<T, A> const& self) noexcept
    {
        T z(0);
        detail::static_check_supported_config<T, A>();
        return select(self, batch<T, A>(T(~z)), batch<T, A>(z));
    }

    template <class T, class A, typename std::enable_if<std::is_floating_point<T>::value, int>::type = 3>
    XSIMD_INLINE batch<T, A> bitwise_cast(batch_bool<T, A> const& self) noexcept
    {
        T z0(0), z1(0);
        using int_type = as_unsigned_integer_t<T>;
        int_type value(~int_type(0));
        std::memcpy(&z1, &value, sizeof(int_type));
        detail::static_check_supported_config<T, A>();
        return select(self, batch<T, A>(z1), batch<T, A>(z0));
    }

    /**
     * @ingroup batch_bool_reducers
     *
     * Returns true if all the boolean values in the batch are true,
     * false otherwise.
     * @param x the batch to reduce.
     * @return a boolean scalar.
     */
    template <class T, class A>
    XSIMD_INLINE bool all(batch_bool<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::all<A>(x, A {});
    }

    /**
     * @ingroup batch_bool_reducers
     *
     * Return true if any of the boolean values in the batch is true,
     * false otherwise.
     * @param x the batch to reduce.
     * @return a boolean scalar.
     */
    template <class T, class A>
    XSIMD_INLINE bool any(batch_bool<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return kernel::any<A>(x, A {});
    }

    /**
     * @ingroup batch_bool_reducers
     *
     * Return true if none of the boolean values in the batch is true,
     * false otherwise.
     * @param x the batch to reduce.
     * @return a boolean scalar.
     */
    template <class T, class A>
    XSIMD_INLINE bool none(batch_bool<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        return !xsimd::any(x);
    }

    /**
     * @ingroup batch_miscellaneous
     *
     * Dump the content of batch \c x to stream \c o
     * @param o the stream where the batch is dumped
     * @param x batch to dump.
     * @return a reference to \c o
     */
    template <class T, class A>
    XSIMD_INLINE std::ostream& operator<<(std::ostream& o, batch<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        constexpr auto size = batch<T, A>::size;
        alignas(A::alignment()) T buffer[size];
        x.store_aligned(&buffer[0]);
        o << '(';
        for (std::size_t i = 0; i < size - 1; ++i)
            o << buffer[i] << ", ";
        return o << buffer[size - 1] << ')';
    }

    /**
     * @ingroup batch_miscellaneous
     *
     * Dump the content of batch \c x to stream \c o
     * @param o the stream where the batch is dumped
     * @param x batch to dump.
     * @return a reference to \c o
     */
    template <class T, class A>
    XSIMD_INLINE std::ostream& operator<<(std::ostream& o, batch_bool<T, A> const& x) noexcept
    {
        detail::static_check_supported_config<T, A>();
        constexpr auto size = batch_bool<T, A>::size;
        alignas(A::alignment()) bool buffer[size];
        x.store_aligned(&buffer[0]);
        o << '(';
        for (std::size_t i = 0; i < size - 1; ++i)
            o << buffer[i] << ", ";
        return o << buffer[size - 1] << ')';
    }
}

#endif
