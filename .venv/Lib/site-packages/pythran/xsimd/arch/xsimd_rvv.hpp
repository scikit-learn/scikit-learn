/***************************************************************************

 * Copyright (c) Rivos Inc.                                                 *
 *                                                                          *
 * Distributed under the terms of the BSD 3-Clause License.                 *
 *                                                                          *
 * The full license is in the file LICENSE, distributed with this software. *
 ****************************************************************************/

#ifndef XSIMD_RVV_HPP
#define XSIMD_RVV_HPP

#include <complex>
#include <type_traits>
#include <utility>

#include "../types/xsimd_rvv_register.hpp"
#include "xsimd_constants.hpp"

// This set of macros allows the synthesis of identifiers using a template and
// variable macro arguments.  A single template can then be used by multiple
// macros, or multiple instances of a macro to define the same logic for
// different data types.
//
// First some logic to paste text together...
//
#define XSIMD_RVV_JOIN_(x, y) x##y
#define XSIMD_RVV_JOIN(x, y) XSIMD_RVV_JOIN_(x, y)
#define XSIMD_RVV_PREFIX_T(T, S, then) XSIMD_RVV_JOIN(T, then)
#define XSIMD_RVV_PREFIX_S(T, S, then) XSIMD_RVV_JOIN(S, then)
#define XSIMD_RVV_PREFIX_M(T, S, then) XSIMD_RVV_JOIN(m1, then)
#define XSIMD_RVV_PREFIX(T, S, then) then
//
// XSIMD_RVV_IDENTIFIER accepts type and size parameters, and a template for
// the identifier.  The template is a comma-separated list of alternating
// literal and parameter segments.  Each parameter is appended to XSIMD_RVV_PREFIX to
// form a new macro name which decides which parameter should be inserted.
// Then a literal segment is inserted after that.  Empty literals are used to
// join two or more variables together.
//
#define XSIMD_RVV_IDENTIFIER9(T, S, t, ...) t
#define XSIMD_RVV_IDENTIFIER8(T, S, t, p, ...) XSIMD_RVV_JOIN(t, XSIMD_RVV_PREFIX##p(T, S, XSIMD_RVV_IDENTIFIER9(T, S, __VA_ARGS__)))
#define XSIMD_RVV_IDENTIFIER7(T, S, t, p, ...) XSIMD_RVV_JOIN(t, XSIMD_RVV_PREFIX##p(T, S, XSIMD_RVV_IDENTIFIER8(T, S, __VA_ARGS__)))
#define XSIMD_RVV_IDENTIFIER6(T, S, t, p, ...) XSIMD_RVV_JOIN(t, XSIMD_RVV_PREFIX##p(T, S, XSIMD_RVV_IDENTIFIER7(T, S, __VA_ARGS__)))
#define XSIMD_RVV_IDENTIFIER5(T, S, t, p, ...) XSIMD_RVV_JOIN(t, XSIMD_RVV_PREFIX##p(T, S, XSIMD_RVV_IDENTIFIER6(T, S, __VA_ARGS__)))
#define XSIMD_RVV_IDENTIFIER4(T, S, t, p, ...) XSIMD_RVV_JOIN(t, XSIMD_RVV_PREFIX##p(T, S, XSIMD_RVV_IDENTIFIER5(T, S, __VA_ARGS__)))
#define XSIMD_RVV_IDENTIFIER3(T, S, t, p, ...) XSIMD_RVV_JOIN(t, XSIMD_RVV_PREFIX##p(T, S, XSIMD_RVV_IDENTIFIER4(T, S, __VA_ARGS__)))
#define XSIMD_RVV_IDENTIFIER2(T, S, t, p, ...) XSIMD_RVV_JOIN(t, XSIMD_RVV_PREFIX##p(T, S, XSIMD_RVV_IDENTIFIER3(T, S, __VA_ARGS__)))
#define XSIMD_RVV_IDENTIFIER1(T, S, t, p, ...) XSIMD_RVV_JOIN(t, XSIMD_RVV_PREFIX##p(T, S, XSIMD_RVV_IDENTIFIER2(T, S, __VA_ARGS__)))
#define XSIMD_RVV_IDENTIFIER0(T, S, t, p, ...) XSIMD_RVV_JOIN(t, XSIMD_RVV_PREFIX##p(T, S, XSIMD_RVV_IDENTIFIER1(T, S, __VA_ARGS__)))
//
// UNBRACKET and REPARSE force the preprocessor to handle expansion in a
// specific order.  XSIMD_RVV_UNBRACKET strips the parentheses from the template
// (which were necessary to keep the template as a single, named macro
// parameter up to this point).  XSIMD_RVV_ARG_LIST then forms the new parameter list
// to pass to XSIMD_RVV_IDENTIFIER0, with trailing commas to ensure the unrolled
// XSIMD_RVV_IDENTIFIER loop runs to completion adding empty strings.
//
// However XSIMD_RVV_IDENTIFIER0 is not expanded immediately because it does not
// match a function-like macro in this pass.  XSIMD_RVV_REPARSE forces another
// evaluation after the expansion of XSIMD_RVV_ARG_LIST, where XSIMD_RVV_IDENTIFIER0 will
// now match as a function-like macro, and the cycle of substitutions and
// insertions can begin.
//
#define XSIMD_RVV_REPARSE(v) (v)
#define XSIMD_RVV_UNBRACKET(...) __VA_ARGS__
#define XSIMD_RVV_ARG_LIST(T, S, name) (T, S, XSIMD_RVV_UNBRACKET name, , , , , , , , , , , , , , , , , , , , , )
#define XSIMD_RVV_IDENTIFIER(T, S, name) XSIMD_RVV_REPARSE(XSIMD_RVV_IDENTIFIER0 XSIMD_RVV_ARG_LIST(T, S, name))
//
// To avoid comma-counting bugs, replace the variable references with macros
// which include enough commas to keep proper phase, and then use no commas at
// all in the templates.
//
#define XSIMD_RVV_T , _T,
#define XSIMD_RVV_S , _S,
#define XSIMD_RVV_M , _M,
#define XSIMD_RVV_TSM XSIMD_RVV_T XSIMD_RVV_S XSIMD_RVV_M

// XSIMD_RVV_OVERLOAD, below, expands to a head section, a number of body sections
// (depending on which types are supported), and a tail section.  Different
// variants of these sections are implemented with different suffixes on the
// three macro names XSIMD_RVV_WRAPPER_HEAD, XSIMD_RVV_WRAPPER, and XSIMD_RVV_WRAPPER_TAIL and
// specified as an argument to XSIMD_RVV_OVERLOAD (the empty string is the default,
// but still needs an extra comma to hold its place).
//
// The default XSIMD_RVV_WRAPPER_HEAD provides a class containing convenient names
// for the function signature argument(s) to XSIMD_RVV_OVERLOAD.  That signature can
// also reference the template argument T, because it's a text substitution
// into the template.
#define XSIMD_RVV_WRAPPER_HEAD(NAME, SIGNATURE, ...)                      \
    namespace NAME##_cruft                                                \
    {                                                                     \
        template <class T>                                                \
        struct ctx                                                        \
        {                                                                 \
            static constexpr size_t width = XSIMD_RVV_BITS;               \
            static constexpr size_t vl = width / (sizeof(T) * 8);         \
            using vec = rvv_reg_t<T, width>;                              \
            using uvec = rvv_reg_t<as_unsigned_relaxed_t<T>, width>;      \
            using svec = rvv_reg_t<as_signed_relaxed_t<T>, width>;        \
            using fvec = rvv_reg_t<as_float_relaxed_t<T>, width>;         \
            using bvec = rvv_bool_t<T, width>;                            \
            using scalar_vec = rvv_reg_t<T, types::detail::rvv_width_m1>; \
            using wide_vec = rvv_reg_t<T, width * 2>;                     \
            using narrow_vec = rvv_reg_t<T, width / 2>;                   \
            using type = SIGNATURE;                                       \
        };                                                                \
        template <class T>                                                \
        using sig_t = typename ctx<T>::type;                              \
        template <class K, class T>                                       \
        struct impl                                                       \
        {                                                                 \
            void operator()() const noexcept {};                          \
        };                                                                \
        template <class K>                                                \
        using impl_t = impl<K, sig_t<K>>;

#define XSIMD_RVV_WRAPPER_HEAD_NOVL(...) XSIMD_RVV_WRAPPER_HEAD(__VA_ARGS__)
#define XSIMD_RVV_WRAPPER_HEAD_DROP_1ST(...) XSIMD_RVV_WRAPPER_HEAD(__VA_ARGS__)
#define XSIMD_RVV_WRAPPER_HEAD_DROP_1ST_CUSTOM_ARGS(...) XSIMD_RVV_WRAPPER_HEAD(__VA_ARGS__)
#define XSIMD_RVV_WRAPPER_HEAD_DROP_1ST_CUSTOM_ARGS_NOVL(...) XSIMD_RVV_WRAPPER_HEAD(__VA_ARGS__)

// The body of the wrapper defines a functor (because partial specialisation of
// functions is not legal) which forwards its arguments to the named intrinsic
// with a few manipulations.  In general, vector types are handled as
// rvv_reg_t<> and rely on the conversion operators in that class for
// compatibility with the intrinsics.
//
// The function signature is not mentioned here.  Instead it's provided in the
// tail code as the template argument for which this is a specialisation, which
// overcomes the problem of converting a function signature type to an argument
// list to pass to another function.
//
#define XSIMD_RVV_WRAPPER(KEY, CALLEE, ...)                   \
    template <class Ret, class... Args>                       \
    struct impl<KEY, Ret(Args...)>                            \
    {                                                         \
        using ctx = ctx<KEY>;                                 \
        constexpr Ret operator()(Args... args) const noexcept \
        {                                                     \
            return CALLEE(args..., ctx::vl);                  \
        };                                                    \
    };
#define XSIMD_RVV_WRAPPER_NOVL(KEY, CALLEE, ...)              \
    template <class Ret, class... Args>                       \
    struct impl<KEY, Ret(Args...)>                            \
    {                                                         \
        constexpr Ret operator()(Args... args) const noexcept \
        {                                                     \
            return CALLEE(args...);                           \
        };                                                    \
    };
#define XSIMD_RVV_WRAPPER_DROP_1ST(KEY, CALLEE, ...)                 \
    template <class Ret, class First, class... Args>                 \
    struct impl<KEY, Ret(First, Args...)>                            \
    {                                                                \
        using ctx = ctx<KEY>;                                        \
        constexpr Ret operator()(First, Args... args) const noexcept \
        {                                                            \
            return CALLEE(args..., ctx::vl);                         \
        };                                                           \
    };
#define XSIMD_RVV_WRAPPER_DROP_1ST_CUSTOM_ARGS(KEY, CALLEE, SIGNATURE, ...) \
    template <class Ret, class First, class... Args>                        \
    struct impl<KEY, Ret(First, Args...)>                                   \
    {                                                                       \
        using ctx = ctx<KEY>;                                               \
        constexpr Ret operator()(First, Args... args) const noexcept        \
        {                                                                   \
            return CALLEE(__VA_ARGS__, ctx::vl);                            \
        };                                                                  \
    };
#define XSIMD_RVV_WRAPPER_DROP_1ST_CUSTOM_ARGS_NOVL(KEY, CALLEE, SIGNATURE, ...) \
    template <class Ret, class First, class... Args>                             \
    struct impl<KEY, Ret(First, Args...)>                                        \
    {                                                                            \
        constexpr Ret operator()(First, Args... args) const noexcept             \
        {                                                                        \
            return CALLEE(__VA_ARGS__);                                          \
        };                                                                       \
    };

// This part folds all the above templates down into a single functor instance
// with all the different function signatures available under the one name.
// Not all of the base classes necessarily contain useful code, but there's a
// default implementation so that filtering them out isn't really necessary.
#define XSIMD_RVV_WRAPPER_TAIL(NAME, ...)                     \
    } /* namespace NAME##_cruft */                            \
    static constexpr struct : NAME##_cruft::impl_t<int8_t>,   \
                              NAME##_cruft::impl_t<uint8_t>,  \
                              NAME##_cruft::impl_t<int16_t>,  \
                              NAME##_cruft::impl_t<uint16_t>, \
                              NAME##_cruft::impl_t<int32_t>,  \
                              NAME##_cruft::impl_t<uint32_t>, \
                              NAME##_cruft::impl_t<int64_t>,  \
                              NAME##_cruft::impl_t<uint64_t>, \
                              NAME##_cruft::impl_t<float>,    \
                              NAME##_cruft::impl_t<double>    \
    {                                                         \
        using NAME##_cruft::impl_t<int8_t>::operator();       \
        using NAME##_cruft::impl_t<uint8_t>::operator();      \
        using NAME##_cruft::impl_t<int16_t>::operator();      \
        using NAME##_cruft::impl_t<uint16_t>::operator();     \
        using NAME##_cruft::impl_t<int32_t>::operator();      \
        using NAME##_cruft::impl_t<uint32_t>::operator();     \
        using NAME##_cruft::impl_t<int64_t>::operator();      \
        using NAME##_cruft::impl_t<uint64_t>::operator();     \
        using NAME##_cruft::impl_t<float>::operator();        \
        using NAME##_cruft::impl_t<double>::operator();       \
    } NAME {};
#define XSIMD_RVV_WRAPPER_TAIL_NOVL(...) XSIMD_RVV_WRAPPER_TAIL(__VA_ARGS__)
#define XSIMD_RVV_WRAPPER_TAIL_DROP_1ST(...) XSIMD_RVV_WRAPPER_TAIL(__VA_ARGS__)
#define XSIMD_RVV_WRAPPER_TAIL_DROP_1ST_CUSTOM_ARGS(...) XSIMD_RVV_WRAPPER_TAIL(__VA_ARGS__)
#define XSIMD_RVV_WRAPPER_TAIL_DROP_1ST_CUSTOM_ARGS_NOVL(...) XSIMD_RVV_WRAPPER_TAIL(__VA_ARGS__)

// clang-format off

#define XSIMD_RVV_OVERLOAD_head(my_name, variant, ...) \
    XSIMD_RVV_WRAPPER_HEAD##variant(my_name, __VA_ARGS__)
#define XSIMD_RVV_OVERLOAD_i(name, variant, ...)                                        \
    XSIMD_RVV_WRAPPER##variant(int8_t, XSIMD_RVV_IDENTIFIER(i, 8, name), __VA_ARGS__)   \
    XSIMD_RVV_WRAPPER##variant(int16_t, XSIMD_RVV_IDENTIFIER(i, 16, name), __VA_ARGS__) \
    XSIMD_RVV_WRAPPER##variant(int32_t, XSIMD_RVV_IDENTIFIER(i, 32, name), __VA_ARGS__) \
    XSIMD_RVV_WRAPPER##variant(int64_t, XSIMD_RVV_IDENTIFIER(i, 64, name), __VA_ARGS__)
#define XSIMD_RVV_OVERLOAD_u(name, variant, ...)                                         \
    XSIMD_RVV_WRAPPER##variant(uint8_t, XSIMD_RVV_IDENTIFIER(u, 8, name), __VA_ARGS__)   \
    XSIMD_RVV_WRAPPER##variant(uint16_t, XSIMD_RVV_IDENTIFIER(u, 16, name), __VA_ARGS__) \
    XSIMD_RVV_WRAPPER##variant(uint32_t, XSIMD_RVV_IDENTIFIER(u, 32, name), __VA_ARGS__) \
    XSIMD_RVV_WRAPPER##variant(uint64_t, XSIMD_RVV_IDENTIFIER(u, 64, name), __VA_ARGS__)
#define XSIMD_RVV_OVERLOAD_f(name, variant, ...)                                      \
    XSIMD_RVV_WRAPPER##variant(float, XSIMD_RVV_IDENTIFIER(f, 32, name), __VA_ARGS__) \
    XSIMD_RVV_WRAPPER##variant(double, XSIMD_RVV_IDENTIFIER(f, 64, name), __VA_ARGS__)
#define XSIMD_RVV_OVERLOAD_tail(my_name, variant, ...) \
    XSIMD_RVV_WRAPPER_TAIL##variant(my_name, __VA_ARGS__)

// Use these to create function (actually functor, sorry) wrappers overloaded
// for whichever types are supported.  Being functors means they can't take a
// template argument (until C++14), so if a type can't be deduced then a junk
// value can be passed as the first argument and discarded by using the
// _DROP_1ST variant, instead.
//
// The wrappers use the rvv_reg_t<> types for template accessibility, and
// because some types (eg., vfloat64mf2_t) don't exist and need extra
// abstraction to emulate.
//
// In many cases the intrinsic names are different for signed, unsigned, or
// float variants, the macros OVERLOAD2 and OVERLOAD3 (depending on whether or
// not a float variant exists) take multiple intrinsic names and bring them
// together under a single overloaded identifier where they can be used within
// templates.
//
#define XSIMD_RVV_OVERLOAD2(my_name, name_i, name_u, variant, ...) \
    XSIMD_RVV_OVERLOAD_head(my_name, variant, __VA_ARGS__)         \
    XSIMD_RVV_OVERLOAD_i(name_i, variant, __VA_ARGS__)         \
    XSIMD_RVV_OVERLOAD_u(name_u, variant, __VA_ARGS__)     \
    XSIMD_RVV_OVERLOAD_tail(my_name, variant, __VA_ARGS__)

#define XSIMD_RVV_OVERLOAD3(my_name, name_i, name_u, name_f, variant, ...) \
    XSIMD_RVV_OVERLOAD_head(my_name, variant, __VA_ARGS__)                 \
    XSIMD_RVV_OVERLOAD_i(name_i, variant, __VA_ARGS__)                     \
    XSIMD_RVV_OVERLOAD_u(name_u, variant, __VA_ARGS__)                     \
    XSIMD_RVV_OVERLOAD_f(name_f, variant, __VA_ARGS__)                     \
    XSIMD_RVV_OVERLOAD_tail(my_name, variant, __VA_ARGS__)

#define XSIMD_RVV_OVERLOAD(my_name, name, ...) XSIMD_RVV_OVERLOAD3(my_name, name, name, name, __VA_ARGS__)
#define XSIMD_RVV_OVERLOAD_INTS(my_name, name, ...) XSIMD_RVV_OVERLOAD2(my_name, name, name, __VA_ARGS__)

#define XSIMD_RVV_OVERLOAD_SINTS(my_name, name, variant, ...) \
    XSIMD_RVV_OVERLOAD_head(my_name, variant, __VA_ARGS__)    \
    XSIMD_RVV_OVERLOAD_i(name, variant, __VA_ARGS__)          \
    XSIMD_RVV_OVERLOAD_tail(my_name, variant, __VA_ARGS__)

#define XSIMD_RVV_OVERLOAD_UINTS(my_name, name, variant, ...) \
    XSIMD_RVV_OVERLOAD_head(my_name, variant, __VA_ARGS__)    \
    XSIMD_RVV_OVERLOAD_u(name, variant, __VA_ARGS__)          \
    XSIMD_RVV_OVERLOAD_tail(my_name, variant, __VA_ARGS__)

#define XSIMD_RVV_OVERLOAD_FLOATS(my_name, name, variant, ...) \
    XSIMD_RVV_OVERLOAD_head(my_name, variant, __VA_ARGS__)     \
    XSIMD_RVV_OVERLOAD_f(name, variant, __VA_ARGS__)           \
    XSIMD_RVV_OVERLOAD_tail(my_name, variant, __VA_ARGS__)

// clang-format on

namespace xsimd
{
    template <typename T, class A, T... Values>
    struct batch_constant;

    namespace kernel
    {
        namespace detail
        {
            template <class T>
            using rvv_fix_char_t = types::detail::rvv_fix_char_t<T>;
            template <class T, size_t Width = XSIMD_RVV_BITS>
            using rvv_reg_t = types::detail::rvv_reg_t<T, Width>;
            template <class T, size_t Width = XSIMD_RVV_BITS>
            using rvv_bool_t = types::detail::rvv_bool_t<T, Width>;

            template <size_t>
            struct as_signed_relaxed;
            template <>
            struct as_signed_relaxed<1>
            {
                using type = int8_t;
            };
            template <>
            struct as_signed_relaxed<2>
            {
                using type = int16_t;
            };
            template <>
            struct as_signed_relaxed<4>
            {
                using type = int32_t;
            };
            template <>
            struct as_signed_relaxed<8>
            {
                using type = int64_t;
            };
            template <class T>
            using as_signed_relaxed_t = typename as_signed_relaxed<sizeof(T)>::type;
            template <size_t>
            struct as_unsigned_relaxed;
            template <>
            struct as_unsigned_relaxed<1>
            {
                using type = uint8_t;
            };
            template <>
            struct as_unsigned_relaxed<2>
            {
                using type = uint16_t;
            };
            template <>
            struct as_unsigned_relaxed<4>
            {
                using type = uint32_t;
            };
            template <>
            struct as_unsigned_relaxed<8>
            {
                using type = uint64_t;
            };
            template <class T>
            using as_unsigned_relaxed_t = typename as_unsigned_relaxed<sizeof(T)>::type;
            template <size_t>
            struct as_float_relaxed;
            template <>
            struct as_float_relaxed<1>
            {
                using type = int8_t;
            };
            template <>
            struct as_float_relaxed<2>
            {
                using type = int16_t;
            };
            template <>
            struct as_float_relaxed<4>
            {
                using type = float;
            };
            template <>
            struct as_float_relaxed<8>
            {
                using type = double;
            };
            template <class T>
            using as_float_relaxed_t = typename as_float_relaxed<sizeof(T)>::type;

            template <class T, class U>
            rvv_reg_t<T, U::width> rvvreinterpret(U const& arg) noexcept
            {
                return rvv_reg_t<T, U::width>(arg, types::detail::XSIMD_RVV_BITCAST);
            }
            template <class T, class A, class U>
            rvv_reg_t<T, A::width> rvvreinterpret(batch<U, A> const& arg) noexcept
            {
                typename batch<U, A>::register_type r = arg;
                return rvvreinterpret<T>(r);
            }

            template <class A, class T, class U = as_unsigned_integer_t<T>>
            XSIMD_INLINE batch<U, A> rvv_to_unsigned_batch(batch<T, A> const& arg) noexcept
            {
                return rvvreinterpret<U>(arg.data);
            }

            XSIMD_RVV_OVERLOAD(rvvid,
                               (__riscv_vid_v_u XSIMD_RVV_S XSIMD_RVV_M), _DROP_1ST, uvec(T))

            XSIMD_RVV_OVERLOAD3(rvvmv_splat,
                                (__riscv_vmv_v_x_ XSIMD_RVV_TSM),
                                (__riscv_vmv_v_x_ XSIMD_RVV_TSM),
                                (__riscv_vfmv_v_f_ XSIMD_RVV_TSM), , vec(T))

            XSIMD_RVV_OVERLOAD3(rvvmv_lane0,
                                (__riscv_vmv_x),
                                (__riscv_vmv_x),
                                (__riscv_vfmv_f), _NOVL, T(vec))

            XSIMD_RVV_OVERLOAD(rvvmerge, (__riscv_vmerge), , vec(vec, vec, bvec))
            XSIMD_RVV_OVERLOAD3(rvvmerge_splat,
                                (__riscv_vmerge),
                                (__riscv_vmerge),
                                (__riscv_vfmerge), , vec(vec, T, bvec))

            // count active lanes in a predicate
            XSIMD_RVV_OVERLOAD(rvvcpop, (__riscv_vcpop),
                               , size_t(bvec));

            template <class T, size_t Width>
            XSIMD_INLINE rvv_bool_t<T, Width> pmask8(uint8_t mask) noexcept
            {
                return rvv_bool_t<T, Width>(mask);
            }
            template <class T, size_t Width>
            XSIMD_INLINE rvv_bool_t<T, Width> pmask(uint64_t mask) noexcept
            {
                return rvv_bool_t<T, Width>(mask);
            }

            template <class A, class T, size_t offset = 0, int shift = 0>
            XSIMD_INLINE rvv_reg_t<T, A::width> vindex() noexcept
            {
                auto index = rvvid(T {});
                if (shift < 0)
                    index = __riscv_vsrl(index, -shift, batch<T, A>::size);
                else
                    index = __riscv_vsll(index, shift, batch<T, A>::size);
                return __riscv_vadd(index, T(offset), batch<T, A>::size);
            }

            // enable for signed integers
            template <class T>
            using rvv_enable_signed_int_t = typename std::enable_if<std::is_integral<T>::value && std::is_signed<T>::value, int>::type;

            // enable for unsigned integers
            template <class T>
            using rvv_enable_unsigned_int_t = typename std::enable_if<std::is_integral<T>::value && std::is_unsigned<T>::value, int>::type;

            // enable for floating points
            template <class T>
            using rvv_enable_floating_point_t = typename std::enable_if<std::is_floating_point<T>::value, int>::type;

            // enable for signed integers or floating points
            template <class T>
            using rvv_enable_signed_int_or_floating_point_t = typename std::enable_if<std::is_signed<T>::value, int>::type;

            // enable for all RVE supported types
            template <class T>
            using rvv_enable_all_t = typename std::enable_if<std::is_arithmetic<T>::value, int>::type;
        } // namespace detail

        /********************
         * Scalar to vector *
         ********************/

        namespace detail
        {
            template <class T, size_t Width>
            XSIMD_INLINE detail::rvv_reg_t<T, Width> broadcast(T arg) noexcept
            {
                // A bit of a dance, here, because rvvmv_splat has no other
                // argument from which to deduce type, and T=char is not
                // supported.
                detail::rvv_fix_char_t<T> arg_not_char(arg);
                const auto splat = detail::rvvmv_splat(arg_not_char);
                return detail::rvv_reg_t<T, Width>(splat.get_bytes(), types::detail::XSIMD_RVV_BITCAST);
            }
        }

        // broadcast
        template <class A, class T>
        XSIMD_INLINE batch<T, A> broadcast(T arg, requires_arch<rvv>) noexcept
        {
            return detail::broadcast<T, A::width>(arg);
        }

        /*********
         * Load *
         *********/

        namespace detail
        {
            XSIMD_RVV_OVERLOAD(rvvle, (__riscv_vle XSIMD_RVV_S _v_ XSIMD_RVV_TSM), , vec(T const*))
            XSIMD_RVV_OVERLOAD(rvvse, (__riscv_vse XSIMD_RVV_S _v_ XSIMD_RVV_TSM), , void(T*, vec))
        }

        template <class A, class T, detail::rvv_enable_all_t<T> = 0>
        XSIMD_INLINE batch<T, A> load_aligned(T const* src, convert<T>, requires_arch<rvv>) noexcept
        {
            return detail::rvvle(reinterpret_cast<detail::rvv_fix_char_t<T> const*>(src));
        }

        template <class A, class T, detail::rvv_enable_all_t<T> = 0>
        XSIMD_INLINE batch<T, A> load_unaligned(T const* src, convert<T>, requires_arch<rvv>) noexcept
        {
            return load_aligned<A>(src, convert<T>(), rvv {});
        }

        // load_complex
        namespace detail
        {
            template <class T, size_t W, typename std::enable_if<W >= types::detail::rvv_width_m1, int>::type = 0>
            XSIMD_INLINE rvv_reg_t<T, W * 2> rvvabut(rvv_reg_t<T, W> const& lo, rvv_reg_t<T, W> const& hi) noexcept
            {
                typename rvv_reg_t<T, W * 2>::register_type tmp;
                tmp = __riscv_vset(tmp, 0, lo);
                return __riscv_vset(tmp, 1, hi);
            }

            template <class T, size_t W, typename std::enable_if<W<types::detail::rvv_width_m1, int>::type = 0> XSIMD_INLINE rvv_reg_t<T, W * 2> rvvabut(rvv_reg_t<T, W> const& lo, rvv_reg_t<T, W> const& hi) noexcept
            {
                return __riscv_vslideup(lo, hi, lo.vl, lo.vl * 2);
            }

            XSIMD_RVV_OVERLOAD(rvvget_lo_, (__riscv_vget_ XSIMD_RVV_TSM), _DROP_1ST_CUSTOM_ARGS_NOVL, vec(T, wide_vec), args..., 0)
            XSIMD_RVV_OVERLOAD(rvvget_hi_, (__riscv_vget_ XSIMD_RVV_TSM), _DROP_1ST_CUSTOM_ARGS_NOVL, vec(T, wide_vec), args..., 1)

            template <class T, size_t W, typename std::enable_if<W >= types::detail::rvv_width_m1, int>::type = 0>
            rvv_reg_t<T, W> rvvget_lo(rvv_reg_t<T, W * 2> const& vv) noexcept
            {
                typename rvv_reg_t<T, W>::register_type tmp = rvvget_lo_(T {}, vv);
                return tmp;
            }
            template <class T, size_t W, typename std::enable_if<W >= types::detail::rvv_width_m1, int>::type = 0>
            rvv_reg_t<T, W> rvvget_hi(rvv_reg_t<T, W * 2> const& vv) noexcept
            {
                typename rvv_reg_t<T, W>::register_type tmp = rvvget_hi_(T {}, vv);
                return tmp;
            }
            template <class T, size_t W, typename std::enable_if<W<types::detail::rvv_width_m1, int>::type = 0> rvv_reg_t<T, W> rvvget_lo(rvv_reg_t<T, W * 2> const& vv) noexcept
            {
                typename rvv_reg_t<T, W>::register_type tmp = vv;
                return tmp;
            }
            template <class T, size_t W, typename std::enable_if<W<types::detail::rvv_width_m1, int>::type = 0> rvv_reg_t<T, W> rvvget_hi(rvv_reg_t<T, W * 2> const& vv) noexcept
            {
                return __riscv_vslidedown(vv, vv.vl / 2, vv.vl);
            }

            template <class A, class T, detail::rvv_enable_floating_point_t<T> = 0>
            XSIMD_INLINE batch<std::complex<T>, A> load_complex(batch<T, A> const& lo, batch<T, A> const& hi, requires_arch<rvv>) noexcept
            {
                const auto real_index = vindex<A, as_unsigned_integer_t<T>, 0, 1>();
                const auto imag_index = vindex<A, as_unsigned_integer_t<T>, 1, 1>();
                const auto index = rvvabut<as_unsigned_integer_t<T>, A::width>(real_index, imag_index);
                const auto input = rvvabut<T, A::width>(lo.data, hi.data);
                const rvv_reg_t<T, A::width * 2> result = __riscv_vrgather(input, index, index.vl);

                return { rvvget_lo<T, A::width>(result), rvvget_hi<T, A::width>(result) };
            }
        }

        /*********
         * Store *
         *********/

        template <class A, class T, detail::rvv_enable_all_t<T> = 0>
        XSIMD_INLINE void store_aligned(T* dst, batch<T, A> const& src, requires_arch<rvv>) noexcept
        {
            detail::rvvse(reinterpret_cast<detail::rvv_fix_char_t<T>*>(dst), src);
        }

        template <class A, class T, detail::rvv_enable_all_t<T> = 0>
        XSIMD_INLINE void store_unaligned(T* dst, batch<T, A> const& src, requires_arch<rvv>) noexcept
        {
            store_aligned<A>(dst, src, rvv {});
        }

        /******************
         * scatter/gather *
         ******************/

        namespace detail
        {
            template <class T, class U>
            using rvv_enable_sg_t = typename std::enable_if<(sizeof(T) == sizeof(U) && (sizeof(T) == 4 || sizeof(T) == 8)), int>::type;
            XSIMD_RVV_OVERLOAD(rvvloxei, (__riscv_vloxei XSIMD_RVV_S), , vec(T const*, uvec))
            XSIMD_RVV_OVERLOAD(rvvsoxei, (__riscv_vsoxei XSIMD_RVV_S), , void(T*, uvec, vec))
            XSIMD_RVV_OVERLOAD3(rvvmul_splat,
                                (__riscv_vmul),
                                (__riscv_vmul),
                                (__riscv_vfmul), , vec(vec, T))
        }

        // scatter
        template <class A, class T, class U, detail::rvv_enable_sg_t<T, U> = 0>
        XSIMD_INLINE void scatter(batch<T, A> const& vals, T* dst, batch<U, A> const& index, kernel::requires_arch<rvv>) noexcept
        {
            using UU = as_unsigned_integer_t<U>;
            const auto uindex = detail::rvv_to_unsigned_batch(index);
            auto* base = reinterpret_cast<detail::rvv_fix_char_t<T>*>(dst);
            // or rvvsuxei
            const auto bi = detail::rvvmul_splat(uindex, sizeof(T));
            detail::rvvsoxei(base, bi, vals);
        }

        // gather
        template <class A, class T, class U, detail::rvv_enable_sg_t<T, U> = 0>
        XSIMD_INLINE batch<T, A> gather(batch<T, A> const&, T const* src, batch<U, A> const& index, kernel::requires_arch<rvv>) noexcept
        {
            using UU = as_unsigned_integer_t<U>;
            const auto uindex = detail::rvv_to_unsigned_batch(index);
            auto const* base = reinterpret_cast<detail::rvv_fix_char_t<T> const*>(src);
            // or rvvluxei
            const auto bi = detail::rvvmul_splat(uindex, sizeof(T));
            return detail::rvvloxei(base, bi);
        }

        /**************
         * Arithmetic *
         **************/

        namespace detail
        {
            XSIMD_RVV_OVERLOAD3(rvvadd,
                                (__riscv_vadd),
                                (__riscv_vadd),
                                (__riscv_vfadd), , vec(vec, vec))
            XSIMD_RVV_OVERLOAD2(rvvsadd,
                                (__riscv_vsadd),
                                (__riscv_vsaddu), , vec(vec, vec))
            XSIMD_RVV_OVERLOAD3(rvvsub,
                                (__riscv_vsub),
                                (__riscv_vsub),
                                (__riscv_vfsub), , vec(vec, vec))
            XSIMD_RVV_OVERLOAD2(rvvssub,
                                (__riscv_vssub),
                                (__riscv_vssubu), , vec(vec, vec))
            XSIMD_RVV_OVERLOAD2(rvvaadd,
                                (__riscv_vaadd),
                                (__riscv_vaaddu), , vec(vec, vec))
            XSIMD_RVV_OVERLOAD3(rvvmul,
                                (__riscv_vmul),
                                (__riscv_vmul),
                                (__riscv_vfmul), , vec(vec, vec))
            XSIMD_RVV_OVERLOAD3(rvvdiv,
                                (__riscv_vdiv),
                                (__riscv_vdivu),
                                (__riscv_vfdiv), , vec(vec, vec))
            XSIMD_RVV_OVERLOAD3(rvvmax,
                                (__riscv_vmax),
                                (__riscv_vmaxu),
                                (__riscv_vfmax), , vec(vec, vec))
            XSIMD_RVV_OVERLOAD3(rvvmin,
                                (__riscv_vmin),
                                (__riscv_vminu),
                                (__riscv_vfmin), , vec(vec, vec))
            XSIMD_RVV_OVERLOAD3(rvvneg,
                                (__riscv_vneg),
                                (abort),
                                (__riscv_vfneg), , vec(vec))
            XSIMD_RVV_OVERLOAD_FLOATS(rvvabs,
                                      (__riscv_vfabs), , vec(vec))
            XSIMD_RVV_OVERLOAD3(rvvmacc,
                                (__riscv_vmacc),
                                (__riscv_vmacc),
                                (__riscv_vfmacc), , vec(vec, vec, vec))
            XSIMD_RVV_OVERLOAD3(rvvnmsac,
                                (__riscv_vnmsac),
                                (__riscv_vnmsac),
                                (__riscv_vfnmsac), , vec(vec, vec, vec))
            XSIMD_RVV_OVERLOAD3(rvvmadd,
                                (__riscv_vmadd),
                                (__riscv_vmadd),
                                (__riscv_vfmadd), , vec(vec, vec, vec))
            XSIMD_RVV_OVERLOAD3(rvvnmsub,
                                (__riscv_vnmsub),
                                (__riscv_vnmsub),
                                (__riscv_vfnmsub), , vec(vec, vec, vec))

#define RISCV_VMSXX(XX)                                      \
    XSIMD_RVV_OVERLOAD3(rvvms##XX,                           \
                        (__riscv_vms##XX),                   \
                        (__riscv_vms##XX##u),                \
                        (__riscv_vmf##XX), , bvec(vec, vec)) \
    XSIMD_RVV_OVERLOAD3(rvvms##XX##_splat,                   \
                        (__riscv_vms##XX),                   \
                        (__riscv_vms##XX##u),                \
                        (__riscv_vmf##XX), , bvec(vec, T))
#define __riscv_vmsequ __riscv_vmseq
#define __riscv_vmsneu __riscv_vmsne
            RISCV_VMSXX(eq)
            RISCV_VMSXX(ne)
            RISCV_VMSXX(lt)
            RISCV_VMSXX(le)
            RISCV_VMSXX(gt)
            RISCV_VMSXX(ge)
#undef __riscv_vmsequ
#undef __riscv_vmsneu
#undef RISCV_VMSXX
        } // namespace detail

        // add
        template <class A, class T, detail::rvv_enable_all_t<T> = 0>
        XSIMD_INLINE batch<T, A> add(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<rvv>) noexcept
        {
            return detail::rvvadd(lhs, rhs);
        }

        // sadd
        template <class A, class T, detail::enable_integral_t<T> = 0>
        XSIMD_INLINE batch<T, A> sadd(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<rvv>) noexcept
        {
            return detail::rvvsadd(lhs, rhs);
        }

        // sub
        template <class A, class T, detail::rvv_enable_all_t<T> = 0>
        XSIMD_INLINE batch<T, A> sub(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<rvv>) noexcept
        {
            return detail::rvvsub(lhs, rhs);
        }

        // ssub
        template <class A, class T, detail::enable_integral_t<T> = 0>
        XSIMD_INLINE batch<T, A> ssub(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<rvv>) noexcept
        {
            return detail::rvvssub(lhs, rhs);
        }

        // mul
        template <class A, class T, detail::rvv_enable_all_t<T> = 0>
        XSIMD_INLINE batch<T, A> mul(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<rvv>) noexcept
        {
            return detail::rvvmul(lhs, rhs);
        }

        // div
        template <class A, class T, typename detail::rvv_enable_all_t<T> = 0>
        XSIMD_INLINE batch<T, A> div(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<rvv>) noexcept
        {
            return detail::rvvdiv(lhs, rhs);
        }

        // max
        template <class A, class T, detail::rvv_enable_all_t<T> = 0>
        XSIMD_INLINE batch<T, A> max(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<rvv>) noexcept
        {
            return detail::rvvmax(lhs, rhs);
        }

        // min
        template <class A, class T, detail::rvv_enable_all_t<T> = 0>
        XSIMD_INLINE batch<T, A> min(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<rvv>) noexcept
        {
            return detail::rvvmin(lhs, rhs);
        }

        // neg
        template <class A, class T, detail::rvv_enable_unsigned_int_t<T> = 0>
        XSIMD_INLINE batch<T, A> neg(batch<T, A> const& arg, requires_arch<rvv>) noexcept
        {
            using S = as_signed_integer_t<T>;
            const auto as_signed = detail::rvvreinterpret<S>(arg);
            const auto result = detail::rvvneg(as_signed);
            return detail::rvvreinterpret<T>(result);
        }

        template <class A, class T, detail::rvv_enable_signed_int_or_floating_point_t<T> = 0>
        XSIMD_INLINE batch<T, A> neg(batch<T, A> const& arg, requires_arch<rvv>) noexcept
        {
            return detail::rvvneg(arg);
        }

        // abs
        template <class A, class T, detail::rvv_enable_unsigned_int_t<T> = 0>
        XSIMD_INLINE batch<T, A> abs(batch<T, A> const& arg, requires_arch<rvv>) noexcept
        {
            return arg;
        }

        template <class A, class T, detail::rvv_enable_floating_point_t<T> = 0>
        XSIMD_INLINE batch<T, A> abs(batch<T, A> const& arg, requires_arch<rvv>) noexcept
        {
            return detail::rvvabs(arg);
        }

        // fma: x * y + z
        template <class A, class T, detail::rvv_enable_all_t<T> = 0>
        XSIMD_INLINE batch<T, A> fma(batch<T, A> const& x, batch<T, A> const& y, batch<T, A> const& z, requires_arch<rvv>) noexcept
        {
            // also detail::rvvmadd(x, y, z);
            return detail::rvvmacc(z, x, y);
        }

        // fnma: z - x * y
        template <class A, class T, detail::rvv_enable_all_t<T> = 0>
        XSIMD_INLINE batch<T, A> fnma(batch<T, A> const& x, batch<T, A> const& y, batch<T, A> const& z, requires_arch<rvv>) noexcept
        {
            // also detail::rvvnmsub(x, y, z);
            return detail::rvvnmsac(z, x, y);
        }

        // fms: x * y - z
        template <class A, class T, detail::rvv_enable_all_t<T> = 0>
        XSIMD_INLINE batch<T, A> fms(batch<T, A> const& x, batch<T, A> const& y, batch<T, A> const& z, requires_arch<rvv>) noexcept
        {
            // also vfmsac(z, x, y), but lacking integer version
            // also vfmsub(x, y, z), but lacking integer version
            return -fnma(x, y, z);
        }

        // fnms: - x * y - z
        template <class A, class T, detail::rvv_enable_all_t<T> = 0>
        XSIMD_INLINE batch<T, A> fnms(batch<T, A> const& x, batch<T, A> const& y, batch<T, A> const& z, requires_arch<rvv>) noexcept
        {
            // also vfnmacc(z, x, y), but lacking integer version
            // also vfnmadd(x, y, z), but lacking integer version
            return -fma(z, x, y);
        }

        /**********************
         * Logical operations *
         **********************/

        namespace detail
        {
            XSIMD_RVV_OVERLOAD_INTS(rvvand, (__riscv_vand), , vec(vec, vec))
            XSIMD_RVV_OVERLOAD_INTS(rvvor, (__riscv_vor), , vec(vec, vec))
            XSIMD_RVV_OVERLOAD_INTS(rvvor_splat, (__riscv_vor), , vec(vec, T))
            XSIMD_RVV_OVERLOAD_INTS(rvvxor, (__riscv_vxor), , vec(vec, vec))
            XSIMD_RVV_OVERLOAD_INTS(rvvnot, (__riscv_vnot), , vec(vec))
            XSIMD_RVV_OVERLOAD(rvvmand, (__riscv_vmand_mm_b XSIMD_RVV_S), , bvec(bvec, bvec))
            XSIMD_RVV_OVERLOAD(rvvmor, (__riscv_vmor_mm_b XSIMD_RVV_S), , bvec(bvec, bvec))
            XSIMD_RVV_OVERLOAD(rvvmxor, (__riscv_vmxor_mm_b XSIMD_RVV_S), , bvec(bvec, bvec))
            XSIMD_RVV_OVERLOAD(rvvmandn, (__riscv_vmandn_mm_b XSIMD_RVV_S), , bvec(bvec, bvec))
            XSIMD_RVV_OVERLOAD(rvvmnot, (__riscv_vmnot), , bvec(bvec))
        }

        // bitwise_and
        template <class A, class T, detail::enable_integral_t<T> = 0>
        XSIMD_INLINE batch<T, A> bitwise_and(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<rvv>) noexcept
        {
            return detail::rvvand(lhs, rhs);
        }

        template <class A, class T, detail::rvv_enable_floating_point_t<T> = 0>
        XSIMD_INLINE batch<T, A> bitwise_and(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<rvv>) noexcept
        {
            const auto lhs_bits = detail::rvv_to_unsigned_batch(lhs);
            const auto rhs_bits = detail::rvv_to_unsigned_batch(rhs);
            const auto result_bits = detail::rvvand(lhs_bits, rhs_bits);
            return detail::rvvreinterpret<T>(result_bits);
        }

        template <class A, class T, detail::rvv_enable_all_t<T> = 0>
        XSIMD_INLINE batch_bool<T, A> bitwise_and(batch_bool<T, A> const& lhs, batch_bool<T, A> const& rhs, requires_arch<rvv>) noexcept
        {
            return detail::rvvmand(lhs, rhs);
        }

        // bitwise_andnot
        template <class A, class T, detail::enable_integral_t<T> = 0>
        XSIMD_INLINE batch<T, A> bitwise_andnot(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<rvv>) noexcept
        {
            const auto not_rhs = detail::rvvnot(rhs);
            return detail::rvvand(lhs, not_rhs);
        }

        template <class A, class T, detail::rvv_enable_floating_point_t<T> = 0>
        XSIMD_INLINE batch<T, A> bitwise_andnot(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<rvv>) noexcept
        {
            const auto lhs_bits = detail::rvv_to_unsigned_batch(lhs);
            const auto rhs_bits = detail::rvv_to_unsigned_batch(rhs);
            const auto not_rhs = detail::rvvnot(rhs_bits);
            const auto result_bits = detail::rvvand(lhs_bits, not_rhs);
            return detail::rvvreinterpret<T>(result_bits);
        }

        template <class A, class T, detail::rvv_enable_all_t<T> = 0>
        XSIMD_INLINE batch_bool<T, A> bitwise_andnot(batch_bool<T, A> const& lhs, batch_bool<T, A> const& rhs, requires_arch<rvv>) noexcept
        {
            return detail::rvvmandn(lhs, rhs);
        }

        // bitwise_or
        template <class A, class T, detail::enable_integral_t<T> = 0>
        XSIMD_INLINE batch<T, A> bitwise_or(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<rvv>) noexcept
        {
            return detail::rvvor(lhs, rhs);
        }

        template <class A, class T, detail::rvv_enable_floating_point_t<T> = 0>
        XSIMD_INLINE batch<T, A> bitwise_or(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<rvv>) noexcept
        {
            const auto lhs_bits = detail::rvv_to_unsigned_batch(lhs);
            const auto rhs_bits = detail::rvv_to_unsigned_batch(rhs);
            const auto result_bits = detail::rvvor(lhs_bits, rhs_bits);
            return detail::rvvreinterpret<T>(result_bits);
        }

        template <class A, class T, detail::rvv_enable_all_t<T> = 0>
        XSIMD_INLINE batch_bool<T, A> bitwise_or(batch_bool<T, A> const& lhs, batch_bool<T, A> const& rhs, requires_arch<rvv>) noexcept
        {
            return detail::rvvmor(lhs, rhs);
        }

        // bitwise_xor
        template <class A, class T, detail::enable_integral_t<T> = 0>
        XSIMD_INLINE batch<T, A> bitwise_xor(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<rvv>) noexcept
        {
            return detail::rvvxor(lhs, rhs);
        }

        template <class A, class T, detail::rvv_enable_floating_point_t<T> = 0>
        XSIMD_INLINE batch<T, A> bitwise_xor(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<rvv>) noexcept
        {
            const auto lhs_bits = detail::rvv_to_unsigned_batch(lhs);
            const auto rhs_bits = detail::rvv_to_unsigned_batch(rhs);
            const auto result_bits = detail::rvvxor(lhs_bits, rhs_bits);
            return detail::rvvreinterpret<T>(result_bits);
        }

        template <class A, class T, detail::rvv_enable_all_t<T> = 0>
        XSIMD_INLINE batch_bool<T, A> bitwise_xor(batch_bool<T, A> const& lhs, batch_bool<T, A> const& rhs, requires_arch<rvv>) noexcept
        {
            return detail::rvvmxor(lhs, rhs);
        }

        // bitwise_not
        template <class A, class T, detail::enable_integral_t<T> = 0>
        XSIMD_INLINE batch<T, A> bitwise_not(batch<T, A> const& arg, requires_arch<rvv>) noexcept
        {
            return detail::rvvnot(arg);
        }

        template <class A, class T, detail::rvv_enable_floating_point_t<T> = 0>
        XSIMD_INLINE batch<T, A> bitwise_not(batch<T, A> const& arg, requires_arch<rvv>) noexcept
        {
            const auto arg_bits = detail::rvv_to_unsigned_batch(arg);
            const auto result_bits = detail::rvvnot(arg_bits);
            return detail::rvvreinterpret<T>(result_bits);
        }

        template <class A, class T, detail::rvv_enable_all_t<T> = 0>
        XSIMD_INLINE batch_bool<T, A> bitwise_not(batch_bool<T, A> const& arg, requires_arch<rvv>) noexcept
        {
            return detail::rvvmnot(arg);
        }

        /**********
         * Shifts *
         **********/

        namespace detail
        {
            XSIMD_RVV_OVERLOAD_INTS(rvvsll_splat, (__riscv_vsll), , vec(vec, size_t))
            XSIMD_RVV_OVERLOAD_INTS(rvvsll, (__riscv_vsll), , vec(vec, uvec))
            XSIMD_RVV_OVERLOAD2(rvvsr_splat,
                                (__riscv_vsra),
                                (__riscv_vsrl), , vec(vec, size_t))
            XSIMD_RVV_OVERLOAD2(rvvsr,
                                (__riscv_vsra),
                                (__riscv_vsrl), , vec(vec, uvec))
        } // namespace detail

        // bitwise_lshift
        template <class A, class T, detail::enable_integral_t<T> = 0>
        XSIMD_INLINE batch<T, A> bitwise_lshift(batch<T, A> const& arg, int n, requires_arch<rvv>) noexcept
        {
            constexpr size_t size = sizeof(typename batch<T, A>::value_type) * 8;
            assert(0 <= n && static_cast<size_t>(n) < size && "index in bounds");
            return detail::rvvsll_splat(arg, n);
        }

        template <class A, class T, detail::enable_integral_t<T> = 0>
        XSIMD_INLINE batch<T, A> bitwise_lshift(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<rvv>) noexcept
        {
            return detail::rvvsll(lhs, detail::rvv_to_unsigned_batch<A, T>(rhs));
        }

        // bitwise_rshift
        template <class A, class T, detail::enable_integral_t<T> = 0>
        XSIMD_INLINE batch<T, A> bitwise_rshift(batch<T, A> const& arg, int n, requires_arch<rvv>) noexcept
        {
            constexpr size_t size = sizeof(typename batch<T, A>::value_type) * 8;
            assert(0 <= n && static_cast<size_t>(n) < size && "index in bounds");
            return detail::rvvsr_splat(arg, n);
        }

        template <class A, class T, detail::enable_integral_t<T> = 0>
        XSIMD_INLINE batch<T, A> bitwise_rshift(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<rvv>) noexcept
        {
            return detail::rvvsr(lhs, detail::rvv_to_unsigned_batch<A, T>(rhs));
        }

        /**************
         * Reductions *
         **************/

        namespace detail
        {
            XSIMD_RVV_OVERLOAD3(rvvredsum,
                                (__riscv_vredsum),
                                (__riscv_vredsum),
                                (__riscv_vfredosum), // or __riscv_vfredusum
                                , scalar_vec(vec, scalar_vec))
            XSIMD_RVV_OVERLOAD3(rvvredmax,
                                (__riscv_vredmax),
                                (__riscv_vredmaxu),
                                (__riscv_vfredmax), , scalar_vec(vec, scalar_vec))
            XSIMD_RVV_OVERLOAD3(rvvredmin,
                                (__riscv_vredmin),
                                (__riscv_vredminu),
                                (__riscv_vfredmin), , scalar_vec(vec, scalar_vec))
            XSIMD_RVV_OVERLOAD3(rvvslide1up,
                                (__riscv_vslide1up),
                                (__riscv_vslide1up),
                                (__riscv_vfslide1up), , vec(vec, vec))
            XSIMD_RVV_OVERLOAD3(rvvslide1down,
                                (__riscv_vslide1down),
                                (__riscv_vslide1down),
                                (__riscv_vfslide1down), , vec(vec, T))

            template <class A, class T>
            XSIMD_INLINE T reduce_scalar(rvv_reg_t<T, types::detail::rvv_width_m1> const& arg)
            {
                return detail::rvvmv_lane0(rvv_reg_t<T, A::width>(arg.get_bytes(), types::detail::XSIMD_RVV_BITCAST));
            }
        }
        // reduce_add
        template <class A, class T, class V = typename batch<T, A>::value_type, detail::rvv_enable_all_t<T> = 0>
        XSIMD_INLINE V reduce_add(batch<T, A> const& arg, requires_arch<rvv>) noexcept
        {
            const auto zero = detail::broadcast<T, types::detail::rvv_width_m1>(T(0));
            const auto r = detail::rvvredsum(arg, zero);
            return detail::reduce_scalar<A, T>(r);
        }

        // reduce_max
        template <class A, class T, detail::rvv_enable_all_t<T> = 0>
        XSIMD_INLINE T reduce_max(batch<T, A> const& arg, requires_arch<rvv>) noexcept
        {
            const auto lowest = detail::broadcast<T, types::detail::rvv_width_m1>(std::numeric_limits<T>::lowest());
            const auto r = detail::rvvredmax(arg, lowest);
            return detail::reduce_scalar<A, T>(r);
        }

        // reduce_min
        template <class A, class T, detail::rvv_enable_all_t<T> = 0>
        XSIMD_INLINE T reduce_min(batch<T, A> const& arg, requires_arch<rvv>) noexcept
        {
            const auto max = detail::broadcast<T, types::detail::rvv_width_m1>(std::numeric_limits<T>::max());
            const auto r = detail::rvvredmin(arg, max);
            return detail::reduce_scalar<A, T>(r);
        }

        // haddp
        template <class A, class T, detail::rvv_enable_floating_point_t<T> = 0>
        XSIMD_INLINE batch<T, A> haddp(const batch<T, A>* row, requires_arch<rvv>) noexcept
        {
            constexpr std::size_t size = batch<T, A>::size;
            T sums[size];
#pragma unroll size
            for (std::size_t i = 0; i < size; ++i)
            {
                sums[i] = reduce_add(row[i], rvv {});
            }
            return load_aligned<A>(sums, convert<T>(), rvv {});
        }

        /***************
         * Comparisons *
         ***************/

        // eq
        template <class A, class T, detail::rvv_enable_all_t<T> = 0>
        XSIMD_INLINE batch_bool<T, A> eq(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<rvv>) noexcept
        {
            return detail::rvvmseq(lhs, rhs);
        }

        template <class A, class T, detail::rvv_enable_all_t<T> = 0>
        XSIMD_INLINE batch_bool<T, A> eq(batch_bool<T, A> const& lhs, batch_bool<T, A> const& rhs, requires_arch<rvv>) noexcept
        {
            const auto neq_result = detail::rvvmxor(lhs, rhs);
            return detail::rvvmnot(neq_result);
        }

        // neq
        template <class A, class T, detail::rvv_enable_all_t<T> = 0>
        XSIMD_INLINE batch_bool<T, A> neq(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<rvv>) noexcept
        {
            return detail::rvvmsne(lhs, rhs);
        }

        template <class A, class T, detail::rvv_enable_all_t<T> = 0>
        XSIMD_INLINE batch_bool<T, A> neq(batch_bool<T, A> const& lhs, batch_bool<T, A> const& rhs, requires_arch<rvv>) noexcept
        {
            return detail::rvvmxor(lhs, rhs);
        }

        // lt
        template <class A, class T, detail::rvv_enable_all_t<T> = 0>
        XSIMD_INLINE batch_bool<T, A> lt(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<rvv>) noexcept
        {
            return detail::rvvmslt(lhs, rhs);
        }

        // le
        template <class A, class T, detail::rvv_enable_all_t<T> = 0>
        XSIMD_INLINE batch_bool<T, A> le(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<rvv>) noexcept
        {
            return detail::rvvmsle(lhs, rhs);
        }

        // gt
        template <class A, class T, detail::rvv_enable_all_t<T> = 0>
        XSIMD_INLINE batch_bool<T, A> gt(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<rvv>) noexcept
        {
            return detail::rvvmsgt(lhs, rhs);
        }

        // ge
        template <class A, class T, detail::rvv_enable_all_t<T> = 0>
        XSIMD_INLINE batch_bool<T, A> ge(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<rvv>) noexcept
        {
            return detail::rvvmsge(lhs, rhs);
        }

        /*************
         * Selection *
         *************/
        namespace detail
        {
            XSIMD_RVV_OVERLOAD(rvvcompress, (__riscv_vcompress), , vec(vec, bvec))
        }
        // compress
        template <class A, class T>
        XSIMD_INLINE batch<T, A> compress(batch<T, A> const& x, batch_bool<T, A> const& mask, requires_arch<rvv>) noexcept
        {
            return detail::rvvcompress(x, mask);
        }

        /***************
         * Permutation *
         ***************/
        namespace detail
        {
            XSIMD_RVV_OVERLOAD(rvvrgather, (__riscv_vrgather), , vec(vec, uvec))
            XSIMD_RVV_OVERLOAD(rvvslideup, (__riscv_vslideup), , vec(vec, vec, size_t))
            XSIMD_RVV_OVERLOAD(rvvslidedown, (__riscv_vslidedown), , vec(vec, size_t))
        }

        // swizzle
        template <class A, class T, class I, I... idx>
        XSIMD_INLINE batch<T, A> swizzle(batch<T, A> const& arg, batch_constant<I, A, idx...>, requires_arch<rvv>) noexcept
        {
            static_assert(batch<T, A>::size == sizeof...(idx), "invalid swizzle indices");
            const batch<I, A> indices { idx... };
            return detail::rvvrgather(arg, indices);
        }

        template <class A, class T, class I, I... idx>
        XSIMD_INLINE batch<std::complex<T>, A> swizzle(batch<std::complex<T>, A> const& self,
                                                       batch_constant<I, A, idx...>,
                                                       requires_arch<rvv>) noexcept
        {
            const auto real = swizzle(self.real(), batch_constant<I, A, idx...> {}, rvv {});
            const auto imag = swizzle(self.imag(), batch_constant<I, A, idx...> {}, rvv {});
            return batch<std::complex<T>>(real, imag);
        }

        /*************
         * Selection *
         *************/

        // extract_pair

        template <class A, class T, detail::rvv_enable_all_t<T> = 0>
        XSIMD_INLINE batch<T, A> extract_pair(batch<T, A> const& lhs, batch<T, A> const& rhs, size_t n, requires_arch<rvv>) noexcept
        {
            const auto tmp = detail::rvvslidedown(rhs, n);
            return detail::rvvslideup(tmp, lhs, lhs.size - n);
        }

        // select
        template <class A, class T, detail::rvv_enable_all_t<T> = 0>
        XSIMD_INLINE batch<T, A> select(batch_bool<T, A> const& cond, batch<T, A> const& a, batch<T, A> const& b, requires_arch<rvv>) noexcept
        {
            return detail::rvvmerge(b, a, cond);
        }

        template <class A, class T, bool... b>
        XSIMD_INLINE batch<T, A> select(batch_bool_constant<T, A, b...> const&, batch<T, A> const& true_br, batch<T, A> const& false_br, requires_arch<rvv>) noexcept
        {
            return select(batch_bool<T, A> { b... }, true_br, false_br, rvv {});
        }

        // zip_lo
        template <class A, class T, detail::rvv_enable_all_t<T> = 0>
        XSIMD_INLINE batch<T, A> zip_lo(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<rvv>) noexcept
        {
            const auto index = detail::vindex<A, as_unsigned_integer_t<T>, 0, -1>();
            const auto mask = detail::pmask8<T, A::width>(0xaa);
            return detail::rvvmerge(detail::rvvrgather(lhs, index),
                                    detail::rvvrgather(rhs, index),
                                    mask);
        }

        // zip_hi
        template <class A, class T, detail::rvv_enable_all_t<T> = 0>
        XSIMD_INLINE batch<T, A> zip_hi(batch<T, A> const& lhs, batch<T, A> const& rhs, requires_arch<rvv>) noexcept
        {
            const auto index = detail::vindex<A, as_unsigned_integer_t<T>, batch<T, A>::size / 2, -1>();
            const auto mask = detail::pmask8<T, A::width>(0xaa);
            return detail::rvvmerge(detail::rvvrgather(lhs, index),
                                    detail::rvvrgather(rhs, index),
                                    mask);
        }

        // store_complex
        template <class A, class T, detail::rvv_enable_floating_point_t<T> = 0>
        XSIMD_INLINE void store_complex_aligned(std::complex<T>* dst, batch<std::complex<T>, A> const& src, requires_arch<rvv>) noexcept
        {
            const auto lo = zip_lo(src.real(), src.imag());
            const auto hi = zip_hi(src.real(), src.imag());
            T* buf = reinterpret_cast<T*>(dst);
            store_aligned(buf, lo, rvv {});
            store_aligned(buf + lo.size, hi, rvv {});
        }

        template <class A, class T, detail::rvv_enable_floating_point_t<T> = 0>
        XSIMD_INLINE void store_complex_unaligned(std::complex<T>* dst, batch<std::complex<T>, A> const& src, requires_arch<rvv>) noexcept
        {
            store_complex_aligned(dst, src, rvv {});
        }

        /*****************************
         * Floating-point arithmetic *
         *****************************/

        namespace detail
        {
            XSIMD_RVV_OVERLOAD_FLOATS(rvvfsqrt, (__riscv_vfsqrt), , vec(vec))
            XSIMD_RVV_OVERLOAD_FLOATS(rvvfrec7, (__riscv_vfrec7), , vec(vec))
            XSIMD_RVV_OVERLOAD_FLOATS(rvvfrsqrt7, (__riscv_vfrsqrt7), , vec(vec))
        }

        // rsqrt
        template <class A, class T, detail::rvv_enable_floating_point_t<T> = 0>
        XSIMD_INLINE batch<T, A> rsqrt(batch<T, A> const& arg, requires_arch<rvv>) noexcept
        {
            auto approx = detail::rvvfrsqrt7(arg);
            approx = approx * (1.5 - (0.5 * arg * approx * approx));
            return approx;
        }

        // sqrt
        template <class A, class T, detail::rvv_enable_floating_point_t<T> = 0>
        XSIMD_INLINE batch<T, A> sqrt(batch<T, A> const& arg, requires_arch<rvv>) noexcept
        {
            return detail::rvvfsqrt(arg);
        }

        // reciprocal
        template <class A, class T, detail::rvv_enable_floating_point_t<T> = 0>
        XSIMD_INLINE batch<T, A> reciprocal(const batch<T, A>& arg, requires_arch<rvv>) noexcept
        {
            return detail::rvvfrec7(arg);
        }

        /******************************
         * Floating-point conversions *
         ******************************/

        // fast_cast
        namespace detail
        {
            XSIMD_RVV_OVERLOAD2(rvvfcvt_rtz, // truncating conversion, like C.
                                (__riscv_vfcvt_rtz_x),
                                (__riscv_vfcvt_rtz_xu), _DROP_1ST, vec(T, fvec))
            XSIMD_RVV_OVERLOAD2(rvvfcvt_rne, // round to nearest, ties to even
                                (__riscv_vfcvt_x),
                                (__riscv_vfcvt_xu), _DROP_1ST_CUSTOM_ARGS, vec(T, fvec), args..., __RISCV_FRM_RNE)
            XSIMD_RVV_OVERLOAD2(rvvfcvt_rmm, // round to nearest, ties to max magnitude
                                (__riscv_vfcvt_x),
                                (__riscv_vfcvt_xu), _DROP_1ST_CUSTOM_ARGS, vec(T, fvec), args..., __RISCV_FRM_RMM)
            XSIMD_RVV_OVERLOAD2(rvvfcvt, // round to current rounding mode.
                                (__riscv_vfcvt_x),
                                (__riscv_vfcvt_xu), _DROP_1ST, vec(T, fvec))
            XSIMD_RVV_OVERLOAD_INTS(rvvfcvt_f, (__riscv_vfcvt_f), , fvec(vec))

            template <class T, class U>
            using rvv_enable_ftoi_t = typename std::enable_if<(sizeof(T) == sizeof(U) && std::is_floating_point<T>::value && !std::is_floating_point<U>::value), int>::type;
            template <class T, class U>
            using rvv_enable_itof_t = typename std::enable_if<(sizeof(T) == sizeof(U) && !std::is_floating_point<T>::value && std::is_floating_point<U>::value), int>::type;

            template <class A, class T, class U, rvv_enable_ftoi_t<T, U> = 0>
            XSIMD_INLINE batch<U, A> fast_cast(batch<T, A> const& arg, batch<U, A> const&, requires_arch<rvv>) noexcept
            {
                return rvvfcvt_rtz(U {}, arg);
            }
            template <class A, class T, class U, rvv_enable_itof_t<T, U> = 0>
            XSIMD_INLINE batch<U, A> fast_cast(batch<T, A> const& arg, batch<U, A> const&, requires_arch<rvv>) noexcept
            {
                return rvvfcvt_f(arg);
            }
        }

        /*********
         * Miscs *
         *********/

        // set
        template <class A, class T, class... Args>
        XSIMD_INLINE batch<T, A> set(batch<T, A> const&, requires_arch<rvv>, Args... args) noexcept
        {
            const std::array<T, batch<T, A>::size> tmp { args... };
            return load_unaligned<A>(tmp.data(), convert<T>(), rvv {});
        }

        template <class A, class T, class... Args>
        XSIMD_INLINE batch<std::complex<T>, A> set(batch<std::complex<T>, A> const&, requires_arch<rvv>,
                                                   Args... args_complex) noexcept
        {
            return batch<std::complex<T>>(set(batch<T, rvv> {}, rvv {}, args_complex.real()...),
                                          set(batch<T, rvv> {}, rvv {}, args_complex.imag()...));
        }

        template <class A, class T, class... Args>
        XSIMD_INLINE batch_bool<T, A> set(batch_bool<T, A> const&, requires_arch<rvv>, Args... args) noexcept
        {
            using U = as_unsigned_integer_t<T>;
            const auto values = set(batch<U, rvv> {}, rvv {}, static_cast<U>(args)...);
            const auto zero = broadcast<A>(U(0), rvv {});
            detail::rvv_bool_t<T> result = detail::rvvmsne(values, zero);
            return result;
        }

        // insert
        template <class A, class T, size_t I, detail::rvv_enable_all_t<T> = 0>
        XSIMD_INLINE batch<T, A> insert(batch<T, A> const& arg, T val, index<I>, requires_arch<rvv>) noexcept
        {
            const auto mask = detail::pmask<T, A::width>(uint64_t(1) << I);
            return detail::rvvmerge_splat(arg, val, mask);
        }

        // get
        template <class A, class T, detail::rvv_enable_all_t<T> = 0>
        XSIMD_INLINE T get(batch<T, A> const& arg, size_t i, requires_arch<rvv>) noexcept
        {
            const auto tmp = detail::rvvslidedown(arg, i);
            return detail::rvvmv_lane0(tmp);
        }

        template <class A, class T, detail::rvv_enable_all_t<T> = 0>
        XSIMD_INLINE std::complex<T> get(batch<std::complex<T>, A> const& arg, size_t i, requires_arch<rvv>) noexcept
        {
            const auto tmpr = detail::rvvslidedown(arg.real(), i);
            const auto tmpi = detail::rvvslidedown(arg.imag(), i);
            return std::complex<T> { detail::rvvmv_lane0(tmpr), detail::rvvmv_lane0(tmpi) };
        }

        // all
        template <class A, class T, detail::rvv_enable_all_t<T> = 0>
        XSIMD_INLINE bool all(batch_bool<T, A> const& arg, requires_arch<rvv>) noexcept
        {
            return detail::rvvcpop(arg) == batch_bool<T, A>::size;
        }

        // any
        template <class A, class T, detail::rvv_enable_all_t<T> = 0>
        XSIMD_INLINE bool any(batch_bool<T, A> const& arg, requires_arch<rvv>) noexcept
        {
            return detail::rvvcpop(arg) > 0;
        }

        // bitwise_cast
        template <class A, class T, class R, detail::rvv_enable_all_t<T> = 0, detail::rvv_enable_all_t<R> = 0>
        XSIMD_INLINE batch<R, A> bitwise_cast(batch<T, A> const& arg, batch<R, A> const&, requires_arch<rvv>) noexcept
        {
            return detail::rvv_reg_t<R, A::width>(arg.data.get_bytes(), types::detail::XSIMD_RVV_BITCAST);
        }

        // batch_bool_cast
        template <class A, class T_out, class T_in, detail::rvv_enable_all_t<T_in> = 0>
        XSIMD_INLINE batch_bool<T_out, A> batch_bool_cast(batch_bool<T_in, A> const& arg, batch_bool<T_out, A> const&, requires_arch<rvv>) noexcept
        {
            using intermediate_t = typename detail::rvv_bool_t<T_out>;
            return intermediate_t(arg.data);
        }

        // from_bool
        template <class A, class T, detail::rvv_enable_all_t<T> = 0>
        XSIMD_INLINE batch<T, A> from_bool(batch_bool<T, A> const& arg, requires_arch<rvv>) noexcept
        {
            const auto zero = broadcast<A>(T(0), rvv {});
            return detail::rvvmerge_splat(zero, T(1), arg);
        }

        namespace detail
        {
            template <size_t Width>
            XSIMD_INLINE vuint8m1_t rvvslidedownbytes(vuint8m1_t arg, size_t i)
            {
                return __riscv_vslidedown(arg, i, types::detail::rvv_width_m1 / 8);
            }
            template <>
            XSIMD_INLINE vuint8m1_t rvvslidedownbytes<types::detail::rvv_width_mf2>(vuint8m1_t arg, size_t i)
            {
                const auto bytes = __riscv_vlmul_trunc_u8mf2(arg);
                const auto result = __riscv_vslidedown(bytes, i, types::detail::rvv_width_mf2 / 8);
                return __riscv_vlmul_ext_u8m1(result);
            }
            template <>
            XSIMD_INLINE vuint8m1_t rvvslidedownbytes<types::detail::rvv_width_mf4>(vuint8m1_t arg, size_t i)
            {
                const auto bytes = __riscv_vlmul_trunc_u8mf4(arg);
                const auto result = __riscv_vslidedown(bytes, i, types::detail::rvv_width_mf4 / 8);
                return __riscv_vlmul_ext_u8m1(result);
            }
            template <>
            XSIMD_INLINE vuint8m1_t rvvslidedownbytes<types::detail::rvv_width_mf8>(vuint8m1_t arg, size_t i)
            {
                const auto bytes = __riscv_vlmul_trunc_u8mf8(arg);
                const auto result = __riscv_vslidedown(bytes, i, types::detail::rvv_width_mf8 / 8);
                return __riscv_vlmul_ext_u8m1(result);
            }
        }

        // slide_left
        template <size_t N, class A, class T, detail::rvv_enable_all_t<T> = 0>
        XSIMD_INLINE batch<T, A> slide_left(batch<T, A> const& arg, requires_arch<rvv>) noexcept
        {
            const auto zero = broadcast<A>(uint8_t(0), rvv {});
            const auto bytes = arg.data.get_bytes();
            return detail::rvvreinterpret<T>(detail::rvvslideup(zero, bytes, N));
        }

        // slide_right
        template <size_t N, class A, class T, detail::rvv_enable_all_t<T> = 0>
        XSIMD_INLINE batch<T, A> slide_right(batch<T, A> const& arg, requires_arch<rvv>) noexcept
        {
            using reg_t = detail::rvv_reg_t<T, A::width>;
            const auto bytes = arg.data.get_bytes();
            return reg_t(detail::rvvslidedownbytes<A::width>(bytes, N), types::detail::XSIMD_RVV_BITCAST);
        }

        // isnan
        template <class A, class T, detail::rvv_enable_floating_point_t<T> = 0>
        XSIMD_INLINE batch_bool<T, A> isnan(batch<T, A> const& arg, requires_arch<rvv>) noexcept
        {
            return !(arg == arg);
        }

        namespace detail
        {
            template <class T>
            using rvv_as_signed_integer_t = as_signed_integer_t<as_unsigned_integer_t<T>>;

            template <class A, class T, class U = rvv_as_signed_integer_t<T>>
            XSIMD_INLINE batch<U, A> rvvfcvt_default(batch<T, A> const& arg) noexcept
            {
                return rvvfcvt_rne(U {}, arg);
            }

            template <class A, class T, class U = rvv_as_signed_integer_t<T>>
            XSIMD_INLINE batch<U, A> rvvfcvt_afz(batch<T, A> const& arg) noexcept
            {
                return rvvfcvt_rmm(U {}, arg);
            }
        }

        // nearbyint_as_int
        template <class A, class T, class U = detail::rvv_as_signed_integer_t<T>>
        XSIMD_INLINE batch<U, A> nearbyint_as_int(batch<T, A> const& arg, requires_arch<rvv>) noexcept
        {
            // Reference rounds ties to nearest even
            return detail::rvvfcvt_default(arg);
        }

        // round
        template <class A, class T, detail::rvv_enable_floating_point_t<T> = 0>
        XSIMD_INLINE batch<T, A> round(batch<T, A> const& arg, requires_arch<rvv>) noexcept
        {
            // Round ties away from zero.
            const auto mask = abs(arg) < constants::maxflint<batch<T, A>>();
            return select(mask, to_float(detail::rvvfcvt_afz(arg)), arg, rvv {});
        }

        // nearbyint
        template <class A, class T, detail::rvv_enable_floating_point_t<T> = 0>
        XSIMD_INLINE batch<T, A> nearbyint(batch<T, A> const& arg, requires_arch<rvv>) noexcept
        {
            // Round according to current rounding mode.
            const auto mask = abs(arg) < constants::maxflint<batch<T, A>>();
            return select(mask, to_float(detail::rvvfcvt_default(arg)), arg, rvv {});
        }
    } // namespace kernel
} // namespace xsimd

#endif
