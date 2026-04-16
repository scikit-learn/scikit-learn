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

#ifndef XSIMD_REGISTER_HPP
#define XSIMD_REGISTER_HPP

#include <type_traits>

namespace xsimd
{
    namespace types
    {
        template <class T, class A>
        struct has_simd_register : std::false_type
        {
        };

        template <class T, class Arch>
        struct simd_register
        {
            struct register_type
            {
            };
        };

#define XSIMD_DECLARE_SIMD_REGISTER(SCALAR_TYPE, ISA, VECTOR_TYPE) \
    template <>                                                    \
    struct simd_register<SCALAR_TYPE, ISA>                         \
    {                                                              \
        using register_type = VECTOR_TYPE;                         \
        register_type data;                                        \
        XSIMD_INLINE operator register_type() const noexcept       \
        {                                                          \
            return data;                                           \
        }                                                          \
    };                                                             \
    template <>                                                    \
    struct has_simd_register<SCALAR_TYPE, ISA> : std::true_type    \
    {                                                              \
    }

#define XSIMD_DECLARE_INVALID_SIMD_REGISTER(SCALAR_TYPE, ISA)    \
    template <>                                                  \
    struct has_simd_register<SCALAR_TYPE, ISA> : std::false_type \
    {                                                            \
    }

#define XSIMD_DECLARE_SIMD_REGISTER_ALIAS(ISA, ISA_BASE)                          \
    template <class T>                                                            \
    struct simd_register<T, ISA> : simd_register<T, ISA_BASE>                     \
    {                                                                             \
        using register_type = typename simd_register<T, ISA_BASE>::register_type; \
        simd_register(register_type reg) noexcept                                 \
            : simd_register<T, ISA_BASE> { reg }                                  \
        {                                                                         \
        }                                                                         \
        simd_register() = default;                                                \
    };                                                                            \
    template <class T>                                                            \
    struct has_simd_register<T, ISA> : has_simd_register<T, ISA_BASE>             \
    {                                                                             \
    }

        template <class T, class Arch>
        struct get_bool_simd_register
        {
            using type = simd_register<T, Arch>;
        };

        template <class T, class Arch>
        using get_bool_simd_register_t = typename get_bool_simd_register<T, Arch>::type;
    }

    namespace kernel
    {
        template <class A>
        // makes requires_arch equal to A const&, using type_traits functions
        using requires_arch = typename std::add_lvalue_reference<typename std::add_const<A>::type>::type;
        template <class T>
        struct convert
        {
        };
    }
}

#endif
