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

#ifndef XSIMD_EMULATED_REGISTER_HPP
#define XSIMD_EMULATED_REGISTER_HPP

#include "./xsimd_generic_arch.hpp"
#include "./xsimd_register.hpp"

namespace xsimd
{
    /**
     * @ingroup architectures
     *
     * emulated instructions
     */
    template <size_t N>
    struct emulated : generic
    {
        static constexpr bool supported() noexcept { return true; }
        static constexpr bool available() noexcept { return true; }
        static constexpr bool requires_alignment() noexcept { return false; }
        static constexpr std::size_t alignment() noexcept { return 8; }
        static constexpr char const* name() noexcept { return "emulated"; }
    };

    namespace types
    {
        template <size_t N>
        struct simd_emulated_bool_register
        {
            using register_type = std::array<bool, N>;
            register_type data;
            simd_emulated_bool_register() = default;
            simd_emulated_bool_register(register_type r) { data = r; }
            operator register_type() const noexcept { return data; }
        };
        template <typename T, size_t N>
        struct get_bool_simd_register<T, emulated<N>>
        {
            using type = simd_emulated_bool_register<N / (8 * sizeof(T))>;
        };

        template <typename T, size_t N>
        struct simd_register<T, emulated<N>>
        {
            static_assert(N % (8 * sizeof(T)) == 0, "bit width must be a multiple of scalar width");
            using register_type = std::array<T, N / (8 * sizeof(T))>;
            register_type data;
            XSIMD_INLINE operator register_type() const noexcept
            {
                return data;
            }
        };
        template <typename T, size_t N>
        struct has_simd_register<T, emulated<N>> : std::is_scalar<T>
        {
        };
        template <typename T, size_t N>
        struct has_simd_register<std::complex<T>, emulated<N>> : std::true_type
        {
        };
#ifdef XSIMD_ENABLE_XTL_COMPLEX
        template <typename T, bool i3ec, size_t N>
        struct has_simd_register<xtl::complex<T, T, i3ec>, emulated<N>> : std::true_type
        {
        };
#endif
    }
}

#endif
