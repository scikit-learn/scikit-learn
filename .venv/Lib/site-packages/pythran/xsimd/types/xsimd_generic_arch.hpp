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

#ifndef XSIMD_GENERIC_ARCH_HPP
#define XSIMD_GENERIC_ARCH_HPP

#include "../config/xsimd_config.hpp"

/**
 * @defgroup architectures Architecture description
 * */
namespace xsimd
{
    /**
     * @ingroup architectures
     *
     * Base class for all architectures.
     */
    struct generic
    {
        /// Whether this architecture is supported at compile-time.
        static constexpr bool supported() noexcept { return true; }
        /// Whether this architecture is available at run-time.
        static constexpr bool available() noexcept { return true; }
        /// If this architectures supports aligned memory accesses, the required
        /// alignment.
        static constexpr std::size_t alignment() noexcept { return 0; }
        /// Whether this architecture requires aligned memory access.
        static constexpr bool requires_alignment() noexcept { return false; }
        /// Name of the architecture.
        static constexpr char const* name() noexcept { return "generic"; }
    };

    struct unsupported
    {
    };
}

#endif
