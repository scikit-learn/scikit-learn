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

#ifndef XSIMD_ALIGNMENT_HPP
#define XSIMD_ALIGNMENT_HPP

#include "../types/xsimd_utils.hpp"
#include "xsimd_aligned_allocator.hpp"

namespace xsimd
{
    /**
     * @struct aligned_mode
     * @brief tag for load and store of aligned memory.
     */
    struct aligned_mode
    {
    };

    /**
     * @struct unaligned_mode
     * @brief tag for load and store of unaligned memory.
     */
    struct unaligned_mode
    {
    };

    /***********************
     * Allocator alignment *
     ***********************/

    template <class A>
    struct allocator_alignment
    {
        using type = unaligned_mode;
    };

    template <class T, size_t N>
    struct allocator_alignment<aligned_allocator<T, N>>
    {
        using type = aligned_mode;
    };

    template <class A>
    using allocator_alignment_t = typename allocator_alignment<A>::type;

    /***********************
     * container alignment *
     ***********************/

    template <class C, class = void>
    struct container_alignment
    {
        using type = unaligned_mode;
    };

    template <class C>
    struct container_alignment<C, detail::void_t<typename C::allocator_type>>
    {
        using type = allocator_alignment_t<typename C::allocator_type>;
    };

    template <class C>
    using container_alignment_t = typename container_alignment<C>::type;

    /*********************
     * alignment checker *
     *********************/

    /**
     * Checks whether pointer \c ptr is aligned according the alignment
     * requirements of \c Arch.
     * @return true if the alignment requirements are met
     */
    template <class Arch = default_arch>
    XSIMD_INLINE bool is_aligned(void const* ptr)
    {
        return (reinterpret_cast<uintptr_t>(ptr) % static_cast<uintptr_t>(Arch::alignment())) == 0;
    }

}

#endif
