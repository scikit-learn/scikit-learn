/*
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 * Copyright (c) 2009 Helge Bahmann
 * Copyright (c) 2013 Tim Blechmann
 * Copyright (c) 2014 Andrey Semashev
 */
/*!
 * \file   atomic/detail/ops_gcc_arm_common.hpp
 *
 * This header contains basic utilities for gcc ARM backend.
 */

#ifndef BOOST_ATOMIC_DETAIL_OPS_GCC_ARM_COMMON_HPP_INCLUDED_
#define BOOST_ATOMIC_DETAIL_OPS_GCC_ARM_COMMON_HPP_INCLUDED_

#include <boost/cstdint.hpp>
#include <boost/memory_order.hpp>
#include <boost/atomic/detail/config.hpp>

#ifdef BOOST_HAS_PRAGMA_ONCE
#pragma once
#endif

namespace boost {
namespace atomics {
namespace detail {

// A memory barrier is effected using a "co-processor 15" instruction,
// though a separate assembler mnemonic is available for it in v7.
//
// "Thumb 1" is a subset of the ARM instruction set that uses a 16-bit encoding.  It
// doesn't include all instructions and in particular it doesn't include the co-processor
// instruction used for the memory barrier or the load-locked/store-conditional
// instructions.  So, if we're compiling in "Thumb 1" mode, we need to wrap all of our
// asm blocks with code to temporarily change to ARM mode.
//
// You can only change between ARM and Thumb modes when branching using the bx instruction.
// bx takes an address specified in a register.  The least significant bit of the address
// indicates the mode, so 1 is added to indicate that the destination code is Thumb.
// A temporary register is needed for the address and is passed as an argument to these
// macros.  It must be one of the "low" registers accessible to Thumb code, specified
// using the "l" attribute in the asm statement.
//
// Architecture v7 introduces "Thumb 2", which does include (almost?) all of the ARM
// instruction set.  (Actually, there was an extension of v6 called v6T2 which supported
// "Thumb 2" mode, but its architecture manual is no longer available, referring to v7.)
// So in v7 we don't need to change to ARM mode; we can write "universal
// assembler" which will assemble to Thumb 2 or ARM code as appropriate.  The only thing
// we need to do to make this "universal" assembler mode work is to insert "IT" instructions
// to annotate the conditional instructions.  These are ignored in other modes (e.g. v6),
// so they can always be present.

// A note about memory_order_consume. Technically, this architecture allows to avoid
// unnecessary memory barrier after consume load since it supports data dependency ordering.
// However, some compiler optimizations may break a seemingly valid code relying on data
// dependency tracking by injecting bogus branches to aid out of order execution.
// This may happen not only in Boost.Atomic code but also in user's code, which we have no
// control of. See this thread: http://lists.boost.org/Archives/boost/2014/06/213890.php.
// For this reason we promote memory_order_consume to memory_order_acquire.

#if defined(__thumb__) && !defined(__thumb2__)
#define BOOST_ATOMIC_DETAIL_ARM_ASM_START(TMPREG) "adr " #TMPREG ", 8f\n" "bx " #TMPREG "\n" ".arm\n" ".align 4\n" "8:\n"
#define BOOST_ATOMIC_DETAIL_ARM_ASM_END(TMPREG)   "adr " #TMPREG ", 9f + 1\n" "bx " #TMPREG "\n" ".thumb\n" ".align 2\n" "9:\n"
#define BOOST_ATOMIC_DETAIL_ARM_ASM_TMPREG_CONSTRAINT(var) "=&l" (var)
#else
// The tmpreg may be wasted in this case, which is non-optimal.
#define BOOST_ATOMIC_DETAIL_ARM_ASM_START(TMPREG)
#define BOOST_ATOMIC_DETAIL_ARM_ASM_END(TMPREG)
#define BOOST_ATOMIC_DETAIL_ARM_ASM_TMPREG_CONSTRAINT(var) "=&r" (var)
#endif

struct gcc_arm_operations_base
{
    static BOOST_CONSTEXPR_OR_CONST bool full_cas_based = false;
    static BOOST_CONSTEXPR_OR_CONST bool is_always_lock_free = true;

    static BOOST_FORCEINLINE void fence_before(memory_order order) BOOST_NOEXCEPT
    {
        if ((static_cast< unsigned int >(order) & static_cast< unsigned int >(memory_order_release)) != 0u)
            hardware_full_fence();
    }

    static BOOST_FORCEINLINE void fence_after(memory_order order) BOOST_NOEXCEPT
    {
        if ((static_cast< unsigned int >(order) & (static_cast< unsigned int >(memory_order_consume) | static_cast< unsigned int >(memory_order_acquire))) != 0u)
            hardware_full_fence();
    }

    static BOOST_FORCEINLINE void fence_after_store(memory_order order) BOOST_NOEXCEPT
    {
        if (order == memory_order_seq_cst)
            hardware_full_fence();
    }

    static BOOST_FORCEINLINE void hardware_full_fence() BOOST_NOEXCEPT
    {
#if defined(BOOST_ATOMIC_DETAIL_ARM_HAS_DMB)
        // Older binutils (supposedly, older than 2.21.1) didn't support symbolic or numeric arguments of the "dmb" instruction such as "ish" or "#11".
        // As a workaround we have to inject encoded bytes of the instruction. There are two encodings for the instruction: ARM and Thumb. See ARM Architecture Reference Manual, A8.8.43.
        // Since we cannot detect binutils version at compile time, we'll have to always use this hack.
        __asm__ __volatile__
        (
#if defined(__thumb2__)
            ".short 0xF3BF, 0x8F5B\n" // dmb ish
#else
            ".word 0xF57FF05B\n" // dmb ish
#endif
            :
            :
            : "memory"
        );
#else
        uint32_t tmp;
        __asm__ __volatile__
        (
            BOOST_ATOMIC_DETAIL_ARM_ASM_START(%0)
            "mcr\tp15, 0, r0, c7, c10, 5\n"
            BOOST_ATOMIC_DETAIL_ARM_ASM_END(%0)
            : "=&l" (tmp)
            :
            : "memory"
        );
#endif
    }
};

} // namespace detail
} // namespace atomics
} // namespace boost

#endif // BOOST_ATOMIC_DETAIL_OPS_GCC_ARM_COMMON_HPP_INCLUDED_
