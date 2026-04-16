/*
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 * Copyright (c) 2017 - 2018 Andrey Semashev
 */
/*!
 * \file   atomic/detail/extra_ops_gcc_arm.hpp
 *
 * This header contains implementation of the extra atomic operations for ARM.
 */

#ifndef BOOST_ATOMIC_DETAIL_EXTRA_OPS_GCC_ARM_HPP_INCLUDED_
#define BOOST_ATOMIC_DETAIL_EXTRA_OPS_GCC_ARM_HPP_INCLUDED_

#include <cstddef>
#include <boost/cstdint.hpp>
#include <boost/memory_order.hpp>
#include <boost/atomic/detail/config.hpp>
#include <boost/atomic/detail/platform.hpp>
#include <boost/atomic/detail/storage_type.hpp>
#include <boost/atomic/detail/extra_operations_fwd.hpp>
#include <boost/atomic/detail/extra_ops_generic.hpp>
#include <boost/atomic/detail/ops_gcc_arm_common.hpp>
#include <boost/atomic/capabilities.hpp>

#ifdef BOOST_HAS_PRAGMA_ONCE
#pragma once
#endif

namespace boost {
namespace atomics {
namespace detail {

template< typename Base >
struct gcc_arm_extra_operations_common :
    public Base
{
    typedef Base base_type;
    typedef typename base_type::storage_type storage_type;

    static BOOST_FORCEINLINE void opaque_negate(storage_type volatile& storage, memory_order order) BOOST_NOEXCEPT
    {
        base_type::fetch_negate(storage, order);
    }

    static BOOST_FORCEINLINE void opaque_complement(storage_type volatile& storage, memory_order order) BOOST_NOEXCEPT
    {
        base_type::fetch_complement(storage, order);
    }

    static BOOST_FORCEINLINE bool negate_and_test(storage_type volatile& storage, memory_order order) BOOST_NOEXCEPT
    {
        return !!base_type::negate(storage, order);
    }

    static BOOST_FORCEINLINE bool add_and_test(storage_type volatile& storage, storage_type v, memory_order order) BOOST_NOEXCEPT
    {
        return !!base_type::add(storage, v, order);
    }

    static BOOST_FORCEINLINE bool sub_and_test(storage_type volatile& storage, storage_type v, memory_order order) BOOST_NOEXCEPT
    {
        return !!base_type::sub(storage, v, order);
    }

    static BOOST_FORCEINLINE bool and_and_test(storage_type volatile& storage, storage_type v, memory_order order) BOOST_NOEXCEPT
    {
        return !!base_type::bitwise_and(storage, v, order);
    }

    static BOOST_FORCEINLINE bool or_and_test(storage_type volatile& storage, storage_type v, memory_order order) BOOST_NOEXCEPT
    {
        return !!base_type::bitwise_or(storage, v, order);
    }

    static BOOST_FORCEINLINE bool xor_and_test(storage_type volatile& storage, storage_type v, memory_order order) BOOST_NOEXCEPT
    {
        return !!base_type::bitwise_xor(storage, v, order);
    }

    static BOOST_FORCEINLINE bool complement_and_test(storage_type volatile& storage, memory_order order) BOOST_NOEXCEPT
    {
        return !!base_type::bitwise_complement(storage, order);
    }
};

template< typename Base, std::size_t Size, bool Signed >
struct gcc_arm_extra_operations;

#if defined(BOOST_ATOMIC_DETAIL_ARM_HAS_LDREXB_STREXB)

template< typename Base, bool Signed >
struct gcc_arm_extra_operations< Base, 1u, Signed > :
    public generic_extra_operations< Base, 1u, Signed >
{
    typedef generic_extra_operations< Base, 1u, Signed > base_type;
    typedef typename base_type::storage_type storage_type;
    typedef typename make_storage_type< 4u >::type extended_storage_type;

    static BOOST_FORCEINLINE storage_type fetch_negate(storage_type volatile& storage, memory_order order) BOOST_NOEXCEPT
    {
        gcc_arm_operations_base::fence_before(order);
        uint32_t tmp;
        extended_storage_type original, result;
        __asm__ __volatile__
        (
            BOOST_ATOMIC_DETAIL_ARM_ASM_START(%[tmp])
            "1:\n"
            "ldrexb   %[original], %[storage]\n"           // original = zero_extend(*(&storage))
            "rsb      %[result], %[original], #0\n"        // result = 0 - original
            "strexb   %[tmp], %[result], %[storage]\n"     // *(&storage) = result, tmp = store failed
            "teq      %[tmp], #0\n"                        // flags = tmp==0
            "bne      1b\n"                                // if (!flags.equal) goto retry
            BOOST_ATOMIC_DETAIL_ARM_ASM_END(%[tmp])
            : [original] "=&r" (original),  // %0
              [result] "=&r" (result),      // %1
              [tmp] "=&l" (tmp),            // %2
              [storage] "+Q" (storage)      // %3
            :
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC
        );
        gcc_arm_operations_base::fence_after(order);
        return static_cast< storage_type >(original);
    }

    static BOOST_FORCEINLINE storage_type negate(storage_type volatile& storage, memory_order order) BOOST_NOEXCEPT
    {
        gcc_arm_operations_base::fence_before(order);
        uint32_t tmp;
        extended_storage_type original, result;
        __asm__ __volatile__
        (
            BOOST_ATOMIC_DETAIL_ARM_ASM_START(%[tmp])
            "1:\n"
            "ldrexb   %[original], %[storage]\n"           // original = zero_extend(*(&storage))
            "rsb      %[result], %[original], #0\n"        // result = 0 - original
            "strexb   %[tmp], %[result], %[storage]\n"     // *(&storage) = result, tmp = store failed
            "teq      %[tmp], #0\n"                        // flags = tmp==0
            "bne      1b\n"                                // if (!flags.equal) goto retry
            BOOST_ATOMIC_DETAIL_ARM_ASM_END(%[tmp])
            : [original] "=&r" (original),  // %0
              [result] "=&r" (result),      // %1
              [tmp] "=&l" (tmp),            // %2
              [storage] "+Q" (storage)      // %3
            :
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC
        );
        gcc_arm_operations_base::fence_after(order);
        return static_cast< storage_type >(result);
    }

    static BOOST_FORCEINLINE storage_type add(storage_type volatile& storage, storage_type v, memory_order order) BOOST_NOEXCEPT
    {
        gcc_arm_operations_base::fence_before(order);
        uint32_t tmp;
        extended_storage_type original, result;
        __asm__ __volatile__
        (
            BOOST_ATOMIC_DETAIL_ARM_ASM_START(%[tmp])
            "1:\n"
            "ldrexb   %[original], %[storage]\n"           // original = zero_extend(*(&storage))
            "add      %[result], %[original], %[value]\n"  // result = original + value
            "strexb   %[tmp], %[result], %[storage]\n"     // *(&storage) = result, tmp = store failed
            "teq      %[tmp], #0\n"                        // flags = tmp==0
            "bne      1b\n"                                // if (!flags.equal) goto retry
            BOOST_ATOMIC_DETAIL_ARM_ASM_END(%[tmp])
            : [original] "=&r" (original),  // %0
              [result] "=&r" (result),      // %1
              [tmp] "=&l" (tmp),            // %2
              [storage] "+Q" (storage)      // %3
            : [value] "Ir" (v)              // %4
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC
        );
        gcc_arm_operations_base::fence_after(order);
        return static_cast< storage_type >(result);
    }

    static BOOST_FORCEINLINE storage_type sub(storage_type volatile& storage, storage_type v, memory_order order) BOOST_NOEXCEPT
    {
        gcc_arm_operations_base::fence_before(order);
        uint32_t tmp;
        extended_storage_type original, result;
        __asm__ __volatile__
        (
            BOOST_ATOMIC_DETAIL_ARM_ASM_START(%[tmp])
            "1:\n"
            "ldrexb   %[original], %[storage]\n"           // original = zero_extend(*(&storage))
            "sub      %[result], %[original], %[value]\n"  // result = original - value
            "strexb   %[tmp], %[result], %[storage]\n"     // *(&storage) = result, tmp = store failed
            "teq      %[tmp], #0\n"                        // flags = tmp==0
            "bne      1b\n"                                // if (!flags.equal) goto retry
            BOOST_ATOMIC_DETAIL_ARM_ASM_END(%[tmp])
            : [original] "=&r" (original),  // %0
              [result] "=&r" (result),      // %1
              [tmp] "=&l" (tmp),            // %2
              [storage] "+Q" (storage)      // %3
            : [value] "Ir" (v)              // %4
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC
        );
        gcc_arm_operations_base::fence_after(order);
        return static_cast< storage_type >(result);
    }

    static BOOST_FORCEINLINE storage_type bitwise_and(storage_type volatile& storage, storage_type v, memory_order order) BOOST_NOEXCEPT
    {
        gcc_arm_operations_base::fence_before(order);
        uint32_t tmp;
        extended_storage_type original, result;
        __asm__ __volatile__
        (
            BOOST_ATOMIC_DETAIL_ARM_ASM_START(%[tmp])
            "1:\n"
            "ldrexb   %[original], %[storage]\n"           // original = zero_extend(*(&storage))
            "and      %[result], %[original], %[value]\n"  // result = original & value
            "strexb   %[tmp], %[result], %[storage]\n"     // *(&storage) = result, tmp = store failed
            "teq      %[tmp], #0\n"                        // flags = tmp==0
            "bne      1b\n"                                // if (!flags.equal) goto retry
            BOOST_ATOMIC_DETAIL_ARM_ASM_END(%[tmp])
            : [original] "=&r" (original),  // %0
              [result] "=&r" (result),      // %1
              [tmp] "=&l" (tmp),            // %2
              [storage] "+Q" (storage)      // %3
            : [value] "Ir" (v)              // %4
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC
        );
        gcc_arm_operations_base::fence_after(order);
        return static_cast< storage_type >(result);
    }

    static BOOST_FORCEINLINE storage_type bitwise_or(storage_type volatile& storage, storage_type v, memory_order order) BOOST_NOEXCEPT
    {
        gcc_arm_operations_base::fence_before(order);
        uint32_t tmp;
        extended_storage_type original, result;
        __asm__ __volatile__
        (
            BOOST_ATOMIC_DETAIL_ARM_ASM_START(%[tmp])
            "1:\n"
            "ldrexb   %[original], %[storage]\n"           // original = zero_extend(*(&storage))
            "orr      %[result], %[original], %[value]\n"  // result = original | value
            "strexb   %[tmp], %[result], %[storage]\n"     // *(&storage) = result, tmp = store failed
            "teq      %[tmp], #0\n"                        // flags = tmp==0
            "bne      1b\n"                                // if (!flags.equal) goto retry
            BOOST_ATOMIC_DETAIL_ARM_ASM_END(%[tmp])
            : [original] "=&r" (original),  // %0
              [result] "=&r" (result),      // %1
              [tmp] "=&l" (tmp),            // %2
              [storage] "+Q" (storage)      // %3
            : [value] "Ir" (v)              // %4
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC
        );
        gcc_arm_operations_base::fence_after(order);
        return static_cast< storage_type >(result);
    }

    static BOOST_FORCEINLINE storage_type bitwise_xor(storage_type volatile& storage, storage_type v, memory_order order) BOOST_NOEXCEPT
    {
        gcc_arm_operations_base::fence_before(order);
        uint32_t tmp;
        extended_storage_type original, result;
        __asm__ __volatile__
        (
            BOOST_ATOMIC_DETAIL_ARM_ASM_START(%[tmp])
            "1:\n"
            "ldrexb   %[original], %[storage]\n"           // original = zero_extend(*(&storage))
            "eor      %[result], %[original], %[value]\n"  // result = original ^ value
            "strexb   %[tmp], %[result], %[storage]\n"     // *(&storage) = result, tmp = store failed
            "teq      %[tmp], #0\n"                        // flags = tmp==0
            "bne      1b\n"                                // if (!flags.equal) goto retry
            BOOST_ATOMIC_DETAIL_ARM_ASM_END(%[tmp])
            : [original] "=&r" (original),  // %0
              [result] "=&r" (result),      // %1
              [tmp] "=&l" (tmp),            // %2
              [storage] "+Q" (storage)      // %3
            : [value] "Ir" (v)              // %4
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC
        );
        gcc_arm_operations_base::fence_after(order);
        return static_cast< storage_type >(result);
    }

    static BOOST_FORCEINLINE storage_type fetch_complement(storage_type volatile& storage, memory_order order) BOOST_NOEXCEPT
    {
        gcc_arm_operations_base::fence_before(order);
        uint32_t tmp;
        extended_storage_type original, result;
        __asm__ __volatile__
        (
            BOOST_ATOMIC_DETAIL_ARM_ASM_START(%[tmp])
            "1:\n"
            "ldrexb   %[original], %[storage]\n"           // original = zero_extend(*(&storage))
            "mvn      %[result], %[original]\n"            // result = NOT original
            "strexb   %[tmp], %[result], %[storage]\n"     // *(&storage) = result, tmp = store failed
            "teq      %[tmp], #0\n"                        // flags = tmp==0
            "bne      1b\n"                                // if (!flags.equal) goto retry
            BOOST_ATOMIC_DETAIL_ARM_ASM_END(%[tmp])
            : [original] "=&r" (original),  // %0
              [result] "=&r" (result),      // %1
              [tmp] "=&l" (tmp),            // %2
              [storage] "+Q" (storage)      // %3
            :
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC
        );
        gcc_arm_operations_base::fence_after(order);
        return static_cast< storage_type >(original);
    }

    static BOOST_FORCEINLINE storage_type bitwise_complement(storage_type volatile& storage, memory_order order) BOOST_NOEXCEPT
    {
        gcc_arm_operations_base::fence_before(order);
        uint32_t tmp;
        extended_storage_type original, result;
        __asm__ __volatile__
        (
            BOOST_ATOMIC_DETAIL_ARM_ASM_START(%[tmp])
            "1:\n"
            "ldrexb   %[original], %[storage]\n"           // original = zero_extend(*(&storage))
            "mvn      %[result], %[original]\n"            // result = NOT original
            "strexb   %[tmp], %[result], %[storage]\n"     // *(&storage) = result, tmp = store failed
            "teq      %[tmp], #0\n"                        // flags = tmp==0
            "bne      1b\n"                                // if (!flags.equal) goto retry
            BOOST_ATOMIC_DETAIL_ARM_ASM_END(%[tmp])
            : [original] "=&r" (original),  // %0
              [result] "=&r" (result),      // %1
              [tmp] "=&l" (tmp),            // %2
              [storage] "+Q" (storage)      // %3
            :
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC
        );
        gcc_arm_operations_base::fence_after(order);
        return static_cast< storage_type >(result);
    }
};

template< typename Base, bool Signed >
struct extra_operations< Base, 1u, Signed, true > :
    public gcc_arm_extra_operations_common< gcc_arm_extra_operations< Base, 1u, Signed > >
{
};

#endif // defined(BOOST_ATOMIC_DETAIL_ARM_HAS_LDREXB_STREXB)

#if defined(BOOST_ATOMIC_DETAIL_ARM_HAS_LDREXH_STREXH)

template< typename Base, bool Signed >
struct gcc_arm_extra_operations< Base, 2u, Signed > :
    public generic_extra_operations< Base, 2u, Signed >
{
    typedef generic_extra_operations< Base, 2u, Signed > base_type;
    typedef typename base_type::storage_type storage_type;
    typedef typename make_storage_type< 4u >::type extended_storage_type;

    static BOOST_FORCEINLINE storage_type fetch_negate(storage_type volatile& storage, memory_order order) BOOST_NOEXCEPT
    {
        gcc_arm_operations_base::fence_before(order);
        uint32_t tmp;
        extended_storage_type original, result;
        __asm__ __volatile__
        (
            BOOST_ATOMIC_DETAIL_ARM_ASM_START(%[tmp])
            "1:\n"
            "ldrexh   %[original], %[storage]\n"           // original = zero_extend(*(&storage))
            "rsb      %[result], %[original], #0\n"        // result = 0 - original
            "strexh   %[tmp], %[result], %[storage]\n"     // *(&storage) = result, tmp = store failed
            "teq      %[tmp], #0\n"                        // flags = tmp==0
            "bne      1b\n"                                // if (!flags.equal) goto retry
            BOOST_ATOMIC_DETAIL_ARM_ASM_END(%[tmp])
            : [original] "=&r" (original),  // %0
              [result] "=&r" (result),      // %1
              [tmp] "=&l" (tmp),            // %2
              [storage] "+Q" (storage)      // %3
            :
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC
        );
        gcc_arm_operations_base::fence_after(order);
        return static_cast< storage_type >(original);
    }

    static BOOST_FORCEINLINE storage_type negate(storage_type volatile& storage, memory_order order) BOOST_NOEXCEPT
    {
        gcc_arm_operations_base::fence_before(order);
        uint32_t tmp;
        extended_storage_type original, result;
        __asm__ __volatile__
        (
            BOOST_ATOMIC_DETAIL_ARM_ASM_START(%[tmp])
            "1:\n"
            "ldrexh   %[original], %[storage]\n"           // original = zero_extend(*(&storage))
            "rsb      %[result], %[original], #0\n"        // result = 0 - original
            "strexh   %[tmp], %[result], %[storage]\n"     // *(&storage) = result, tmp = store failed
            "teq      %[tmp], #0\n"                        // flags = tmp==0
            "bne      1b\n"                                // if (!flags.equal) goto retry
            BOOST_ATOMIC_DETAIL_ARM_ASM_END(%[tmp])
            : [original] "=&r" (original),  // %0
              [result] "=&r" (result),      // %1
              [tmp] "=&l" (tmp),            // %2
              [storage] "+Q" (storage)      // %3
            :
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC
        );
        gcc_arm_operations_base::fence_after(order);
        return static_cast< storage_type >(result);
    }

    static BOOST_FORCEINLINE storage_type add(storage_type volatile& storage, storage_type v, memory_order order) BOOST_NOEXCEPT
    {
        gcc_arm_operations_base::fence_before(order);
        uint32_t tmp;
        extended_storage_type original, result;
        __asm__ __volatile__
        (
            BOOST_ATOMIC_DETAIL_ARM_ASM_START(%[tmp])
            "1:\n"
            "ldrexh   %[original], %[storage]\n"           // original = zero_extend(*(&storage))
            "add      %[result], %[original], %[value]\n"  // result = original + value
            "strexh   %[tmp], %[result], %[storage]\n"     // *(&storage) = result, tmp = store failed
            "teq      %[tmp], #0\n"                        // flags = tmp==0
            "bne      1b\n"                                // if (!flags.equal) goto retry
            BOOST_ATOMIC_DETAIL_ARM_ASM_END(%[tmp])
            : [original] "=&r" (original),  // %0
              [result] "=&r" (result),      // %1
              [tmp] "=&l" (tmp),            // %2
              [storage] "+Q" (storage)      // %3
            : [value] "Ir" (v)              // %4
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC
        );
        gcc_arm_operations_base::fence_after(order);
        return static_cast< storage_type >(result);
    }

    static BOOST_FORCEINLINE storage_type sub(storage_type volatile& storage, storage_type v, memory_order order) BOOST_NOEXCEPT
    {
        gcc_arm_operations_base::fence_before(order);
        uint32_t tmp;
        extended_storage_type original, result;
        __asm__ __volatile__
        (
            BOOST_ATOMIC_DETAIL_ARM_ASM_START(%[tmp])
            "1:\n"
            "ldrexh   %[original], %[storage]\n"           // original = zero_extend(*(&storage))
            "sub      %[result], %[original], %[value]\n"  // result = original - value
            "strexh   %[tmp], %[result], %[storage]\n"     // *(&storage) = result, tmp = store failed
            "teq      %[tmp], #0\n"                        // flags = tmp==0
            "bne      1b\n"                                // if (!flags.equal) goto retry
            BOOST_ATOMIC_DETAIL_ARM_ASM_END(%[tmp])
            : [original] "=&r" (original),  // %0
              [result] "=&r" (result),      // %1
              [tmp] "=&l" (tmp),            // %2
              [storage] "+Q" (storage)      // %3
            : [value] "Ir" (v)              // %4
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC
        );
        gcc_arm_operations_base::fence_after(order);
        return static_cast< storage_type >(result);
    }

    static BOOST_FORCEINLINE storage_type bitwise_and(storage_type volatile& storage, storage_type v, memory_order order) BOOST_NOEXCEPT
    {
        gcc_arm_operations_base::fence_before(order);
        uint32_t tmp;
        extended_storage_type original, result;
        __asm__ __volatile__
        (
            BOOST_ATOMIC_DETAIL_ARM_ASM_START(%[tmp])
            "1:\n"
            "ldrexh   %[original], %[storage]\n"           // original = zero_extend(*(&storage))
            "and      %[result], %[original], %[value]\n"  // result = original & value
            "strexh   %[tmp], %[result], %[storage]\n"     // *(&storage) = result, tmp = store failed
            "teq      %[tmp], #0\n"                        // flags = tmp==0
            "bne      1b\n"                                // if (!flags.equal) goto retry
            BOOST_ATOMIC_DETAIL_ARM_ASM_END(%[tmp])
            : [original] "=&r" (original),  // %0
              [result] "=&r" (result),      // %1
              [tmp] "=&l" (tmp),            // %2
              [storage] "+Q" (storage)      // %3
            : [value] "Ir" (v)              // %4
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC
        );
        gcc_arm_operations_base::fence_after(order);
        return static_cast< storage_type >(result);
    }

    static BOOST_FORCEINLINE storage_type bitwise_or(storage_type volatile& storage, storage_type v, memory_order order) BOOST_NOEXCEPT
    {
        gcc_arm_operations_base::fence_before(order);
        uint32_t tmp;
        extended_storage_type original, result;
        __asm__ __volatile__
        (
            BOOST_ATOMIC_DETAIL_ARM_ASM_START(%[tmp])
            "1:\n"
            "ldrexh   %[original], %[storage]\n"           // original = zero_extend(*(&storage))
            "orr      %[result], %[original], %[value]\n"  // result = original | value
            "strexh   %[tmp], %[result], %[storage]\n"     // *(&storage) = result, tmp = store failed
            "teq      %[tmp], #0\n"                        // flags = tmp==0
            "bne      1b\n"                                // if (!flags.equal) goto retry
            BOOST_ATOMIC_DETAIL_ARM_ASM_END(%[tmp])
            : [original] "=&r" (original),  // %0
              [result] "=&r" (result),      // %1
              [tmp] "=&l" (tmp),            // %2
              [storage] "+Q" (storage)      // %3
            : [value] "Ir" (v)              // %4
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC
        );
        gcc_arm_operations_base::fence_after(order);
        return static_cast< storage_type >(result);
    }

    static BOOST_FORCEINLINE storage_type bitwise_xor(storage_type volatile& storage, storage_type v, memory_order order) BOOST_NOEXCEPT
    {
        gcc_arm_operations_base::fence_before(order);
        uint32_t tmp;
        extended_storage_type original, result;
        __asm__ __volatile__
        (
            BOOST_ATOMIC_DETAIL_ARM_ASM_START(%[tmp])
            "1:\n"
            "ldrexh   %[original], %[storage]\n"           // original = zero_extend(*(&storage))
            "eor      %[result], %[original], %[value]\n"  // result = original ^ value
            "strexh   %[tmp], %[result], %[storage]\n"     // *(&storage) = result, tmp = store failed
            "teq      %[tmp], #0\n"                        // flags = tmp==0
            "bne      1b\n"                                // if (!flags.equal) goto retry
            BOOST_ATOMIC_DETAIL_ARM_ASM_END(%[tmp])
            : [original] "=&r" (original),  // %0
              [result] "=&r" (result),      // %1
              [tmp] "=&l" (tmp),            // %2
              [storage] "+Q" (storage)      // %3
            : [value] "Ir" (v)              // %4
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC
        );
        gcc_arm_operations_base::fence_after(order);
        return static_cast< storage_type >(result);
    }

    static BOOST_FORCEINLINE storage_type fetch_complement(storage_type volatile& storage, memory_order order) BOOST_NOEXCEPT
    {
        gcc_arm_operations_base::fence_before(order);
        uint32_t tmp;
        extended_storage_type original, result;
        __asm__ __volatile__
        (
            BOOST_ATOMIC_DETAIL_ARM_ASM_START(%[tmp])
            "1:\n"
            "ldrexh   %[original], %[storage]\n"           // original = zero_extend(*(&storage))
            "mvn      %[result], %[original]\n"            // result = NOT original
            "strexh   %[tmp], %[result], %[storage]\n"     // *(&storage) = result, tmp = store failed
            "teq      %[tmp], #0\n"                        // flags = tmp==0
            "bne      1b\n"                                // if (!flags.equal) goto retry
            BOOST_ATOMIC_DETAIL_ARM_ASM_END(%[tmp])
            : [original] "=&r" (original),  // %0
              [result] "=&r" (result),      // %1
              [tmp] "=&l" (tmp),            // %2
              [storage] "+Q" (storage)      // %3
            :
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC
        );
        gcc_arm_operations_base::fence_after(order);
        return static_cast< storage_type >(original);
    }

    static BOOST_FORCEINLINE storage_type bitwise_complement(storage_type volatile& storage, memory_order order) BOOST_NOEXCEPT
    {
        gcc_arm_operations_base::fence_before(order);
        uint32_t tmp;
        extended_storage_type original, result;
        __asm__ __volatile__
        (
            BOOST_ATOMIC_DETAIL_ARM_ASM_START(%[tmp])
            "1:\n"
            "ldrexh   %[original], %[storage]\n"           // original = zero_extend(*(&storage))
            "mvn      %[result], %[original]\n"            // result = NOT original
            "strexh   %[tmp], %[result], %[storage]\n"     // *(&storage) = result, tmp = store failed
            "teq      %[tmp], #0\n"                        // flags = tmp==0
            "bne      1b\n"                                // if (!flags.equal) goto retry
            BOOST_ATOMIC_DETAIL_ARM_ASM_END(%[tmp])
            : [original] "=&r" (original),  // %0
              [result] "=&r" (result),      // %1
              [tmp] "=&l" (tmp),            // %2
              [storage] "+Q" (storage)      // %3
            :
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC
        );
        gcc_arm_operations_base::fence_after(order);
        return static_cast< storage_type >(result);
    }
};

template< typename Base, bool Signed >
struct extra_operations< Base, 2u, Signed, true > :
    public gcc_arm_extra_operations_common< gcc_arm_extra_operations< Base, 2u, Signed > >
{
};

#endif // defined(BOOST_ATOMIC_DETAIL_ARM_HAS_LDREXH_STREXH)

template< typename Base, bool Signed >
struct gcc_arm_extra_operations< Base, 4u, Signed > :
    public generic_extra_operations< Base, 4u, Signed >
{
    typedef generic_extra_operations< Base, 4u, Signed > base_type;
    typedef typename base_type::storage_type storage_type;

    static BOOST_FORCEINLINE storage_type fetch_negate(storage_type volatile& storage, memory_order order) BOOST_NOEXCEPT
    {
        gcc_arm_operations_base::fence_before(order);
        uint32_t tmp;
        storage_type original, result;
        __asm__ __volatile__
        (
            BOOST_ATOMIC_DETAIL_ARM_ASM_START(%[tmp])
            "1:\n"
            "ldrex    %[original], %[storage]\n"           // original = *(&storage)
            "rsb      %[result], %[original], #0\n"        // result = 0 - original
            "strex    %[tmp], %[result], %[storage]\n"     // *(&storage) = result, tmp = store failed
            "teq      %[tmp], #0\n"                        // flags = tmp==0
            "bne      1b\n"                                // if (!flags.equal) goto retry
            BOOST_ATOMIC_DETAIL_ARM_ASM_END(%[tmp])
            : [original] "=&r" (original),  // %0
              [result] "=&r" (result),      // %1
              [tmp] "=&l" (tmp),            // %2
              [storage] "+Q" (storage)      // %3
            :
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC
        );
        gcc_arm_operations_base::fence_after(order);
        return original;
    }

    static BOOST_FORCEINLINE storage_type negate(storage_type volatile& storage, memory_order order) BOOST_NOEXCEPT
    {
        gcc_arm_operations_base::fence_before(order);
        uint32_t tmp;
        storage_type original, result;
        __asm__ __volatile__
        (
            BOOST_ATOMIC_DETAIL_ARM_ASM_START(%[tmp])
            "1:\n"
            "ldrex    %[original], %[storage]\n"           // original = *(&storage)
            "rsb      %[result], %[original], #0\n"        // result = 0 - original
            "strex    %[tmp], %[result], %[storage]\n"     // *(&storage) = result, tmp = store failed
            "teq      %[tmp], #0\n"                        // flags = tmp==0
            "bne      1b\n"                                // if (!flags.equal) goto retry
            BOOST_ATOMIC_DETAIL_ARM_ASM_END(%[tmp])
            : [original] "=&r" (original),  // %0
              [result] "=&r" (result),      // %1
              [tmp] "=&l" (tmp),            // %2
              [storage] "+Q" (storage)      // %3
            :
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC
        );
        gcc_arm_operations_base::fence_after(order);
        return result;
    }

    static BOOST_FORCEINLINE storage_type add(storage_type volatile& storage, storage_type v, memory_order order) BOOST_NOEXCEPT
    {
        gcc_arm_operations_base::fence_before(order);
        uint32_t tmp;
        storage_type original, result;
        __asm__ __volatile__
        (
            BOOST_ATOMIC_DETAIL_ARM_ASM_START(%[tmp])
            "1:\n"
            "ldrex   %[original], %[storage]\n"           // original = *(&storage)
            "add     %[result], %[original], %[value]\n"  // result = original + value
            "strex   %[tmp], %[result], %[storage]\n"     // *(&storage) = result, tmp = store failed
            "teq     %[tmp], #0\n"                        // flags = tmp==0
            "bne     1b\n"                                // if (!flags.equal) goto retry
            BOOST_ATOMIC_DETAIL_ARM_ASM_END(%[tmp])
            : [original] "=&r" (original),  // %0
              [result] "=&r" (result),      // %1
              [tmp] "=&l" (tmp),            // %2
              [storage] "+Q" (storage)      // %3
            : [value] "Ir" (v)              // %4
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC
        );
        gcc_arm_operations_base::fence_after(order);
        return result;
    }

    static BOOST_FORCEINLINE storage_type sub(storage_type volatile& storage, storage_type v, memory_order order) BOOST_NOEXCEPT
    {
        gcc_arm_operations_base::fence_before(order);
        uint32_t tmp;
        storage_type original, result;
        __asm__ __volatile__
        (
            BOOST_ATOMIC_DETAIL_ARM_ASM_START(%[tmp])
            "1:\n"
            "ldrex   %[original], %[storage]\n"           // original = *(&storage)
            "sub     %[result], %[original], %[value]\n"  // result = original - value
            "strex   %[tmp], %[result], %[storage]\n"     // *(&storage) = result, tmp = store failed
            "teq     %[tmp], #0\n"                        // flags = tmp==0
            "bne     1b\n"                                // if (!flags.equal) goto retry
            BOOST_ATOMIC_DETAIL_ARM_ASM_END(%[tmp])
            : [original] "=&r" (original),  // %0
              [result] "=&r" (result),      // %1
              [tmp] "=&l" (tmp),            // %2
              [storage] "+Q" (storage)      // %3
            : [value] "Ir" (v)              // %4
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC
        );
        gcc_arm_operations_base::fence_after(order);
        return result;
    }

    static BOOST_FORCEINLINE storage_type bitwise_and(storage_type volatile& storage, storage_type v, memory_order order) BOOST_NOEXCEPT
    {
        gcc_arm_operations_base::fence_before(order);
        uint32_t tmp;
        storage_type original, result;
        __asm__ __volatile__
        (
            BOOST_ATOMIC_DETAIL_ARM_ASM_START(%[tmp])
            "1:\n"
            "ldrex   %[original], %[storage]\n"           // original = *(&storage)
            "and     %[result], %[original], %[value]\n"  // result = original & value
            "strex   %[tmp], %[result], %[storage]\n"     // *(&storage) = result, tmp = store failed
            "teq     %[tmp], #0\n"                        // flags = tmp==0
            "bne     1b\n"                                // if (!flags.equal) goto retry
            BOOST_ATOMIC_DETAIL_ARM_ASM_END(%[tmp])
            : [original] "=&r" (original),  // %0
              [result] "=&r" (result),      // %1
              [tmp] "=&l" (tmp),            // %2
              [storage] "+Q" (storage)      // %3
            : [value] "Ir" (v)              // %4
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC
        );
        gcc_arm_operations_base::fence_after(order);
        return result;
    }

    static BOOST_FORCEINLINE storage_type bitwise_or(storage_type volatile& storage, storage_type v, memory_order order) BOOST_NOEXCEPT
    {
        gcc_arm_operations_base::fence_before(order);
        uint32_t tmp;
        storage_type original, result;
        __asm__ __volatile__
        (
            BOOST_ATOMIC_DETAIL_ARM_ASM_START(%[tmp])
            "1:\n"
            "ldrex   %[original], %[storage]\n"           // original = *(&storage)
            "orr     %[result], %[original], %[value]\n"  // result = original | value
            "strex   %[tmp], %[result], %[storage]\n"     // *(&storage) = result, tmp = store failed
            "teq     %[tmp], #0\n"                        // flags = tmp==0
            "bne     1b\n"                                // if (!flags.equal) goto retry
            BOOST_ATOMIC_DETAIL_ARM_ASM_END(%[tmp])
            : [original] "=&r" (original),  // %0
              [result] "=&r" (result),      // %1
              [tmp] "=&l" (tmp),            // %2
              [storage] "+Q" (storage)      // %3
            : [value] "Ir" (v)              // %4
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC
        );
        gcc_arm_operations_base::fence_after(order);
        return result;
    }

    static BOOST_FORCEINLINE storage_type bitwise_xor(storage_type volatile& storage, storage_type v, memory_order order) BOOST_NOEXCEPT
    {
        gcc_arm_operations_base::fence_before(order);
        uint32_t tmp;
        storage_type original, result;
        __asm__ __volatile__
        (
            BOOST_ATOMIC_DETAIL_ARM_ASM_START(%[tmp])
            "1:\n"
            "ldrex   %[original], %[storage]\n"           // original = *(&storage)
            "eor     %[result], %[original], %[value]\n"  // result = original ^ value
            "strex   %[tmp], %[result], %[storage]\n"     // *(&storage) = result, tmp = store failed
            "teq     %[tmp], #0\n"                        // flags = tmp==0
            "bne     1b\n"                                // if (!flags.equal) goto retry
            BOOST_ATOMIC_DETAIL_ARM_ASM_END(%[tmp])
            : [original] "=&r" (original),  // %0
              [result] "=&r" (result),      // %1
              [tmp] "=&l" (tmp),            // %2
              [storage] "+Q" (storage)      // %3
            : [value] "Ir" (v)              // %4
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC
        );
        gcc_arm_operations_base::fence_after(order);
        return result;
    }

    static BOOST_FORCEINLINE storage_type fetch_complement(storage_type volatile& storage, memory_order order) BOOST_NOEXCEPT
    {
        gcc_arm_operations_base::fence_before(order);
        uint32_t tmp;
        storage_type original, result;
        __asm__ __volatile__
        (
            BOOST_ATOMIC_DETAIL_ARM_ASM_START(%[tmp])
            "1:\n"
            "ldrex    %[original], %[storage]\n"           // original = *(&storage)
            "mvn      %[result], %[original]\n"            // result = NOT original
            "strex    %[tmp], %[result], %[storage]\n"     // *(&storage) = result, tmp = store failed
            "teq      %[tmp], #0\n"                        // flags = tmp==0
            "bne      1b\n"                                // if (!flags.equal) goto retry
            BOOST_ATOMIC_DETAIL_ARM_ASM_END(%[tmp])
            : [original] "=&r" (original),  // %0
              [result] "=&r" (result),      // %1
              [tmp] "=&l" (tmp),            // %2
              [storage] "+Q" (storage)      // %3
            :
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC
        );
        gcc_arm_operations_base::fence_after(order);
        return original;
    }

    static BOOST_FORCEINLINE storage_type bitwise_complement(storage_type volatile& storage, memory_order order) BOOST_NOEXCEPT
    {
        gcc_arm_operations_base::fence_before(order);
        uint32_t tmp;
        storage_type original, result;
        __asm__ __volatile__
        (
            BOOST_ATOMIC_DETAIL_ARM_ASM_START(%[tmp])
            "1:\n"
            "ldrex    %[original], %[storage]\n"           // original = *(&storage)
            "mvn      %[result], %[original]\n"            // result = NOT original
            "strex    %[tmp], %[result], %[storage]\n"     // *(&storage) = result, tmp = store failed
            "teq      %[tmp], #0\n"                        // flags = tmp==0
            "bne      1b\n"                                // if (!flags.equal) goto retry
            BOOST_ATOMIC_DETAIL_ARM_ASM_END(%[tmp])
            : [original] "=&r" (original),  // %0
              [result] "=&r" (result),      // %1
              [tmp] "=&l" (tmp),            // %2
              [storage] "+Q" (storage)      // %3
            :
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC
        );
        gcc_arm_operations_base::fence_after(order);
        return result;
    }
};

template< typename Base, bool Signed >
struct extra_operations< Base, 4u, Signed, true > :
    public gcc_arm_extra_operations_common< gcc_arm_extra_operations< Base, 4u, Signed > >
{
};

#if defined(BOOST_ATOMIC_DETAIL_ARM_HAS_LDREXD_STREXD)

template< typename Base, bool Signed >
struct gcc_arm_extra_operations< Base, 8u, Signed > :
    public generic_extra_operations< Base, 8u, Signed >
{
    typedef generic_extra_operations< Base, 8u, Signed > base_type;
    typedef typename base_type::storage_type storage_type;

    static BOOST_FORCEINLINE storage_type fetch_negate(storage_type volatile& storage, memory_order order) BOOST_NOEXCEPT
    {
        gcc_arm_operations_base::fence_before(order);
        storage_type original, result;
        uint32_t tmp;
        __asm__ __volatile__
        (
            BOOST_ATOMIC_DETAIL_ARM_ASM_START(%0)
            "1:\n"
            "ldrexd  %1, %H1, [%3]\n"               // original = *(&storage)
            "mvn     %2, %1\n"                      // result = NOT original
            "mvn     %H2, %H1\n"
            "adds    %2, %2, #1\n"                  // result = result + 1
            "adc     %H2, %H2, #0\n"
            "strexd  %0, %2, %H2, [%3]\n"           // *(&storage) = result, tmp = store failed
            "teq     %0, #0\n"                      // flags = tmp==0
            "bne     1b\n"                          // if (!flags.equal) goto retry
            BOOST_ATOMIC_DETAIL_ARM_ASM_END(%0)
            : BOOST_ATOMIC_DETAIL_ARM_ASM_TMPREG_CONSTRAINT(tmp), // %0
              "=&r" (original),  // %1
              "=&r" (result)     // %2
            : "r" (&storage)     // %3
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC_COMMA "memory"
        );
        gcc_arm_operations_base::fence_after(order);
        return original;
    }

    static BOOST_FORCEINLINE storage_type negate(storage_type volatile& storage, memory_order order) BOOST_NOEXCEPT
    {
        gcc_arm_operations_base::fence_before(order);
        storage_type original, result;
        uint32_t tmp;
        __asm__ __volatile__
        (
            BOOST_ATOMIC_DETAIL_ARM_ASM_START(%0)
            "1:\n"
            "ldrexd  %1, %H1, [%3]\n"               // original = *(&storage)
            "mvn     %2, %1\n"                      // result = NOT original
            "mvn     %H2, %H1\n"
            "adds    %2, %2, #1\n"                  // result = result + 1
            "adc     %H2, %H2, #0\n"
            "strexd  %0, %2, %H2, [%3]\n"           // *(&storage) = result, tmp = store failed
            "teq     %0, #0\n"                      // flags = tmp==0
            "bne     1b\n"                          // if (!flags.equal) goto retry
            BOOST_ATOMIC_DETAIL_ARM_ASM_END(%0)
            : BOOST_ATOMIC_DETAIL_ARM_ASM_TMPREG_CONSTRAINT(tmp), // %0
              "=&r" (original),  // %1
              "=&r" (result)     // %2
            : "r" (&storage)     // %3
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC_COMMA "memory"
        );
        gcc_arm_operations_base::fence_after(order);
        return result;
    }

    static BOOST_FORCEINLINE storage_type add(storage_type volatile& storage, storage_type v, memory_order order) BOOST_NOEXCEPT
    {
        gcc_arm_operations_base::fence_before(order);
        storage_type original, result;
        uint32_t tmp;
        __asm__ __volatile__
        (
            BOOST_ATOMIC_DETAIL_ARM_ASM_START(%0)
            "1:\n"
            "ldrexd  %1, %H1, [%3]\n"               // original = *(&storage)
            "adds    %2, %1, %4\n"                  // result = original + value
            "adc     %H2, %H1, %H4\n"
            "strexd  %0, %2, %H2, [%3]\n"           // *(&storage) = result, tmp = store failed
            "teq     %0, #0\n"                      // flags = tmp==0
            "bne     1b\n"                          // if (!flags.equal) goto retry
            BOOST_ATOMIC_DETAIL_ARM_ASM_END(%0)
            : BOOST_ATOMIC_DETAIL_ARM_ASM_TMPREG_CONSTRAINT(tmp), // %0
              "=&r" (original),  // %1
              "=&r" (result)     // %2
            : "r" (&storage),    // %3
              "r" (v)            // %4
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC_COMMA "memory"
        );
        gcc_arm_operations_base::fence_after(order);
        return result;
    }

    static BOOST_FORCEINLINE storage_type sub(storage_type volatile& storage, storage_type v, memory_order order) BOOST_NOEXCEPT
    {
        gcc_arm_operations_base::fence_before(order);
        storage_type original, result;
        uint32_t tmp;
        __asm__ __volatile__
        (
            BOOST_ATOMIC_DETAIL_ARM_ASM_START(%0)
            "1:\n"
            "ldrexd  %1, %H1, [%3]\n"               // original = *(&storage)
            "subs    %2, %1, %4\n"                  // result = original - value
            "sbc     %H2, %H1, %H4\n"
            "strexd  %0, %2, %H2, [%3]\n"           // *(&storage) = result, tmp = store failed
            "teq     %0, #0\n"                      // flags = tmp==0
            "bne     1b\n"                          // if (!flags.equal) goto retry
            BOOST_ATOMIC_DETAIL_ARM_ASM_END(%0)
            : BOOST_ATOMIC_DETAIL_ARM_ASM_TMPREG_CONSTRAINT(tmp), // %0
              "=&r" (original),  // %1
              "=&r" (result)     // %2
            : "r" (&storage),    // %3
              "r" (v)            // %4
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC_COMMA "memory"
        );
        gcc_arm_operations_base::fence_after(order);
        return result;
    }

    static BOOST_FORCEINLINE storage_type bitwise_and(storage_type volatile& storage, storage_type v, memory_order order) BOOST_NOEXCEPT
    {
        gcc_arm_operations_base::fence_before(order);
        storage_type original, result;
        uint32_t tmp;
        __asm__ __volatile__
        (
            BOOST_ATOMIC_DETAIL_ARM_ASM_START(%0)
            "1:\n"
            "ldrexd  %1, %H1, [%3]\n"               // original = *(&storage)
            "and     %2, %1, %4\n"                  // result = original & value
            "and     %H2, %H1, %H4\n"
            "strexd  %0, %2, %H2, [%3]\n"           // *(&storage) = result, tmp = store failed
            "teq     %0, #0\n"                      // flags = tmp==0
            "bne     1b\n"                          // if (!flags.equal) goto retry
            BOOST_ATOMIC_DETAIL_ARM_ASM_END(%0)
            : BOOST_ATOMIC_DETAIL_ARM_ASM_TMPREG_CONSTRAINT(tmp), // %0
              "=&r" (original),  // %1
              "=&r" (result)     // %2
            : "r" (&storage),    // %3
              "r" (v)            // %4
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC_COMMA "memory"
        );
        gcc_arm_operations_base::fence_after(order);
        return result;
    }

    static BOOST_FORCEINLINE storage_type bitwise_or(storage_type volatile& storage, storage_type v, memory_order order) BOOST_NOEXCEPT
    {
        gcc_arm_operations_base::fence_before(order);
        storage_type original, result;
        uint32_t tmp;
        __asm__ __volatile__
        (
            BOOST_ATOMIC_DETAIL_ARM_ASM_START(%0)
            "1:\n"
            "ldrexd  %1, %H1, [%3]\n"               // original = *(&storage)
            "orr     %2, %1, %4\n"                  // result = original | value
            "orr     %H2, %H1, %H4\n"
            "strexd  %0, %2, %H2, [%3]\n"           // *(&storage) = result, tmp = store failed
            "teq     %0, #0\n"                      // flags = tmp==0
            "bne     1b\n"                          // if (!flags.equal) goto retry
            BOOST_ATOMIC_DETAIL_ARM_ASM_END(%0)
            : BOOST_ATOMIC_DETAIL_ARM_ASM_TMPREG_CONSTRAINT(tmp), // %0
              "=&r" (original),  // %1
              "=&r" (result)     // %2
            : "r" (&storage),    // %3
              "r" (v)            // %4
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC_COMMA "memory"
        );
        gcc_arm_operations_base::fence_after(order);
        return result;
    }

    static BOOST_FORCEINLINE storage_type bitwise_xor(storage_type volatile& storage, storage_type v, memory_order order) BOOST_NOEXCEPT
    {
        gcc_arm_operations_base::fence_before(order);
        storage_type original, result;
        uint32_t tmp;
        __asm__ __volatile__
        (
            BOOST_ATOMIC_DETAIL_ARM_ASM_START(%0)
            "1:\n"
            "ldrexd  %1, %H1, [%3]\n"               // original = *(&storage)
            "eor     %2, %1, %4\n"                  // result = original ^ value
            "eor     %H2, %H1, %H4\n"
            "strexd  %0, %2, %H2, [%3]\n"           // *(&storage) = result, tmp = store failed
            "teq     %0, #0\n"                      // flags = tmp==0
            "bne     1b\n"                          // if (!flags.equal) goto retry
            BOOST_ATOMIC_DETAIL_ARM_ASM_END(%0)
            : BOOST_ATOMIC_DETAIL_ARM_ASM_TMPREG_CONSTRAINT(tmp), // %0
              "=&r" (original),  // %1
              "=&r" (result)     // %2
            : "r" (&storage),    // %3
              "r" (v)            // %4
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC_COMMA "memory"
        );
        gcc_arm_operations_base::fence_after(order);
        return result;
    }

    static BOOST_FORCEINLINE storage_type fetch_complement(storage_type volatile& storage, memory_order order) BOOST_NOEXCEPT
    {
        gcc_arm_operations_base::fence_before(order);
        storage_type original, result;
        uint32_t tmp;
        __asm__ __volatile__
        (
            BOOST_ATOMIC_DETAIL_ARM_ASM_START(%0)
            "1:\n"
            "ldrexd  %1, %H1, [%3]\n"               // original = *(&storage)
            "mvn     %2, %1\n"                      // result = NOT original
            "mvn     %H2, %H1\n"
            "strexd  %0, %2, %H2, [%3]\n"           // *(&storage) = result, tmp = store failed
            "teq     %0, #0\n"                      // flags = tmp==0
            "bne     1b\n"                          // if (!flags.equal) goto retry
            BOOST_ATOMIC_DETAIL_ARM_ASM_END(%0)
            : BOOST_ATOMIC_DETAIL_ARM_ASM_TMPREG_CONSTRAINT(tmp), // %0
              "=&r" (original),  // %1
              "=&r" (result)     // %2
            : "r" (&storage)     // %3
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC_COMMA "memory"
        );
        gcc_arm_operations_base::fence_after(order);
        return original;
    }

    static BOOST_FORCEINLINE storage_type bitwise_complement(storage_type volatile& storage, memory_order order) BOOST_NOEXCEPT
    {
        gcc_arm_operations_base::fence_before(order);
        storage_type original, result;
        uint32_t tmp;
        __asm__ __volatile__
        (
            BOOST_ATOMIC_DETAIL_ARM_ASM_START(%0)
            "1:\n"
            "ldrexd  %1, %H1, [%3]\n"               // original = *(&storage)
            "mvn     %2, %1\n"                      // result = NOT original
            "mvn     %H2, %H1\n"
            "strexd  %0, %2, %H2, [%3]\n"           // *(&storage) = result, tmp = store failed
            "teq     %0, #0\n"                      // flags = tmp==0
            "bne     1b\n"                          // if (!flags.equal) goto retry
            BOOST_ATOMIC_DETAIL_ARM_ASM_END(%0)
            : BOOST_ATOMIC_DETAIL_ARM_ASM_TMPREG_CONSTRAINT(tmp), // %0
              "=&r" (original),  // %1
              "=&r" (result)     // %2
            : "r" (&storage)     // %3
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC_COMMA "memory"
        );
        gcc_arm_operations_base::fence_after(order);
        return result;
    }
};

template< typename Base, bool Signed >
struct extra_operations< Base, 8u, Signed, true > :
    public gcc_arm_extra_operations_common< gcc_arm_extra_operations< Base, 8u, Signed > >
{
};

#endif // defined(BOOST_ATOMIC_DETAIL_ARM_HAS_LDREXD_STREXD)

} // namespace detail
} // namespace atomics
} // namespace boost

#endif // BOOST_ATOMIC_DETAIL_EXTRA_OPS_GCC_ARM_HPP_INCLUDED_
