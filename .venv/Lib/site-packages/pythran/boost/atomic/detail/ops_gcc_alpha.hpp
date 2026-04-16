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
 * \file   atomic/detail/ops_gcc_alpha.hpp
 *
 * This header contains implementation of the \c operations template.
 */

#ifndef BOOST_ATOMIC_DETAIL_OPS_GCC_ALPHA_HPP_INCLUDED_
#define BOOST_ATOMIC_DETAIL_OPS_GCC_ALPHA_HPP_INCLUDED_

#include <cstddef>
#include <boost/memory_order.hpp>
#include <boost/atomic/detail/config.hpp>
#include <boost/atomic/detail/storage_type.hpp>
#include <boost/atomic/detail/operations_fwd.hpp>
#include <boost/atomic/capabilities.hpp>

#ifdef BOOST_HAS_PRAGMA_ONCE
#pragma once
#endif

namespace boost {
namespace atomics {
namespace detail {

/*
  Refer to http://h71000.www7.hp.com/doc/82final/5601/5601pro_004.html
  (HP OpenVMS systems documentation) and the Alpha Architecture Reference Manual.
 */

/*
    NB: The most natural thing would be to write the increment/decrement
    operators along the following lines:

    __asm__ __volatile__
    (
        "1: ldl_l %0,%1 \n"
        "addl %0,1,%0 \n"
        "stl_c %0,%1 \n"
        "beq %0,1b\n"
        : "=&b" (tmp)
        : "m" (value)
        : "cc"
    );

    However according to the comments on the HP website and matching
    comments in the Linux kernel sources this defies branch prediction,
    as the cpu assumes that backward branches are always taken; so
    instead copy the trick from the Linux kernel, introduce a forward
    branch and back again.

    I have, however, had a hard time measuring the difference between
    the two versions in microbenchmarks -- I am leaving it in nevertheless
    as it apparently does not hurt either.
*/

struct gcc_alpha_operations_base
{
    static BOOST_CONSTEXPR_OR_CONST bool full_cas_based = false;
    static BOOST_CONSTEXPR_OR_CONST bool is_always_lock_free = true;

    static BOOST_FORCEINLINE void fence_before(memory_order order) BOOST_NOEXCEPT
    {
        if ((static_cast< unsigned int >(order) & static_cast< unsigned int >(memory_order_release)) != 0u)
            __asm__ __volatile__ ("mb" ::: "memory");
    }

    static BOOST_FORCEINLINE void fence_after(memory_order order) BOOST_NOEXCEPT
    {
        if ((static_cast< unsigned int >(order) & (static_cast< unsigned int >(memory_order_consume) | static_cast< unsigned int >(memory_order_acquire))) != 0u)
            __asm__ __volatile__ ("mb" ::: "memory");
    }

    static BOOST_FORCEINLINE void fence_after_store(memory_order order) BOOST_NOEXCEPT
    {
        if (order == memory_order_seq_cst)
            __asm__ __volatile__ ("mb" ::: "memory");
    }
};


template< bool Signed >
struct operations< 4u, Signed > :
    public gcc_alpha_operations_base
{
    typedef typename make_storage_type< 4u >::type storage_type;
    typedef typename make_storage_type< 4u >::aligned aligned_storage_type;

    static BOOST_CONSTEXPR_OR_CONST std::size_t storage_size = 4u;
    static BOOST_CONSTEXPR_OR_CONST bool is_signed = Signed;

    static BOOST_FORCEINLINE void store(storage_type volatile& storage, storage_type v, memory_order order) BOOST_NOEXCEPT
    {
        fence_before(order);
        storage = v;
        fence_after_store(order);
    }

    static BOOST_FORCEINLINE storage_type load(storage_type const volatile& storage, memory_order order) BOOST_NOEXCEPT
    {
        storage_type v = storage;
        fence_after(order);
        return v;
    }

    static BOOST_FORCEINLINE storage_type exchange(storage_type volatile& storage, storage_type v, memory_order order) BOOST_NOEXCEPT
    {
        storage_type original, tmp;
        fence_before(order);
        __asm__ __volatile__
        (
            "1:\n"
            "mov %3, %1\n"
            "ldl_l %0, %2\n"
            "stl_c %1, %2\n"
            "beq %1, 2f\n"

            ".subsection 2\n"
            "2: br 1b\n"
            ".previous\n"

            : "=&r" (original),  // %0
              "=&r" (tmp)        // %1
            : "m" (storage),     // %2
              "r" (v)            // %3
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC
        );
        fence_after(order);
        return original;
    }

    static BOOST_FORCEINLINE bool compare_exchange_weak(
        storage_type volatile& storage, storage_type& expected, storage_type desired, memory_order success_order, memory_order failure_order) BOOST_NOEXCEPT
    {
        fence_before(success_order);
        int success;
        storage_type current;
        __asm__ __volatile__
        (
            "1:\n"
            "ldl_l %2, %4\n"                // current = *(&storage)
            "cmpeq %2, %0, %3\n"            // success = current == expected
            "mov %2, %0\n"                  // expected = current
            "beq %3, 2f\n"                  // if (success == 0) goto end
            "stl_c %1, %4\n"                // storage = desired; desired = store succeeded
            "mov %1, %3\n"                  // success = desired
            "2:\n"
            : "+&r" (expected),  // %0
              "+&r" (desired),   // %1
              "=&r" (current),   // %2
              "=&r" (success)    // %3
            : "m" (storage)      // %4
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC
        );
        if (success)
            fence_after(success_order);
        else
            fence_after(failure_order);
        return !!success;
    }

    static BOOST_FORCEINLINE bool compare_exchange_strong(
        storage_type volatile& storage, storage_type& expected, storage_type desired, memory_order success_order, memory_order failure_order) BOOST_NOEXCEPT
    {
        int success;
        storage_type current, tmp;
        fence_before(success_order);
        __asm__ __volatile__
        (
            "1:\n"
            "mov %5, %1\n"                  // tmp = desired
            "ldl_l %2, %4\n"                // current = *(&storage)
            "cmpeq %2, %0, %3\n"            // success = current == expected
            "mov %2, %0\n"                  // expected = current
            "beq %3, 2f\n"                  // if (success == 0) goto end
            "stl_c %1, %4\n"                // storage = tmp; tmp = store succeeded
            "beq %1, 3f\n"                  // if (tmp == 0) goto retry
            "mov %1, %3\n"                  // success = tmp
            "2:\n"

            ".subsection 2\n"
            "3: br 1b\n"
            ".previous\n"

            : "+&r" (expected),  // %0
              "=&r" (tmp),       // %1
              "=&r" (current),   // %2
              "=&r" (success)    // %3
            : "m" (storage),     // %4
              "r" (desired)      // %5
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC
        );
        if (success)
            fence_after(success_order);
        else
            fence_after(failure_order);
        return !!success;
    }

    static BOOST_FORCEINLINE storage_type fetch_add(storage_type volatile& storage, storage_type v, memory_order order) BOOST_NOEXCEPT
    {
        storage_type original, modified;
        fence_before(order);
        __asm__ __volatile__
        (
            "1:\n"
            "ldl_l %0, %2\n"
            "addl %0, %3, %1\n"
            "stl_c %1, %2\n"
            "beq %1, 2f\n"

            ".subsection 2\n"
            "2: br 1b\n"
            ".previous\n"

            : "=&r" (original),  // %0
              "=&r" (modified)   // %1
            : "m" (storage),     // %2
              "r" (v)            // %3
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC
        );
        fence_after(order);
        return original;
    }

    static BOOST_FORCEINLINE storage_type fetch_sub(storage_type volatile& storage, storage_type v, memory_order order) BOOST_NOEXCEPT
    {
        storage_type original, modified;
        fence_before(order);
        __asm__ __volatile__
        (
            "1:\n"
            "ldl_l %0, %2\n"
            "subl %0, %3, %1\n"
            "stl_c %1, %2\n"
            "beq %1, 2f\n"

            ".subsection 2\n"
            "2: br 1b\n"
            ".previous\n"

            : "=&r" (original),  // %0
              "=&r" (modified)   // %1
            : "m" (storage),     // %2
              "r" (v)            // %3
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC
        );
        fence_after(order);
        return original;
    }

    static BOOST_FORCEINLINE storage_type fetch_and(storage_type volatile& storage, storage_type v, memory_order order) BOOST_NOEXCEPT
    {
        storage_type original, modified;
        fence_before(order);
        __asm__ __volatile__
        (
            "1:\n"
            "ldl_l %0, %2\n"
            "and %0, %3, %1\n"
            "stl_c %1, %2\n"
            "beq %1, 2f\n"

            ".subsection 2\n"
            "2: br 1b\n"
            ".previous\n"

            : "=&r" (original),  // %0
              "=&r" (modified)   // %1
            : "m" (storage),     // %2
              "r" (v)            // %3
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC
        );
        fence_after(order);
        return original;
    }

    static BOOST_FORCEINLINE storage_type fetch_or(storage_type volatile& storage, storage_type v, memory_order order) BOOST_NOEXCEPT
    {
        storage_type original, modified;
        fence_before(order);
        __asm__ __volatile__
        (
            "1:\n"
            "ldl_l %0, %2\n"
            "bis %0, %3, %1\n"
            "stl_c %1, %2\n"
            "beq %1, 2f\n"

            ".subsection 2\n"
            "2: br 1b\n"
            ".previous\n"

            : "=&r" (original),  // %0
              "=&r" (modified)   // %1
            : "m" (storage),     // %2
              "r" (v)            // %3
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC
        );
        fence_after(order);
        return original;
    }

    static BOOST_FORCEINLINE storage_type fetch_xor(storage_type volatile& storage, storage_type v, memory_order order) BOOST_NOEXCEPT
    {
        storage_type original, modified;
        fence_before(order);
        __asm__ __volatile__
        (
            "1:\n"
            "ldl_l %0, %2\n"
            "xor %0, %3, %1\n"
            "stl_c %1, %2\n"
            "beq %1, 2f\n"

            ".subsection 2\n"
            "2: br 1b\n"
            ".previous\n"

            : "=&r" (original),  // %0
              "=&r" (modified)   // %1
            : "m" (storage),     // %2
              "r" (v)            // %3
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC
        );
        fence_after(order);
        return original;
    }

    static BOOST_FORCEINLINE bool test_and_set(storage_type volatile& storage, memory_order order) BOOST_NOEXCEPT
    {
        return !!exchange(storage, (storage_type)1, order);
    }

    static BOOST_FORCEINLINE void clear(storage_type volatile& storage, memory_order order) BOOST_NOEXCEPT
    {
        store(storage, 0, order);
    }
};


template< >
struct operations< 1u, false > :
    public operations< 4u, false >
{
    typedef operations< 4u, false > base_type;
    typedef base_type::storage_type storage_type;

    static BOOST_FORCEINLINE storage_type fetch_add(storage_type volatile& storage, storage_type v, memory_order order) BOOST_NOEXCEPT
    {
        storage_type original, modified;
        fence_before(order);
        __asm__ __volatile__
        (
            "1:\n"
            "ldl_l %0, %2\n"
            "addl %0, %3, %1\n"
            "zapnot %1, #1, %1\n"
            "stl_c %1, %2\n"
            "beq %1, 2f\n"

            ".subsection 2\n"
            "2: br 1b\n"
            ".previous\n"

            : "=&r" (original),  // %0
              "=&r" (modified)   // %1
            : "m" (storage),     // %2
              "r" (v)            // %3
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC
        );
        fence_after(order);
        return original;
    }

    static BOOST_FORCEINLINE storage_type fetch_sub(storage_type volatile& storage, storage_type v, memory_order order) BOOST_NOEXCEPT
    {
        storage_type original, modified;
        fence_before(order);
        __asm__ __volatile__
        (
            "1:\n"
            "ldl_l %0, %2\n"
            "subl %0, %3, %1\n"
            "zapnot %1, #1, %1\n"
            "stl_c %1, %2\n"
            "beq %1, 2f\n"

            ".subsection 2\n"
            "2: br 1b\n"
            ".previous\n"

            : "=&r" (original),  // %0
              "=&r" (modified)   // %1
            : "m" (storage),     // %2
              "r" (v)            // %3
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC
        );
        fence_after(order);
        return original;
    }
};

template< >
struct operations< 1u, true > :
    public operations< 4u, true >
{
    typedef operations< 4u, true > base_type;
    typedef base_type::storage_type storage_type;

    static BOOST_FORCEINLINE storage_type fetch_add(storage_type volatile& storage, storage_type v, memory_order order) BOOST_NOEXCEPT
    {
        storage_type original, modified;
        fence_before(order);
        __asm__ __volatile__
        (
            "1:\n"
            "ldl_l %0, %2\n"
            "addl %0, %3, %1\n"
            "sextb %1, %1\n"
            "stl_c %1, %2\n"
            "beq %1, 2f\n"

            ".subsection 2\n"
            "2: br 1b\n"
            ".previous\n"

            : "=&r" (original),  // %0
              "=&r" (modified)   // %1
            : "m" (storage),     // %2
              "r" (v)            // %3
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC
        );
        fence_after(order);
        return original;
    }

    static BOOST_FORCEINLINE storage_type fetch_sub(storage_type volatile& storage, storage_type v, memory_order order) BOOST_NOEXCEPT
    {
        storage_type original, modified;
        fence_before(order);
        __asm__ __volatile__
        (
            "1:\n"
            "ldl_l %0, %2\n"
            "subl %0, %3, %1\n"
            "sextb %1, %1\n"
            "stl_c %1, %2\n"
            "beq %1, 2f\n"

            ".subsection 2\n"
            "2: br 1b\n"
            ".previous\n"

            : "=&r" (original),  // %0
              "=&r" (modified)   // %1
            : "m" (storage),     // %2
              "r" (v)            // %3
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC
        );
        fence_after(order);
        return original;
    }
};


template< >
struct operations< 2u, false > :
    public operations< 4u, false >
{
    typedef operations< 4u, false > base_type;
    typedef base_type::storage_type storage_type;

    static BOOST_FORCEINLINE storage_type fetch_add(storage_type volatile& storage, storage_type v, memory_order order) BOOST_NOEXCEPT
    {
        storage_type original, modified;
        fence_before(order);
        __asm__ __volatile__
        (
            "1:\n"
            "ldl_l %0, %2\n"
            "addl %0, %3, %1\n"
            "zapnot %1, #3, %1\n"
            "stl_c %1, %2\n"
            "beq %1, 2f\n"

            ".subsection 2\n"
            "2: br 1b\n"
            ".previous\n"

            : "=&r" (original),  // %0
              "=&r" (modified)   // %1
            : "m" (storage),     // %2
              "r" (v)            // %3
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC
        );
        fence_after(order);
        return original;
    }

    static BOOST_FORCEINLINE storage_type fetch_sub(storage_type volatile& storage, storage_type v, memory_order order) BOOST_NOEXCEPT
    {
        storage_type original, modified;
        fence_before(order);
        __asm__ __volatile__
        (
            "1:\n"
            "ldl_l %0, %2\n"
            "subl %0, %3, %1\n"
            "zapnot %1, #3, %1\n"
            "stl_c %1, %2\n"
            "beq %1, 2f\n"

            ".subsection 2\n"
            "2: br 1b\n"
            ".previous\n"

            : "=&r" (original),  // %0
              "=&r" (modified)   // %1
            : "m" (storage),     // %2
              "r" (v)            // %3
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC
        );
        fence_after(order);
        return original;
    }
};

template< >
struct operations< 2u, true > :
    public operations< 4u, true >
{
    typedef operations< 4u, true > base_type;
    typedef base_type::storage_type storage_type;

    static BOOST_FORCEINLINE storage_type fetch_add(storage_type volatile& storage, storage_type v, memory_order order) BOOST_NOEXCEPT
    {
        storage_type original, modified;
        fence_before(order);
        __asm__ __volatile__
        (
            "1:\n"
            "ldl_l %0, %2\n"
            "addl %0, %3, %1\n"
            "sextw %1, %1\n"
            "stl_c %1, %2\n"
            "beq %1, 2f\n"

            ".subsection 2\n"
            "2: br 1b\n"
            ".previous\n"

            : "=&r" (original),  // %0
              "=&r" (modified)   // %1
            : "m" (storage),     // %2
              "r" (v)            // %3
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC
        );
        fence_after(order);
        return original;
    }

    static BOOST_FORCEINLINE storage_type fetch_sub(storage_type volatile& storage, storage_type v, memory_order order) BOOST_NOEXCEPT
    {
        storage_type original, modified;
        fence_before(order);
        __asm__ __volatile__
        (
            "1:\n"
            "ldl_l %0, %2\n"
            "subl %0, %3, %1\n"
            "sextw %1, %1\n"
            "stl_c %1, %2\n"
            "beq %1, 2f\n"

            ".subsection 2\n"
            "2: br 1b\n"
            ".previous\n"

            : "=&r" (original),  // %0
              "=&r" (modified)   // %1
            : "m" (storage),     // %2
              "r" (v)            // %3
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC
        );
        fence_after(order);
        return original;
    }
};


template< bool Signed >
struct operations< 8u, Signed > :
    public gcc_alpha_operations_base
{
    typedef typename make_storage_type< 8u >::type storage_type;
    typedef typename make_storage_type< 8u >::aligned aligned_storage_type;

    static BOOST_CONSTEXPR_OR_CONST std::size_t storage_size = 8u;
    static BOOST_CONSTEXPR_OR_CONST bool is_signed = Signed;

    static BOOST_FORCEINLINE void store(storage_type volatile& storage, storage_type v, memory_order order) BOOST_NOEXCEPT
    {
        fence_before(order);
        storage = v;
        fence_after_store(order);
    }

    static BOOST_FORCEINLINE storage_type load(storage_type const volatile& storage, memory_order order) BOOST_NOEXCEPT
    {
        storage_type v = storage;
        fence_after(order);
        return v;
    }

    static BOOST_FORCEINLINE storage_type exchange(storage_type volatile& storage, storage_type v, memory_order order) BOOST_NOEXCEPT
    {
        storage_type original, tmp;
        fence_before(order);
        __asm__ __volatile__
        (
            "1:\n"
            "mov %3, %1\n"
            "ldq_l %0, %2\n"
            "stq_c %1, %2\n"
            "beq %1, 2f\n"

            ".subsection 2\n"
            "2: br 1b\n"
            ".previous\n"

            : "=&r" (original),  // %0
              "=&r" (tmp)        // %1
            : "m" (storage),     // %2
              "r" (v)            // %3
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC
        );
        fence_after(order);
        return original;
    }

    static BOOST_FORCEINLINE bool compare_exchange_weak(
        storage_type volatile& storage, storage_type& expected, storage_type desired, memory_order success_order, memory_order failure_order) BOOST_NOEXCEPT
    {
        fence_before(success_order);
        int success;
        storage_type current;
        __asm__ __volatile__
        (
            "1:\n"
            "ldq_l %2, %4\n"                // current = *(&storage)
            "cmpeq %2, %0, %3\n"            // success = current == expected
            "mov %2, %0\n"                  // expected = current
            "beq %3, 2f\n"                  // if (success == 0) goto end
            "stq_c %1, %4\n"                // storage = desired; desired = store succeeded
            "mov %1, %3\n"                  // success = desired
            "2:\n"
            : "+&r" (expected),  // %0
              "+&r" (desired),   // %1
              "=&r" (current),   // %2
              "=&r" (success)    // %3
            : "m" (storage)      // %4
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC
        );
        if (success)
            fence_after(success_order);
        else
            fence_after(failure_order);
        return !!success;
    }

    static BOOST_FORCEINLINE bool compare_exchange_strong(
        storage_type volatile& storage, storage_type& expected, storage_type desired, memory_order success_order, memory_order failure_order) BOOST_NOEXCEPT
    {
        int success;
        storage_type current, tmp;
        fence_before(success_order);
        __asm__ __volatile__
        (
            "1:\n"
            "mov %5, %1\n"                  // tmp = desired
            "ldq_l %2, %4\n"                // current = *(&storage)
            "cmpeq %2, %0, %3\n"            // success = current == expected
            "mov %2, %0\n"                  // expected = current
            "beq %3, 2f\n"                  // if (success == 0) goto end
            "stq_c %1, %4\n"                // storage = tmp; tmp = store succeeded
            "beq %1, 3f\n"                  // if (tmp == 0) goto retry
            "mov %1, %3\n"                  // success = tmp
            "2:\n"

            ".subsection 2\n"
            "3: br 1b\n"
            ".previous\n"

            : "+&r" (expected),  // %0
              "=&r" (tmp),       // %1
              "=&r" (current),   // %2
              "=&r" (success)    // %3
            : "m" (storage),     // %4
              "r" (desired)      // %5
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC
        );
        if (success)
            fence_after(success_order);
        else
            fence_after(failure_order);
        return !!success;
    }

    static BOOST_FORCEINLINE storage_type fetch_add(storage_type volatile& storage, storage_type v, memory_order order) BOOST_NOEXCEPT
    {
        storage_type original, modified;
        fence_before(order);
        __asm__ __volatile__
        (
            "1:\n"
            "ldq_l %0, %2\n"
            "addq %0, %3, %1\n"
            "stq_c %1, %2\n"
            "beq %1, 2f\n"

            ".subsection 2\n"
            "2: br 1b\n"
            ".previous\n"

            : "=&r" (original),  // %0
              "=&r" (modified)   // %1
            : "m" (storage),     // %2
              "r" (v)            // %3
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC
        );
        fence_after(order);
        return original;
    }

    static BOOST_FORCEINLINE storage_type fetch_sub(storage_type volatile& storage, storage_type v, memory_order order) BOOST_NOEXCEPT
    {
        storage_type original, modified;
        fence_before(order);
        __asm__ __volatile__
        (
            "1:\n"
            "ldq_l %0, %2\n"
            "subq %0, %3, %1\n"
            "stq_c %1, %2\n"
            "beq %1, 2f\n"

            ".subsection 2\n"
            "2: br 1b\n"
            ".previous\n"

            : "=&r" (original),  // %0
              "=&r" (modified)   // %1
            : "m" (storage),     // %2
              "r" (v)            // %3
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC
        );
        fence_after(order);
        return original;
    }

    static BOOST_FORCEINLINE storage_type fetch_and(storage_type volatile& storage, storage_type v, memory_order order) BOOST_NOEXCEPT
    {
        storage_type original, modified;
        fence_before(order);
        __asm__ __volatile__
        (
            "1:\n"
            "ldq_l %0, %2\n"
            "and %0, %3, %1\n"
            "stq_c %1, %2\n"
            "beq %1, 2f\n"

            ".subsection 2\n"
            "2: br 1b\n"
            ".previous\n"

            : "=&r" (original),  // %0
              "=&r" (modified)   // %1
            : "m" (storage),     // %2
              "r" (v)            // %3
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC
        );
        fence_after(order);
        return original;
    }

    static BOOST_FORCEINLINE storage_type fetch_or(storage_type volatile& storage, storage_type v, memory_order order) BOOST_NOEXCEPT
    {
        storage_type original, modified;
        fence_before(order);
        __asm__ __volatile__
        (
            "1:\n"
            "ldq_l %0, %2\n"
            "bis %0, %3, %1\n"
            "stq_c %1, %2\n"
            "beq %1, 2f\n"

            ".subsection 2\n"
            "2: br 1b\n"
            ".previous\n"

            : "=&r" (original),  // %0
              "=&r" (modified)   // %1
            : "m" (storage),     // %2
              "r" (v)            // %3
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC
        );
        fence_after(order);
        return original;
    }

    static BOOST_FORCEINLINE storage_type fetch_xor(storage_type volatile& storage, storage_type v, memory_order order) BOOST_NOEXCEPT
    {
        storage_type original, modified;
        fence_before(order);
        __asm__ __volatile__
        (
            "1:\n"
            "ldq_l %0, %2\n"
            "xor %0, %3, %1\n"
            "stq_c %1, %2\n"
            "beq %1, 2f\n"

            ".subsection 2\n"
            "2: br 1b\n"
            ".previous\n"

            : "=&r" (original),  // %0
              "=&r" (modified)   // %1
            : "m" (storage),     // %2
              "r" (v)            // %3
            : BOOST_ATOMIC_DETAIL_ASM_CLOBBER_CC
        );
        fence_after(order);
        return original;
    }

    static BOOST_FORCEINLINE bool test_and_set(storage_type volatile& storage, memory_order order) BOOST_NOEXCEPT
    {
        return !!exchange(storage, (storage_type)1, order);
    }

    static BOOST_FORCEINLINE void clear(storage_type volatile& storage, memory_order order) BOOST_NOEXCEPT
    {
        store(storage, 0, order);
    }
};


BOOST_FORCEINLINE void thread_fence(memory_order order) BOOST_NOEXCEPT
{
    if (order != memory_order_relaxed)
        __asm__ __volatile__ ("mb" ::: "memory");
}

BOOST_FORCEINLINE void signal_fence(memory_order order) BOOST_NOEXCEPT
{
    if (order != memory_order_relaxed)
        __asm__ __volatile__ ("" ::: "memory");
}

} // namespace detail
} // namespace atomics
} // namespace boost

#endif // BOOST_ATOMIC_DETAIL_OPS_GCC_ALPHA_HPP_INCLUDED_
