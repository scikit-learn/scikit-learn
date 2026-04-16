/*
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 * Copyright (c) 2014 Andrey Semashev
 */
/*!
 * \file   atomic/detail/atomic_flag.hpp
 *
 * This header contains interface definition of \c atomic_flag.
 */

#ifndef BOOST_ATOMIC_DETAIL_ATOMIC_FLAG_HPP_INCLUDED_
#define BOOST_ATOMIC_DETAIL_ATOMIC_FLAG_HPP_INCLUDED_

#include <boost/assert.hpp>
#include <boost/memory_order.hpp>
#include <boost/atomic/detail/config.hpp>
#include <boost/atomic/detail/operations_lockfree.hpp>

#ifdef BOOST_HAS_PRAGMA_ONCE
#pragma once
#endif

/*
 * IMPLEMENTATION NOTE: All interface functions MUST be declared with BOOST_FORCEINLINE,
 *                      see comment for convert_memory_order_to_gcc in ops_gcc_atomic.hpp.
 */

namespace boost {
namespace atomics {

#if defined(BOOST_NO_CXX11_CONSTEXPR) || defined(BOOST_NO_CXX11_UNIFIED_INITIALIZATION_SYNTAX)
#define BOOST_ATOMIC_NO_ATOMIC_FLAG_INIT
#else
#define BOOST_ATOMIC_FLAG_INIT {}
#endif

struct atomic_flag
{
    typedef atomics::detail::operations< 1u, false > operations;
    typedef operations::storage_type storage_type;

    operations::aligned_storage_type m_storage;

    BOOST_FORCEINLINE BOOST_CONSTEXPR atomic_flag() BOOST_NOEXCEPT : m_storage(0)
    {
    }

    BOOST_FORCEINLINE bool test_and_set(memory_order order = memory_order_seq_cst) volatile BOOST_NOEXCEPT
    {
        return operations::test_and_set(m_storage.value, order);
    }

    BOOST_FORCEINLINE void clear(memory_order order = memory_order_seq_cst) volatile BOOST_NOEXCEPT
    {
        BOOST_ASSERT(order != memory_order_consume);
        BOOST_ASSERT(order != memory_order_acquire);
        BOOST_ASSERT(order != memory_order_acq_rel);
        operations::clear(m_storage.value, order);
    }

    BOOST_DELETED_FUNCTION(atomic_flag(atomic_flag const&))
    BOOST_DELETED_FUNCTION(atomic_flag& operator= (atomic_flag const&))
};

} // namespace atomics
} // namespace boost

#endif // BOOST_ATOMIC_DETAIL_ATOMIC_FLAG_HPP_INCLUDED_
