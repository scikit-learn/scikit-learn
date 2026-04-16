/*
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 * Copyright (c) 2014 Andrey Semashev
 */
/*!
 * \file   atomic/detail/ops_extending_cas_based.hpp
 *
 * This header contains a boilerplate of the \c operations template implementation that requires sign/zero extension in arithmetic operations.
 */

#ifndef BOOST_ATOMIC_DETAIL_OPS_EXTENDING_CAS_BASED_HPP_INCLUDED_
#define BOOST_ATOMIC_DETAIL_OPS_EXTENDING_CAS_BASED_HPP_INCLUDED_

#include <cstddef>
#include <boost/memory_order.hpp>
#include <boost/atomic/detail/config.hpp>
#include <boost/atomic/detail/storage_type.hpp>
#include <boost/atomic/detail/integral_extend.hpp>

#ifdef BOOST_HAS_PRAGMA_ONCE
#pragma once
#endif

namespace boost {
namespace atomics {
namespace detail {

template< typename Base, std::size_t Size, bool Signed >
struct extending_cas_based_operations :
    public Base
{
    typedef typename Base::storage_type storage_type;
    typedef typename make_storage_type< Size >::type emulated_storage_type;

    static BOOST_FORCEINLINE storage_type fetch_add(storage_type volatile& storage, storage_type v, memory_order order) BOOST_NOEXCEPT
    {
        storage_type old_val;
        atomics::detail::non_atomic_load(storage, old_val);
        storage_type new_val;
        do
        {
            new_val = atomics::detail::integral_extend< Signed, storage_type >(static_cast< emulated_storage_type >(old_val + v));
        }
        while (!Base::compare_exchange_weak(storage, old_val, new_val, order, memory_order_relaxed));
        return old_val;
    }

    static BOOST_FORCEINLINE storage_type fetch_sub(storage_type volatile& storage, storage_type v, memory_order order) BOOST_NOEXCEPT
    {
        storage_type old_val;
        atomics::detail::non_atomic_load(storage, old_val);
        storage_type new_val;
        do
        {
            new_val = atomics::detail::integral_extend< Signed, storage_type >(static_cast< emulated_storage_type >(old_val - v));
        }
        while (!Base::compare_exchange_weak(storage, old_val, new_val, order, memory_order_relaxed));
        return old_val;
    }
};

} // namespace detail
} // namespace atomics
} // namespace boost

#endif // BOOST_ATOMIC_DETAIL_OPS_EXTENDING_CAS_BASED_HPP_INCLUDED_
