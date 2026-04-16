/*
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 * Copyright (c) 2009 Helge Bahmann
 * Copyright (c) 2012 Tim Blechmann
 * Copyright (c) 2013 - 2014 Andrey Semashev
 */
/*!
 * \file   atomic/detail/storage_type.hpp
 *
 * This header defines underlying types used as storage
 */

#ifndef BOOST_ATOMIC_DETAIL_STORAGE_TYPE_HPP_INCLUDED_
#define BOOST_ATOMIC_DETAIL_STORAGE_TYPE_HPP_INCLUDED_

#include <cstddef>
#include <boost/cstdint.hpp>
#include <boost/atomic/detail/config.hpp>
#include <boost/atomic/detail/string_ops.hpp>

#ifdef BOOST_HAS_PRAGMA_ONCE
#pragma once
#endif

namespace boost {
namespace atomics {
namespace detail {

template< typename T >
BOOST_FORCEINLINE void non_atomic_load(T const volatile& from, T& to) BOOST_NOEXCEPT
{
    to = from;
}

template< std::size_t Size >
struct BOOST_ATOMIC_DETAIL_MAY_ALIAS buffer_storage
{
    BOOST_ALIGNMENT(16) unsigned char data[Size];

    BOOST_FORCEINLINE bool operator! () const BOOST_NOEXCEPT
    {
        return (data[0] == 0u && BOOST_ATOMIC_DETAIL_MEMCMP(data, data + 1, Size - 1) == 0);
    }

    BOOST_FORCEINLINE bool operator== (buffer_storage const& that) const BOOST_NOEXCEPT
    {
        return BOOST_ATOMIC_DETAIL_MEMCMP(data, that.data, Size) == 0;
    }

    BOOST_FORCEINLINE bool operator!= (buffer_storage const& that) const BOOST_NOEXCEPT
    {
        return BOOST_ATOMIC_DETAIL_MEMCMP(data, that.data, Size) != 0;
    }
};

template< std::size_t Size >
BOOST_FORCEINLINE void non_atomic_load(buffer_storage< Size > const volatile& from, buffer_storage< Size >& to) BOOST_NOEXCEPT
{
    BOOST_ATOMIC_DETAIL_MEMCPY(to.data, const_cast< unsigned char const* >(from.data), Size);
}

template< std::size_t Size >
struct make_storage_type
{
    typedef buffer_storage< Size > type;

    struct BOOST_ATOMIC_DETAIL_MAY_ALIAS aligned
    {
        type value;

        BOOST_DEFAULTED_FUNCTION(aligned(), {})
        BOOST_FORCEINLINE BOOST_CONSTEXPR explicit aligned(type const& v) BOOST_NOEXCEPT : value(v) {}
    };
};

template< >
struct make_storage_type< 1u >
{
    typedef boost::uint8_t BOOST_ATOMIC_DETAIL_MAY_ALIAS type;

    struct BOOST_ATOMIC_DETAIL_MAY_ALIAS aligned
    {
        type value;

        BOOST_DEFAULTED_FUNCTION(aligned(), {})
        BOOST_FORCEINLINE BOOST_CONSTEXPR explicit aligned(type v) BOOST_NOEXCEPT : value(v) {}
    };
};

template< >
struct make_storage_type< 2u >
{
    typedef boost::uint16_t BOOST_ATOMIC_DETAIL_MAY_ALIAS type;

    struct BOOST_ATOMIC_DETAIL_MAY_ALIAS aligned
    {
        BOOST_ALIGNMENT(2) type value;

        BOOST_DEFAULTED_FUNCTION(aligned(), {})
        BOOST_FORCEINLINE BOOST_CONSTEXPR explicit aligned(type v) BOOST_NOEXCEPT : value(v) {}
    };
};

template< >
struct make_storage_type< 4u >
{
    typedef boost::uint32_t BOOST_ATOMIC_DETAIL_MAY_ALIAS type;

    struct BOOST_ATOMIC_DETAIL_MAY_ALIAS aligned
    {
        BOOST_ALIGNMENT(4) type value;

        BOOST_DEFAULTED_FUNCTION(aligned(), {})
        BOOST_FORCEINLINE BOOST_CONSTEXPR explicit aligned(type v) BOOST_NOEXCEPT : value(v) {}
    };
};

template< >
struct make_storage_type< 8u >
{
    typedef boost::uint64_t BOOST_ATOMIC_DETAIL_MAY_ALIAS type;

    struct BOOST_ATOMIC_DETAIL_MAY_ALIAS aligned
    {
        BOOST_ALIGNMENT(8) type value;

        BOOST_DEFAULTED_FUNCTION(aligned(), {})
        BOOST_FORCEINLINE BOOST_CONSTEXPR explicit aligned(type v) BOOST_NOEXCEPT : value(v) {}
    };
};

#if defined(BOOST_HAS_INT128)

template< >
struct make_storage_type< 16u >
{
    typedef boost::uint128_type BOOST_ATOMIC_DETAIL_MAY_ALIAS type;

    struct BOOST_ATOMIC_DETAIL_MAY_ALIAS aligned
    {
        BOOST_ALIGNMENT(16) type value;

        BOOST_DEFAULTED_FUNCTION(aligned(), {})
        BOOST_FORCEINLINE BOOST_CONSTEXPR explicit aligned(type v) BOOST_NOEXCEPT : value(v) {}
    };
};

#elif !defined(BOOST_NO_ALIGNMENT)

struct BOOST_ATOMIC_DETAIL_MAY_ALIAS storage128_t
{
    typedef boost::uint64_t BOOST_ATOMIC_DETAIL_MAY_ALIAS element_type;

    element_type data[2];

    BOOST_FORCEINLINE bool operator! () const BOOST_NOEXCEPT
    {
        return (data[0] | data[1]) == 0u;
    }
};

BOOST_FORCEINLINE bool operator== (storage128_t const& left, storage128_t const& right) BOOST_NOEXCEPT
{
    return ((left.data[0] ^ right.data[0]) | (left.data[1] ^ right.data[1])) == 0u;
}
BOOST_FORCEINLINE bool operator!= (storage128_t const& left, storage128_t const& right) BOOST_NOEXCEPT
{
    return !(left == right);
}

BOOST_FORCEINLINE void non_atomic_load(storage128_t const volatile& from, storage128_t& to) BOOST_NOEXCEPT
{
    to.data[0] = from.data[0];
    to.data[1] = from.data[1];
}

template< >
struct make_storage_type< 16u >
{
    typedef storage128_t type;

    struct BOOST_ATOMIC_DETAIL_MAY_ALIAS aligned
    {
        BOOST_ALIGNMENT(16) type value;

        BOOST_DEFAULTED_FUNCTION(aligned(), {})
        BOOST_FORCEINLINE BOOST_CONSTEXPR explicit aligned(type const& v) BOOST_NOEXCEPT : value(v) {}
    };
};

#endif

template< typename T >
struct storage_size_of
{
    static BOOST_CONSTEXPR_OR_CONST std::size_t size = sizeof(T);
    static BOOST_CONSTEXPR_OR_CONST std::size_t value = (size == 3u ? 4u : (size >= 5u && size <= 7u ? 8u : (size >= 9u && size <= 15u ? 16u : size)));
};

} // namespace detail
} // namespace atomics
} // namespace boost

#endif // BOOST_ATOMIC_DETAIL_STORAGE_TYPE_HPP_INCLUDED_
