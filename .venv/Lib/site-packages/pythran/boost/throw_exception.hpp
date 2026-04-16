#ifndef UUID_AA15E74A856F11E08B8D93F24824019B
#define UUID_AA15E74A856F11E08B8D93F24824019B

// MS compatible compilers support #pragma once

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

//
//  boost/throw_exception.hpp
//
//  Copyright (c) 2002 Peter Dimov and Multi Media Ltd.
//  Copyright (c) 2008-2009 Emil Dotchevski and Reverge Studios, Inc.
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)
//
//  http://www.boost.org/libs/utility/throw_exception.html
//

#include <boost/config.hpp>
#include <boost/detail/workaround.hpp>
#include <exception>

#if !defined( BOOST_EXCEPTION_DISABLE ) && defined( __BORLANDC__ ) && BOOST_WORKAROUND( __BORLANDC__, BOOST_TESTED_AT(0x593) )
# define BOOST_EXCEPTION_DISABLE
#endif

#if !defined( BOOST_EXCEPTION_DISABLE ) && defined( BOOST_MSVC ) && BOOST_WORKAROUND( BOOST_MSVC, < 1310 )
# define BOOST_EXCEPTION_DISABLE
#endif

#if !defined( BOOST_EXCEPTION_DISABLE )
# include <boost/exception/exception.hpp>
#if !defined(BOOST_THROW_EXCEPTION_CURRENT_FUNCTION)
# include <boost/current_function.hpp>
# define BOOST_THROW_EXCEPTION_CURRENT_FUNCTION BOOST_CURRENT_FUNCTION
#endif
# define BOOST_THROW_EXCEPTION(x) ::boost::exception_detail::throw_exception_(x,BOOST_THROW_EXCEPTION_CURRENT_FUNCTION,__FILE__,__LINE__)
#else
# define BOOST_THROW_EXCEPTION(x) ::boost::throw_exception(x)
#endif

#if defined(__GNUC__) && (__GNUC__*100+__GNUC_MINOR__>301) && !defined(BOOST_EXCEPTION_ENABLE_WARNINGS)
#pragma GCC system_header
#endif
#if defined(_MSC_VER) && !defined(BOOST_EXCEPTION_ENABLE_WARNINGS)
#pragma warning(push,1)
#endif

namespace boost
{
#ifdef BOOST_NO_EXCEPTIONS

BOOST_NORETURN void throw_exception( std::exception const & e ); // user defined

#else

inline void throw_exception_assert_compatibility( std::exception const & ) { }

template<class E> BOOST_NORETURN inline void throw_exception( E const & e )
{
    //All boost exceptions are required to derive from std::exception,
    //to ensure compatibility with BOOST_NO_EXCEPTIONS.
    throw_exception_assert_compatibility(e);

#ifndef BOOST_EXCEPTION_DISABLE
    throw exception_detail::enable_both( e );
#else
    throw e;
#endif
}

#endif

#if !defined( BOOST_EXCEPTION_DISABLE )
    namespace
    exception_detail
    {
        template <class E>
        BOOST_NORETURN
        void
        throw_exception_( E const & x, char const * current_function, char const * file, int line )
        {
            boost::throw_exception(
                set_info(
                    set_info(
                        set_info(
                            enable_error_info(x),
                            throw_function(current_function)),
                        throw_file(file)),
                    throw_line(line)));
        }
    }
#endif
} // namespace boost

#if defined(_MSC_VER) && !defined(BOOST_EXCEPTION_ENABLE_WARNINGS)
#pragma warning(pop)
#endif
#endif
