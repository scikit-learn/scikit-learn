#ifndef BOOST_SMART_PTR_DETAIL_YIELD_K_HPP_INCLUDED
#define BOOST_SMART_PTR_DETAIL_YIELD_K_HPP_INCLUDED

// MS compatible compilers support #pragma once

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

//
//  yield_k.hpp
//
//  Copyright (c) 2008 Peter Dimov
//  Copyright (c) Microsoft Corporation 2014
//
//  void yield( unsigned k );
//
//  Typical use:
//
//  for( unsigned k = 0; !try_lock(); ++k ) yield( k );
//
//  Distributed under the Boost Software License, Version 1.0.
//  See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt
//

#include <boost/config.hpp>
#include <boost/predef/platform/windows_runtime.h>

#if BOOST_PLAT_WINDOWS_RUNTIME
#include <thread>
#endif

// BOOST_SMT_PAUSE

#if defined(_MSC_VER) && _MSC_VER >= 1310 && ( defined(_M_IX86) || defined(_M_X64) ) && !defined(__c2__)

extern "C" void _mm_pause();

#define BOOST_SMT_PAUSE _mm_pause();

#elif defined(__GNUC__) && ( defined(__i386__) || defined(__x86_64__) )

#define BOOST_SMT_PAUSE __asm__ __volatile__( "rep; nop" : : : "memory" );

#endif

//

#if defined( WIN32 ) || defined( _WIN32 ) || defined( __WIN32__ ) || defined( __CYGWIN__ )

#if defined( BOOST_USE_WINDOWS_H )
# include <windows.h>
#endif

namespace boost
{

namespace detail
{

#if !defined( BOOST_USE_WINDOWS_H ) && !BOOST_PLAT_WINDOWS_RUNTIME

#if defined(__clang__) && defined(__x86_64__)
// clang x64 warns that __stdcall is ignored
# define BOOST_SP_STDCALL
#else
# define BOOST_SP_STDCALL __stdcall
#endif

#if defined(__LP64__) // Cygwin 64
  extern "C" __declspec(dllimport) void BOOST_SP_STDCALL Sleep( unsigned int ms );
#else
  extern "C" __declspec(dllimport) void BOOST_SP_STDCALL Sleep( unsigned long ms );
#endif

#undef BOOST_SP_STDCALL

#endif // !defined( BOOST_USE_WINDOWS_H ) && !BOOST_PLAT_WINDOWS_RUNTIME

inline void yield( unsigned k ) BOOST_NOEXCEPT
{
    if( k < 4 )
    {
    }
#if defined( BOOST_SMT_PAUSE )
    else if( k < 16 )
    {
        BOOST_SMT_PAUSE
    }
#endif
#if !BOOST_PLAT_WINDOWS_RUNTIME
    else if( k < 32 )
    {
        Sleep( 0 );
    }
    else
    {
        Sleep( 1 );
    }
#else
    else
    {
        // Sleep isn't supported on the Windows Runtime.
        std::this_thread::yield();
    }
#endif
}

} // namespace detail

} // namespace boost

#elif defined( BOOST_HAS_PTHREADS )

#ifndef _AIX
#include <sched.h>
#else
   // AIX's sched.h defines ::var which sometimes conflicts with Lambda's var
       extern "C" int sched_yield(void);
#endif

#include <time.h>

namespace boost
{

namespace detail
{

inline void yield( unsigned k )
{
    if( k < 4 )
    {
    }
#if defined( BOOST_SMT_PAUSE )
    else if( k < 16 )
    {
        BOOST_SMT_PAUSE
    }
#endif
    else if( k < 32 || k & 1 )
    {
        sched_yield();
    }
    else
    {
        // g++ -Wextra warns on {} or {0}
        struct timespec rqtp = { 0, 0 };

        // POSIX says that timespec has tv_sec and tv_nsec
        // But it doesn't guarantee order or placement

        rqtp.tv_sec = 0;
        rqtp.tv_nsec = 1000;

        nanosleep( &rqtp, 0 );
    }
}

} // namespace detail

} // namespace boost

#else

namespace boost
{

namespace detail
{

inline void yield( unsigned )
{
}

} // namespace detail

} // namespace boost

#endif

#endif // #ifndef BOOST_SMART_PTR_DETAIL_YIELD_K_HPP_INCLUDED
