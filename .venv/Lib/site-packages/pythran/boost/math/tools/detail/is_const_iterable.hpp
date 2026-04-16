//  (C) Copyright John Maddock 2018.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_TOOLS_IS_CONST_ITERABLE_HPP
#define BOOST_MATH_TOOLS_IS_CONST_ITERABLE_HPP

#if !defined(BOOST_NO_CXX14_VARIABLE_TEMPLATES) && !defined(BOOST_NO_CXX11_DECLTYPE) && !defined(BOOST_NO_CXX11_SFINAE_EXPR)

#define BOOST_MATH_HAS_IS_CONST_ITERABLE

#include <boost/type_traits/is_detected.hpp>
#include <utility>

namespace boost {
   namespace math {
      namespace tools {
         namespace detail {

            template<class T>
            using begin_t = decltype(std::declval<const T&>().begin());
            template<class T>
            using end_t = decltype(std::declval<const T&>().end());
            template<class T>
            using const_iterator_t = typename T::const_iterator;

            template <class T>
            struct is_const_iterable
               : public boost::integral_constant<bool,
               boost::is_detected<begin_t, T>::value
               && boost::is_detected<end_t, T>::value
               && boost::is_detected<const_iterator_t, T>::value
               > {};

} } } }

#endif

#endif // BOOST_MATH_TOOLS_IS_CONST_ITERABLE_HPP
