//  (C) Copyright John Maddock 2005-2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_TOOLS_FRACTION_INCLUDED
#define BOOST_MATH_TOOLS_FRACTION_INCLUDED

#ifdef _MSC_VER
#pragma once
#endif

#include <boost/config/no_tr1/cmath.hpp>
#include <boost/cstdint.hpp>
#include <boost/type_traits/integral_constant.hpp>
#include <boost/mpl/if.hpp>
#include <boost/math/tools/precision.hpp>
#include <boost/math/tools/complex.hpp>

namespace boost{ namespace math{ namespace tools{

namespace detail
{

   template <class T>
   struct is_pair : public boost::false_type{};

   template <class T, class U>
   struct is_pair<std::pair<T,U> > : public boost::true_type{};

   template <class Gen>
   struct fraction_traits_simple
   {
       typedef typename Gen::result_type result_type;
       typedef typename Gen::result_type value_type;

       static result_type a(const value_type&) BOOST_MATH_NOEXCEPT(value_type)
       {
          return 1;
       }
       static result_type b(const value_type& v) BOOST_MATH_NOEXCEPT(value_type)
       {
          return v;
       }
   };

   template <class Gen>
   struct fraction_traits_pair
   {
       typedef typename Gen::result_type value_type;
       typedef typename value_type::first_type result_type;

       static result_type a(const value_type& v) BOOST_MATH_NOEXCEPT(value_type)
       {
          return v.first;
       }
       static result_type b(const value_type& v) BOOST_MATH_NOEXCEPT(value_type)
       {
          return v.second;
       }
   };

   template <class Gen>
   struct fraction_traits
       : public boost::mpl::if_c<
         is_pair<typename Gen::result_type>::value,
         fraction_traits_pair<Gen>,
         fraction_traits_simple<Gen> >::type
   {
   };

   template <class T, bool = is_complex_type<T>::value>
   struct tiny_value
   {
      static T get() {
         return tools::min_value<T>(); 
      }
   };
   template <class T>
   struct tiny_value<T, true>
   {
      typedef typename T::value_type value_type;
      static T get() {
         return tools::min_value<value_type>();
      }
   };

} // namespace detail

//
// continued_fraction_b
// Evaluates:
//
// b0 +       a1
//      ---------------
//      b1 +     a2
//           ----------
//           b2 +   a3
//                -----
//                b3 + ...
//
// Note that the first a0 returned by generator Gen is disarded.
//
template <class Gen, class U>
inline typename detail::fraction_traits<Gen>::result_type continued_fraction_b(Gen& g, const U& factor, boost::uintmax_t& max_terms) 
      BOOST_NOEXCEPT_IF(BOOST_MATH_IS_FLOAT(typename detail::fraction_traits<Gen>::result_type) && noexcept(std::declval<Gen>()()))
{
   BOOST_MATH_STD_USING // ADL of std names

   typedef detail::fraction_traits<Gen> traits;
   typedef typename traits::result_type result_type;
   typedef typename traits::value_type value_type;
   typedef typename integer_scalar_type<result_type>::type integer_type;
   typedef typename scalar_type<result_type>::type scalar_type;

   integer_type const zero(0), one(1);

   result_type tiny = detail::tiny_value<result_type>::get();
   scalar_type terminator = abs(factor);

   value_type v = g();

   result_type f, C, D, delta;
   f = traits::b(v);
   if(f == zero)
      f = tiny;
   C = f;
   D = 0;

   boost::uintmax_t counter(max_terms);

   do{
      v = g();
      D = traits::b(v) + traits::a(v) * D;
      if(D == result_type(0))
         D = tiny;
      C = traits::b(v) + traits::a(v) / C;
      if(C == zero)
         C = tiny;
      D = one/D;
      delta = C*D;
      f = f * delta;
   }while((abs(delta - one) > terminator) && --counter);

   max_terms = max_terms - counter;

   return f;
}

template <class Gen, class U>
inline typename detail::fraction_traits<Gen>::result_type continued_fraction_b(Gen& g, const U& factor)
   BOOST_NOEXCEPT_IF(BOOST_MATH_IS_FLOAT(typename detail::fraction_traits<Gen>::result_type) && noexcept(std::declval<Gen>()()))
{
   boost::uintmax_t max_terms = (std::numeric_limits<boost::uintmax_t>::max)();
   return continued_fraction_b(g, factor, max_terms);
}

template <class Gen>
inline typename detail::fraction_traits<Gen>::result_type continued_fraction_b(Gen& g, int bits)
   BOOST_NOEXCEPT_IF(BOOST_MATH_IS_FLOAT(typename detail::fraction_traits<Gen>::result_type) && noexcept(std::declval<Gen>()()))
{
   BOOST_MATH_STD_USING // ADL of std names

   typedef detail::fraction_traits<Gen> traits;
   typedef typename traits::result_type result_type;

   result_type factor = ldexp(1.0f, 1 - bits); // 1 / pow(result_type(2), bits);
   boost::uintmax_t max_terms = (std::numeric_limits<boost::uintmax_t>::max)();
   return continued_fraction_b(g, factor, max_terms);
}

template <class Gen>
inline typename detail::fraction_traits<Gen>::result_type continued_fraction_b(Gen& g, int bits, boost::uintmax_t& max_terms)
   BOOST_NOEXCEPT_IF(BOOST_MATH_IS_FLOAT(typename detail::fraction_traits<Gen>::result_type) && noexcept(std::declval<Gen>()()))
{
   BOOST_MATH_STD_USING // ADL of std names

   typedef detail::fraction_traits<Gen> traits;
   typedef typename traits::result_type result_type;

   result_type factor = ldexp(1.0f, 1 - bits); // 1 / pow(result_type(2), bits);
   return continued_fraction_b(g, factor, max_terms);
}

//
// continued_fraction_a
// Evaluates:
//
//            a1
//      ---------------
//      b1 +     a2
//           ----------
//           b2 +   a3
//                -----
//                b3 + ...
//
// Note that the first a1 and b1 returned by generator Gen are both used.
//
template <class Gen, class U>
inline typename detail::fraction_traits<Gen>::result_type continued_fraction_a(Gen& g, const U& factor, boost::uintmax_t& max_terms)
   BOOST_NOEXCEPT_IF(BOOST_MATH_IS_FLOAT(typename detail::fraction_traits<Gen>::result_type) && noexcept(std::declval<Gen>()()))
{
   BOOST_MATH_STD_USING // ADL of std names

   typedef detail::fraction_traits<Gen> traits;
   typedef typename traits::result_type result_type;
   typedef typename traits::value_type value_type;
   typedef typename integer_scalar_type<result_type>::type integer_type;
   typedef typename scalar_type<result_type>::type scalar_type;

   integer_type const zero(0), one(1);

   result_type tiny = detail::tiny_value<result_type>::get();
   scalar_type terminator = abs(factor);

   value_type v = g();

   result_type f, C, D, delta, a0;
   f = traits::b(v);
   a0 = traits::a(v);
   if(f == zero)
      f = tiny;
   C = f;
   D = 0;

   boost::uintmax_t counter(max_terms);

   do{
      v = g();
      D = traits::b(v) + traits::a(v) * D;
      if(D == zero)
         D = tiny;
      C = traits::b(v) + traits::a(v) / C;
      if(C == zero)
         C = tiny;
      D = one/D;
      delta = C*D;
      f = f * delta;
   }while((abs(delta - one) > terminator) && --counter);

   max_terms = max_terms - counter;

   return a0/f;
}

template <class Gen, class U>
inline typename detail::fraction_traits<Gen>::result_type continued_fraction_a(Gen& g, const U& factor)
   BOOST_NOEXCEPT_IF(BOOST_MATH_IS_FLOAT(typename detail::fraction_traits<Gen>::result_type) && noexcept(std::declval<Gen>()()))
{
   boost::uintmax_t max_iter = (std::numeric_limits<boost::uintmax_t>::max)();
   return continued_fraction_a(g, factor, max_iter);
}

template <class Gen>
inline typename detail::fraction_traits<Gen>::result_type continued_fraction_a(Gen& g, int bits)
   BOOST_NOEXCEPT_IF(BOOST_MATH_IS_FLOAT(typename detail::fraction_traits<Gen>::result_type) && noexcept(std::declval<Gen>()()))
{
   BOOST_MATH_STD_USING // ADL of std names

   typedef detail::fraction_traits<Gen> traits;
   typedef typename traits::result_type result_type;

   result_type factor = ldexp(1.0f, 1-bits); // 1 / pow(result_type(2), bits);
   boost::uintmax_t max_iter = (std::numeric_limits<boost::uintmax_t>::max)();

   return continued_fraction_a(g, factor, max_iter);
}

template <class Gen>
inline typename detail::fraction_traits<Gen>::result_type continued_fraction_a(Gen& g, int bits, boost::uintmax_t& max_terms)
   BOOST_NOEXCEPT_IF(BOOST_MATH_IS_FLOAT(typename detail::fraction_traits<Gen>::result_type) && noexcept(std::declval<Gen>()()))
{
   BOOST_MATH_STD_USING // ADL of std names

   typedef detail::fraction_traits<Gen> traits;
   typedef typename traits::result_type result_type;

   result_type factor = ldexp(1.0f, 1-bits); // 1 / pow(result_type(2), bits);
   return continued_fraction_a(g, factor, max_terms);
}

} // namespace tools
} // namespace math
} // namespace boost

#endif // BOOST_MATH_TOOLS_FRACTION_INCLUDED

