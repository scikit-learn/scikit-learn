#ifndef PYTHONIC_NUMPY_AROUND_HPP
#define PYTHONIC_NUMPY_AROUND_HPP

#include "pythonic/include/numpy/around.hpp"

#include "pythonic/numpy/asarray.hpp"
#include "pythonic/numpy/float64.hpp"
#include "pythonic/numpy/floor_divide.hpp"
#include "pythonic/numpy/multiply.hpp"
#include "pythonic/numpy/power.hpp"
#include "pythonic/numpy/rint.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  // fast path
  template <class E>
  auto around(E &&a) -> decltype(functor::rint{}(std::forward<E>(a)))
  {
    return functor::rint{}(std::forward<E>(a));
  }

  // generic floating point version, pure numpy_expr
  template <class E>
  auto around(E &&a, long decimals) -> std::enable_if_t<
      !std::is_integral<typename types::dtype_of<std::decay_t<E>>::type>::value,
      decltype(functor::rint{}(functor::multiply{}(
                   std::forward<E>(a),
                   std::declval<typename types::dtype_of<std::decay_t<E>>::type>())) /
               std::declval<typename types::dtype_of<std::decay_t<E>>::type>())>
  {
    typename types::dtype_of<std::decay_t<E>>::type const fact = functor::power{}(10., decimals);
    return functor::rint{}(functor::multiply{}(std::forward<E>(a), fact)) / fact;
  }

  // the integer version is only relevant when decimals < 0
  template <class E>
  auto around(E &&a, long decimals) -> std::enable_if_t<
      std::is_integral<typename types::dtype_of<std::decay_t<E>>::type>::value,
      decltype(numpy::functor::floor_divide{}(
                   functor::float64{}(std::forward<E>(a)),
                   std::declval<typename types::dtype_of<std::decay_t<E>>::type>()) *
               std::declval<typename types::dtype_of<std::decay_t<E>>::type>())>
  {
    typename types::dtype_of<std::decay_t<E>>::type const fact =
        functor::power{}(10L, std::max(0L, -decimals));
    return pythonic::numpy::functor::floor_divide{}(functor::float64{}(std::forward<E>(a)), fact) *
           fact;
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
