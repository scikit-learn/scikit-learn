#ifndef PYTHONIC_INCLUDE_NUMPY_AROUND_HPP
#define PYTHONIC_INCLUDE_NUMPY_AROUND_HPP

#include "pythonic/include/numpy/asarray.hpp"
#include "pythonic/include/numpy/float64.hpp"
#include "pythonic/include/numpy/floor_divide.hpp"
#include "pythonic/include/numpy/multiply.hpp"
#include "pythonic/include/numpy/rint.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  // fast path
  template <class E>
  auto around(E &&a) -> decltype(functor::rint{}(std::forward<E>(a)));

  // generic floating point version, pure numpy_expr
  template <class E>
  auto around(E &&a, long decimals) -> std::enable_if_t<
      !std::is_integral<typename types::dtype_of<std::decay_t<E>>::type>::value,
      decltype(functor::rint{}(functor::multiply{}(
                   std::forward<E>(a),
                   std::declval<typename types::dtype_of<std::decay_t<E>>::type>())) /
               std::declval<typename types::dtype_of<std::decay_t<E>>::type>())>;

  // the integer version is only relevant when decimals < 0
  template <class E>
  auto around(E &&a, long decimals) -> std::enable_if_t<
      std::is_integral<typename types::dtype_of<std::decay_t<E>>::type>::value,
      decltype(numpy::functor::floor_divide{}(
                   functor::float64{}(std::forward<E>(a)),
                   std::declval<typename types::dtype_of<std::decay_t<E>>::type>()) *
               std::declval<typename types::dtype_of<std::decay_t<E>>::type>())>;

  DEFINE_FUNCTOR(pythonic::numpy, around);
} // namespace numpy
PYTHONIC_NS_END

#endif
