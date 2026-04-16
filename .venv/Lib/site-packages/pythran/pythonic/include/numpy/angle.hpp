#ifndef PYTHONIC_INCLUDE_NUMPY_ANGLE_HPP
#define PYTHONIC_INCLUDE_NUMPY_ANGLE_HPP

#include "pythonic/include/numpy/angle_in_deg.hpp"
#include "pythonic/include/numpy/angle_in_rad.hpp"
#include "pythonic/include/types/assignable.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class T>
  auto angle(T const &t, bool in_deg) ->
      typename assignable<decltype(functor::angle_in_rad()(t))>::type;

  // Numpy_expr can be use if only the first argument is given.
  template <class T>
  auto angle(T const &t) -> decltype(functor::angle_in_rad()(t));

  DEFINE_FUNCTOR(pythonic::numpy, angle);
} // namespace numpy
PYTHONIC_NS_END

#endif
