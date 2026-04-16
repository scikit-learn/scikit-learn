#ifndef PYTHONIC_NUMPY_ANGLE_HPP
#define PYTHONIC_NUMPY_ANGLE_HPP

#include "pythonic/include/numpy/angle.hpp"

#include "pythonic/numpy/angle_in_deg.hpp"
#include "pythonic/numpy/angle_in_rad.hpp"
#include "pythonic/types/assignable.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class T>
  auto angle(T const &t, bool in_deg) ->
      typename assignable<decltype(functor::angle_in_rad()(t))>::type
  // assignable to find a common type between the two expression templates
  {
    if (in_deg)
      return functor::angle_in_deg()(t);
    else
      return functor::angle_in_rad()(t);
  }

  // Numpy_expr can be use if only the first argument is given.
  template <class T>
  auto angle(T const &t) -> decltype(functor::angle_in_rad()(t))
  {
    return functor::angle_in_rad()(t);
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
