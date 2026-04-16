#ifndef PYTHONIC_NUMPY_ARGMIN_HPP
#define PYTHONIC_NUMPY_ARGMIN_HPP

#include "pythonic/include/numpy/argmin.hpp"

#include "pythonic/numpy/argminmax.hpp"
#include "pythonic/numpy/minimum.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class E>
  struct argmin_op {
    using op = functor::minimum;
    using expr_type = E;
    static typename E::dtype constexpr limit()
    {
      return std::numeric_limits<typename E::dtype>::max();
    }
    template <class T>
    static T elements(T first, T last)
    {
      return std::min_element(first, last);
    }
    template <class T>
    static bool value(T self, T other)
    {
      return self < other;
    }
  };

  template <class E>
  long argmin(E const &expr)
  {
    return argminmax<argmin_op<E>>(expr);
  }

  template <class E>
  types::ndarray<long, types::array_tuple<long, E::value - 1>> argmin(E const &expr, long axis)
  {
    return argminmax<argmin_op<E>>(expr, axis);
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
