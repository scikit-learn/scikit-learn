#ifndef PYTHONIC_INCLUDE_NUMPY_DEG2RAD_HPP
#define PYTHONIC_INCLUDE_NUMPY_DEG2RAD_HPP

#include "pythonic/include/numpy/pi.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace wrapper
  {
    template <class T>
    auto deg2rad(T const &val) -> decltype(val * pi / 180)
    {
      return val * pi / 180;
    }
  } // namespace wrapper
#define NUMPY_NARY_FUNC_NAME deg2rad
#define NUMPY_NARY_FUNC_SYM wrapper::deg2rad
#include "pythonic/include/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
