#ifndef PYTHONIC_INCLUDE_NUMPY_RAD2DEG_HPP
#define PYTHONIC_INCLUDE_NUMPY_RAD2DEG_HPP

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
    auto rad2deg(T const &val) -> decltype(val * 180 / pi)
    {
      return val * 180 / pi;
    }
  } // namespace wrapper
#define NUMPY_NARY_FUNC_NAME rad2deg
#define NUMPY_NARY_FUNC_SYM wrapper::rad2deg
#include "pythonic/include/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
