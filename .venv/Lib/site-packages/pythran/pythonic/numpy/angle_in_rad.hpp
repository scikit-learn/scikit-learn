#ifndef PYTHONIC_NUMPY_ANGLEINRAD_HPP
#define PYTHONIC_NUMPY_ANGLEINRAD_HPP

#include "pythonic/include/numpy/angle_in_rad.hpp"

#include "pythonic/numpy/arctan.hpp"
#include "pythonic/numpy/pi.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/numpy_traits.hpp"

/* NOTE: angle_in_rad is not part of the official Numpy API,
 * this file is here only to split the angle function in two parts
 */

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace wrapper
  {
    template <class T>
    auto angle_in_rad(T const &t) -> decltype(std::atan2(std::imag(t), std::real(t)))
    {
      return std::atan2(std::imag(t), std::real(t));
    }
  } // namespace wrapper

#define NUMPY_NARY_FUNC_NAME angle_in_rad
#define NUMPY_NARY_FUNC_SYM wrapper::angle_in_rad
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
