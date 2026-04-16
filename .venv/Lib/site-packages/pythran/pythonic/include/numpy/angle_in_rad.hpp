#ifndef PYTHONIC_INCLUDE_NUMPY_ANGLEINRAD_HPP
#define PYTHONIC_INCLUDE_NUMPY_ANGLEINRAD_HPP

#include "pythonic/include/numpy/arctan.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"

/* NOTE: angle_in_rad is not part of the official Numpy API,
 * this file is here only to split the angle function in two parts
 */

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace wrapper
  {
    template <class T>
    auto angle_in_rad(T const &t) -> decltype(std::atan2(std::imag(t), std::real(t)));
  }
#define NUMPY_NARY_FUNC_NAME angle_in_rad
#define NUMPY_NARY_FUNC_SYM wrapper::angle_in_rad
#include "pythonic/include/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
