#ifndef PYTHONIC_INCLUDE_NUMPY_HEAVISIDE_HPP
#define PYTHONIC_INCLUDE_NUMPY_HEAVISIDE_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace details
  {
    template <class T0, class T1>
    T1 heaviside(T0 x0, T1 x1);
  }
#define NUMPY_NARY_FUNC_NAME heaviside
#define NUMPY_NARY_FUNC_SYM details::heaviside
#include "pythonic/include/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
