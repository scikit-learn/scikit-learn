#ifndef PYTHONIC_INCLUDE_NUMPY_ISPOSINF_HPP
#define PYTHONIC_INCLUDE_NUMPY_ISPOSINF_HPP

#include "pythonic/include/numpy/isinf.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace wrapper
  {
    template <class T>
    auto isposinf(T const &t) -> decltype(functor::isinf{}(t) && t >= 0);
  }
#define NUMPY_NARY_FUNC_NAME isposinf
#define NUMPY_NARY_FUNC_SYM wrapper::isposinf
#include "pythonic/include/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
