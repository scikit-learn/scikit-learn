#ifndef PYTHONIC_INCLUDE_NUMPY_RECIPROCAL_HPP
#define PYTHONIC_INCLUDE_NUMPY_RECIPROCAL_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace wrapper
  {
    template <class T>
    auto reciprocal(T const &val) -> decltype(static_cast<T>(1.) / val)
    {
      return static_cast<T>(1.) / val;
    }
  } // namespace wrapper
#define NUMPY_NARY_FUNC_NAME reciprocal
#define NUMPY_NARY_FUNC_SYM wrapper::reciprocal
#include "pythonic/include/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
