#ifndef PYTHONIC_INCLUDE_NUMPY_REMAINDER_HPP
#define PYTHONIC_INCLUDE_NUMPY_REMAINDER_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/types/numpy_broadcast.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"

#include <xsimd/xsimd.hpp>

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace wrapper
  {
    template <class T0, class T1>
    auto remainder(T0 const &x, T1 const &y) -> decltype(x - y * xsimd::floor(x / y))
    {
      return x - y * xsimd::floor(x / y);
    }
  } // namespace wrapper

#define NUMPY_NARY_FUNC_NAME remainder
#define NUMPY_NARY_FUNC_SYM wrapper::remainder
#include "pythonic/include/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
