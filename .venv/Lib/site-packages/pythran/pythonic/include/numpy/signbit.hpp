#ifndef PYTHONIC_INCLUDE_NUMPY_SIGNBIT_HPP
#define PYTHONIC_INCLUDE_NUMPY_SIGNBIT_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"

#include <xsimd/xsimd.hpp>

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace details
  {
    template <class T>
    auto signbit(T val) -> decltype(xsimd::signbit(val) == T(1))
    {
      return xsimd::signbit(val) == T(1);
    }
  } // namespace details

#define NUMPY_NARY_FUNC_NAME signbit
#define NUMPY_NARY_FUNC_SYM details::signbit
#include "pythonic/include/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
