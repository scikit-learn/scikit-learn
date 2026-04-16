#ifndef PYTHONIC_INCLUDE_NUMPY_ISNEGINF_HPP
#define PYTHONIC_INCLUDE_NUMPY_ISNEGINF_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"

#include "pythonic/include/numpy/isinf.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace wrapper
  {
    template <class T>
    auto isneginf(T const &t) -> decltype(functor::isinf{}(t) && (t < 0));
  }

#define NUMPY_NARY_FUNC_NAME isneginf
#define NUMPY_NARY_FUNC_SYM wrapper::isneginf
#include "pythonic/include/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
