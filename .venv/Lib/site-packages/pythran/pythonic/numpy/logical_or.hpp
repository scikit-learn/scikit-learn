#ifndef PYTHONIC_NUMPY_LOGICALOR_HPP
#define PYTHONIC_NUMPY_LOGICALOR_HPP

#include "pythonic/include/numpy/logical_or.hpp"

#include "pythonic/types/ndarray.hpp"
#include "pythonic/types/numpy_broadcast.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace wrapper
  {
    template <class T0, class T1>
    auto logical_or(T0 const &t0, T1 const &t1) -> decltype(t0 || t1)
    {
      return t0 || t1;
    }
  } // namespace wrapper

#define NUMPY_NARY_FUNC_NAME logical_or
#define NUMPY_NARY_FUNC_SYM wrapper::logical_or
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
