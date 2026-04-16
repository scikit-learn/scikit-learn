#ifndef PYTHONIC_INCLUDE_NUMPY_LOGADDEXP_HPP
#define PYTHONIC_INCLUDE_NUMPY_LOGADDEXP_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"

#include "pythonic/include/numpy/exp.hpp"
#include "pythonic/include/numpy/log.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace wrapper
  {
    template <class T0, class T1>
    auto logaddexp(T0 const &t0, T1 const &t1)
        -> decltype(functor::log{}(functor::exp{}(t0) + functor::exp{}(t1)));
  }

#define NUMPY_NARY_FUNC_NAME logaddexp
#define NUMPY_NARY_FUNC_SYM wrapper::logaddexp
#include "pythonic/include/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
