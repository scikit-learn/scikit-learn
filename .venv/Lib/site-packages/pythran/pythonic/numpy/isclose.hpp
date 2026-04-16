#ifndef PYTHONIC_NUMPY_ISCLOSE_HPP
#define PYTHONIC_NUMPY_ISCLOSE_HPP

#include "pythonic/include/numpy/isclose.hpp"

#include "pythonic/numpy/abs.hpp"
#include "pythonic/numpy/isfinite.hpp"
#include "pythonic/numpy/isnan.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN
namespace numpy
{

  namespace wrapper
  {
    template <class T0, class T1>
    bool isclose(T0 const &u, T1 const &v, double rtol, double atol, bool equal_nan)
    {
      if (functor::isfinite()(u) && functor::isfinite()(v))
        return functor::abs()(u - v) <= (atol + rtol * functor::abs()(v));
      else if (functor::isnan()(u) && functor::isnan()(v))
        return equal_nan;
      else
        return (u == v);
    }
  } // namespace wrapper
#define NUMPY_NARY_FUNC_NAME isclose
#define NUMPY_NARY_FUNC_SYM wrapper::isclose
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
