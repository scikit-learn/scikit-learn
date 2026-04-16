#ifndef PYTHONIC_INCLUDE_NUMPY_ISCLOSE_HPP
#define PYTHONIC_INCLUDE_NUMPY_ISCLOSE_HPP

#include "pythonic/include/numpy/abs.hpp"
#include "pythonic/include/numpy/isfinite.hpp"
#include "pythonic/include/numpy/isnan.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN
namespace numpy
{

  namespace wrapper
  {
    template <class T0, class T1>
    bool isclose(T0 const &u, T1 const &v, double rtol = 1e-5, double atol = 1e-8,
                 bool equal_nan = false);
  }
#define NUMPY_NARY_FUNC_NAME isclose
#define NUMPY_NARY_FUNC_SYM wrapper::isclose
#include "pythonic/include/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
