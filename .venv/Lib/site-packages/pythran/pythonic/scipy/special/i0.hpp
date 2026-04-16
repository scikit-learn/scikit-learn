#ifndef PYTHONIC_SCIPY_SPECIAL_I0_HPP
#define PYTHONIC_SCIPY_SPECIAL_I0_HPP

#include "pythonic/include/scipy/special/i0.hpp"
#include "pythonic/scipy/special/chbevl.hpp"

#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace scipy
{
  namespace special
  {
    namespace details
    {
      template <class T>
      double i0(T x_)
      {
        double y;
        double x = x_;

        if (x < 0)
          x = -x;
        if (x <= 8.0) {
          y = (x / 2.0) - 2.0;
          return (exp(x) * chbevl(y, A));
        }

        return (exp(x) * chbevl(32.0 / x - 2.0, B) / sqrt(x));
      }
    } // namespace details

#define NUMPY_NARY_FUNC_NAME i0
#define NUMPY_NARY_FUNC_SYM details::i0
#include "pythonic/types/numpy_nary_expr.hpp"
  } // namespace special
} // namespace scipy
PYTHONIC_NS_END

#endif
