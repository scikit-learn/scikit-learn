#ifndef PYTHONIC_NUMPY_HEAVISIDE_HPP
#define PYTHONIC_NUMPY_HEAVISIDE_HPP

#include "pythonic/include/numpy/cos.hpp"

#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace details
  {

    template <class T0, class T1>
    T1 heaviside(T0 x0, T1 x1)
    {
      if (x0 == 0)
        return x1;
      if (x0 < 0)
        return 0;
      if (x0 > 0)
        return 1;
      return x0; // NaN
    }
  } // namespace details

#define NUMPY_NARY_FUNC_NAME heaviside
#define NUMPY_NARY_FUNC_SYM details::heaviside
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
