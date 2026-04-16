#ifndef PYTHONIC_INCLUDE_SCIPY_SPECIAL_NDTR_HPP
#define PYTHONIC_INCLUDE_SCIPY_SPECIAL_NDTR_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"

#include <xsimd/xsimd.hpp>

PYTHONIC_NS_BEGIN

namespace scipy
{
  namespace special
  {
    namespace details
    {
      template <class T>
      double ndtr(T x);
    }

#define NUMPY_NARY_FUNC_NAME ndtr
#define NUMPY_NARY_FUNC_SYM details::ndtr
#include "pythonic/include/types/numpy_nary_expr.hpp"
  } // namespace special
} // namespace scipy
PYTHONIC_NS_END

#endif
