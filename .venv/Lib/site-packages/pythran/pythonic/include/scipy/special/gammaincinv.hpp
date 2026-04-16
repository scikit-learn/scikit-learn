#ifndef PYTHONIC_INCLUDE_SCIPY_SPECIAL_GAMMAINCINV_HPP
#define PYTHONIC_INCLUDE_SCIPY_SPECIAL_GAMMAINCINV_HPP

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
      template <class T0, class T1>
      double gammaincinv(T0 x, T1 y);
    }

#define NUMPY_NARY_FUNC_NAME gammaincinv
#define NUMPY_NARY_FUNC_SYM details::gammaincinv
#include "pythonic/include/types/numpy_nary_expr.hpp"
  } // namespace special
} // namespace scipy
PYTHONIC_NS_END

#endif
