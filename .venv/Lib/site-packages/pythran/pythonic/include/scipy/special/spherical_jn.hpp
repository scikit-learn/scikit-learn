#ifndef PYTHONIC_INCLUDE_SCIPY_SPECIAL_SPHERICAL_JN_HPP
#define PYTHONIC_INCLUDE_SCIPY_SPECIAL_SPHERICAL_JN_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace scipy
{
  namespace special
  {

    namespace details
    {
      template <class T0, class T1>
      double spherical_jn(T0 v, T1 x, bool derivative = false);
    }

#define NUMPY_NARY_FUNC_NAME spherical_jn
#define NUMPY_NARY_FUNC_SYM details::spherical_jn
#include "pythonic/include/types/numpy_nary_expr.hpp"
  } // namespace special
} // namespace scipy
PYTHONIC_NS_END

#endif
