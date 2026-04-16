#ifndef PYTHONIC_INCLUDE_SCIPY_SPECIAL_SPHERICAL_YN_HPP
#define PYTHONIC_INCLUDE_SCIPY_SPECIAL_SPHERICAL_YN_HPP

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
      double spherical_yn(T0 v, T1 x, bool derivative = false);
    }

#define NUMPY_NARY_FUNC_NAME spherical_yn
#define NUMPY_NARY_FUNC_SYM details::spherical_yn
#include "pythonic/include/types/numpy_nary_expr.hpp"
  } // namespace special
} // namespace scipy
PYTHONIC_NS_END

#endif
