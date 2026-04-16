#ifndef PYTHONIC_INCLUDE_SCIPY_SPECIAL_JV_HPP
#define PYTHONIC_INCLUDE_SCIPY_SPECIAL_JV_HPP

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
      double jv(T0 x, T1 y);
    }

#define NUMPY_NARY_FUNC_NAME jv
#define NUMPY_NARY_FUNC_SYM details::jv
#include "pythonic/include/types/numpy_nary_expr.hpp"
  } // namespace special
} // namespace scipy
PYTHONIC_NS_END

#endif
