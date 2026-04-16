#ifndef PYTHONIC_INCLUDE_SCIPY_SPECIAL_BINOM_HPP
#define PYTHONIC_INCLUDE_SCIPY_SPECIAL_BINOM_HPP

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
      double binom(T0 n, T1 k);
    }

#define NUMPY_NARY_FUNC_NAME binom
#define NUMPY_NARY_FUNC_SYM details::binom
#include "pythonic/include/types/numpy_nary_expr.hpp"
  } // namespace special
} // namespace scipy
PYTHONIC_NS_END

#endif
