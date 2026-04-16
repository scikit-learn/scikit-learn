#ifndef PYTHONIC_INCLUDE_SCIPY_SPECIAL_HANKEL2_HPP
#define PYTHONIC_INCLUDE_SCIPY_SPECIAL_HANKEL2_HPP

#include "pythonic/include/types/complex.hpp"
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
      std::complex<double> hankel2(T0 x, T1 y);
    }

#define NUMPY_NARY_FUNC_NAME hankel2
#define NUMPY_NARY_FUNC_SYM details::hankel2
#include "pythonic/include/types/numpy_nary_expr.hpp"
  } // namespace special
} // namespace scipy
PYTHONIC_NS_END

#endif
