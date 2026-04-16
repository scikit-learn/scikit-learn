#ifndef PYTHONIC_SCIPY_SPECIAL_GAMMALN_HPP
#define PYTHONIC_SCIPY_SPECIAL_GAMMALN_HPP

#include "pythonic/include/scipy/special/gammaln.hpp"

#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace scipy
{
  namespace special
  {

#define NUMPY_NARY_FUNC_NAME gammaln
#define NUMPY_NARY_FUNC_SYM xsimd::lgamma
#include "pythonic/types/numpy_nary_expr.hpp"
  } // namespace special
} // namespace scipy
PYTHONIC_NS_END

#endif
