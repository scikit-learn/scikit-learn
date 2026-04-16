#ifndef PYTHONIC_INCLUDE_SCIPY_SPECIAL_I0E_HPP
#define PYTHONIC_INCLUDE_SCIPY_SPECIAL_I0E_HPP
/*
 * Adapted from Cephes Math Library Release 2.8:  June, 2000
 * Copyright 1984, 1987, 2000 by Stephen L. Moshier
 */

#include "pythonic/include/scipy/special/i0.hpp"

PYTHONIC_NS_BEGIN

namespace scipy
{
  namespace special
  {

    namespace details
    {
      template <class T>
      double i0e(T x);
    }

#define NUMPY_NARY_FUNC_NAME i0e
#define NUMPY_NARY_FUNC_SYM details::i0e
#include "pythonic/include/types/numpy_nary_expr.hpp"
  } // namespace special
} // namespace scipy
PYTHONIC_NS_END

#endif
