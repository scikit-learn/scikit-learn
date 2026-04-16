#ifndef PYTHONIC_INCLUDE_NUMPY_MOD_HPP
#define PYTHONIC_INCLUDE_NUMPY_MOD_HPP

#include "pythonic/include/operator_/mod.hpp"
#include "pythonic/include/types/assignable.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  /* this is still a numpy_expr, because operator::mod_ forwards to
   * operator% which is correctly overloaded
   */

  USING_FUNCTOR(mod, operator_::functor::mod);
} // namespace numpy
PYTHONIC_NS_END

#endif
