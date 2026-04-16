#ifndef PYTHONIC_INCLUDE_OPERATOR_ABS__HPP
#define PYTHONIC_INCLUDE_OPERATOR_ABS__HPP

#include "pythonic/include/builtins/abs.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{
  USING_FUNCTOR(__abs__, builtins::functor::abs);
}
PYTHONIC_NS_END

#endif
