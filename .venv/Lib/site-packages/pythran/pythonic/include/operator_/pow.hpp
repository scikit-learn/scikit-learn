#ifndef PYTHONIC_INCLUDE_OPERATOR_POW_HPP
#define PYTHONIC_INCLUDE_OPERATOR_POW_HPP

#include "pythonic/include/builtins/pow.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{
  USING_FUNCTOR(pow, pythonic::builtins::functor::pow);
}

PYTHONIC_NS_END

#endif
