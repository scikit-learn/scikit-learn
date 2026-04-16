#ifndef PYTHONIC_INCLUDE_FUNCTOOLS_REDUCE_HPP
#define PYTHONIC_INCLUDE_FUNCTOOLS_REDUCE_HPP

#include "pythonic/include/builtins/reduce.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace functools
{
  USING_FUNCTOR(reduce, builtins::functor::reduce);
}
PYTHONIC_NS_END

#endif
