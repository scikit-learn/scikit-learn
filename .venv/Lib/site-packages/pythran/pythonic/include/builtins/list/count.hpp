#ifndef PYTHONIC_INCLUDE_BUILTIN_LIST_COUNT_HPP
#define PYTHONIC_INCLUDE_BUILTIN_LIST_COUNT_HPP

#include "pythonic/include/__dispatch__/count.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace list
  {
    USING_FUNCTOR(count, pythonic::__dispatch__::functor::count);
  }
} // namespace builtins
PYTHONIC_NS_END

#endif
