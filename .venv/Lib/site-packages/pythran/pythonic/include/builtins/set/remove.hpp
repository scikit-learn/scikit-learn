#ifndef PYTHONIC_INCLUDE_BUILTIN_SET_REMOVE_HPP
#define PYTHONIC_INCLUDE_BUILTIN_SET_REMOVE_HPP

#include "pythonic/include/__dispatch__/remove.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN
namespace builtins
{
  namespace set
  {
    USING_FUNCTOR(remove, pythonic::__dispatch__::functor::remove);
  }
} // namespace builtins
PYTHONIC_NS_END
#endif
