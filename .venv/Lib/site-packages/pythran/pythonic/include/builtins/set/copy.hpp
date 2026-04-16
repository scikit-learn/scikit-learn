#ifndef PYTHONIC_INCLUDE_BUILTIN_SET_COPY_HPP
#define PYTHONIC_INCLUDE_BUILTIN_SET_COPY_HPP

#include "pythonic/include/__dispatch__/copy.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN
namespace builtins
{
  namespace set
  {
    USING_FUNCTOR(copy, pythonic::__dispatch__::functor::copy);
  }
} // namespace builtins
PYTHONIC_NS_END
#endif
