#ifndef PYTHONIC_INCLUDE_BUILTIN_SET_UPDATE_HPP
#define PYTHONIC_INCLUDE_BUILTIN_SET_UPDATE_HPP

#include "pythonic/include/__dispatch__/update.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN
namespace builtins
{
  namespace set
  {
    USING_FUNCTOR(update, pythonic::__dispatch__::functor::update);
  }
} // namespace builtins
PYTHONIC_NS_END
#endif
