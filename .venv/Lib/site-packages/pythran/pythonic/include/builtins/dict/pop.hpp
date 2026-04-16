#ifndef PYTHONIC_INCLUDE_BUILTIN_DICT_POP_HPP
#define PYTHONIC_INCLUDE_BUILTIN_DICT_POP_HPP

#include "pythonic/include/__dispatch__/pop.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN
namespace builtins
{
  namespace dict
  {
    USING_FUNCTOR(pop, pythonic::__dispatch__::functor::pop);
  }
} // namespace builtins
PYTHONIC_NS_END

#endif
