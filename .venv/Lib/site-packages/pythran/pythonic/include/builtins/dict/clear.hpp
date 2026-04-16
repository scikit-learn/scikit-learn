#ifndef PYTHONIC_INCLUDE_BUILTIN_DICT_CLEAR_HPP
#define PYTHONIC_INCLUDE_BUILTIN_DICT_CLEAR_HPP

#include "pythonic/include/__dispatch__/clear.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN
namespace builtins
{
  namespace dict
  {
    USING_FUNCTOR(clear, pythonic::__dispatch__::functor::clear);
  }
} // namespace builtins
PYTHONIC_NS_END

#endif
