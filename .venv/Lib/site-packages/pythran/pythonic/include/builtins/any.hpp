#ifndef PYTHONIC_INCLUDE_BUILTIN_ANY_HPP
#define PYTHONIC_INCLUDE_BUILTIN_ANY_HPP

#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{
  template <class Iterable>
  bool any(Iterable &&s);

  DEFINE_FUNCTOR(pythonic::builtins, any);
} // namespace builtins
PYTHONIC_NS_END

#endif
