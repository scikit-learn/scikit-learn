#ifndef PYTHONIC_INCLUDE_BUILTIN_NEXT_HPP
#define PYTHONIC_INCLUDE_BUILTIN_NEXT_HPP

#include "pythonic/include/utils/functor.hpp"

#include <utility>

PYTHONIC_NS_BEGIN

namespace builtins
{

  template <class T>
  auto next(T &&y) -> decltype(*y);

  DEFINE_FUNCTOR(pythonic::builtins, next);
} // namespace builtins
PYTHONIC_NS_END

#endif
