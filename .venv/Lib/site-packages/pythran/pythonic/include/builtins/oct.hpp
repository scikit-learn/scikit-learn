#ifndef PYTHONIC_INCLUDE_BUILTIN_OCT_HPP
#define PYTHONIC_INCLUDE_BUILTIN_OCT_HPP

#include "pythonic/include/types/str.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{
  template <class T>
  types::str oct(T const &v);

  DEFINE_FUNCTOR(pythonic::builtins, oct);
} // namespace builtins
PYTHONIC_NS_END

#endif
