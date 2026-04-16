#ifndef PYTHONIC_INCLUDE_BUILTIN_ORD_HPP
#define PYTHONIC_INCLUDE_BUILTIN_ORD_HPP

#include "pythonic/include/types/str.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{
  long ord(types::str const &v);

  DEFINE_FUNCTOR(pythonic::builtins, ord);
} // namespace builtins
PYTHONIC_NS_END

#endif
