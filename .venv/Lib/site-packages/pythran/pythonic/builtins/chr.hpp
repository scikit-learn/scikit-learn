#ifndef PYTHONIC_BUILTIN_CHR_HPP
#define PYTHONIC_BUILTIN_CHR_HPP

#include "pythonic/include/builtins/chr.hpp"

#include "pythonic/types/str.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{
  template <class T>
  types::str chr(T const &v)
  {
    return types::str((char)v);
  }
} // namespace builtins
PYTHONIC_NS_END

#endif
