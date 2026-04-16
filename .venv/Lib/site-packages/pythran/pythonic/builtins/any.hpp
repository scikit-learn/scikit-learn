#ifndef PYTHONIC_BUILTIN_ANY_HPP
#define PYTHONIC_BUILTIN_ANY_HPP

#include "pythonic/include/builtins/any.hpp"

#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{
  template <class Iterable>
  bool any(Iterable &&s)
  {
    auto iend = s.end();
    for (auto iter = s.begin(); iter != iend; ++iter)
      if (*iter)
        return true;
    return false;
  }
} // namespace builtins
PYTHONIC_NS_END

#endif
