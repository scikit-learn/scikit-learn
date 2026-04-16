#ifndef PYTHONIC_BUILTIN_ALL_HPP
#define PYTHONIC_BUILTIN_ALL_HPP

#include "pythonic/include/builtins/all.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  template <class Iterable>
  bool all(Iterable &&s)
  {
    auto iend = s.end();
    for (auto iter = s.begin(); iter != iend; ++iter)
      if (!*iter)
        return false;
    return true;
  }
} // namespace builtins
PYTHONIC_NS_END

#endif
