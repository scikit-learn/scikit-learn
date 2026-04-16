#ifndef PYTHONIC_BUILTIN_NEXT_HPP
#define PYTHONIC_BUILTIN_NEXT_HPP

#include "pythonic/include/builtins/next.hpp"

#include "pythonic/builtins/StopIteration.hpp"
#include "pythonic/utils/functor.hpp"

#include <utility>

PYTHONIC_NS_BEGIN

namespace builtins
{

  template <class T>
  auto next(T &&y) -> decltype(*y)
  {
    if ((decltype(y.begin()) &)y != y.end()) {
      auto tmp = *y;
      ++y;
      return tmp;
    } else
      throw types::StopIteration();
  }
} // namespace builtins
PYTHONIC_NS_END

#endif
