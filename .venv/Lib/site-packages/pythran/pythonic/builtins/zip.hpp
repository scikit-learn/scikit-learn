#ifndef PYTHONIC_BUILTIN_ZIP_HPP
#define PYTHONIC_BUILTIN_ZIP_HPP

#include "pythonic/include/builtins/zip.hpp"

#include "pythonic/builtins/None.hpp"
#include "pythonic/builtins/map.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  template <typename... Iter>
  auto zip(Iter &&...iters) -> decltype(map(builtins::None, std::forward<Iter>(iters)...))
  {
    return map(builtins::None, std::forward<Iter>(iters)...);
  }
} // namespace builtins
PYTHONIC_NS_END

#endif
