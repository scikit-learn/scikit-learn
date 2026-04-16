#ifndef PYTHONIC_INCLUDE_BUILTIN_ZIP_HPP
#define PYTHONIC_INCLUDE_BUILTIN_ZIP_HPP

#include "pythonic/include/builtins/None.hpp"
#include "pythonic/include/builtins/map.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  template <typename... Iter>
  auto zip(Iter &&...iters) -> decltype(map(builtins::None, std::forward<Iter>(iters)...));

  DEFINE_FUNCTOR(pythonic::builtins, zip);
} // namespace builtins
PYTHONIC_NS_END

#endif
