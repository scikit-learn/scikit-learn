#ifndef PYTHONIC_DISPATCH_COPY_HPP
#define PYTHONIC_DISPATCH_COPY_HPP

#include "pythonic/include/__dispatch__/copy.hpp"

#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace __dispatch__
{
  template <class Any>
  auto copy(Any const &any) -> decltype(any.copy())
  {
    return any.copy();
  }
} // namespace __dispatch__
PYTHONIC_NS_END

#endif
