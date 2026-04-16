#ifndef PYTHONIC_DISPATCH_UPDATE_HPP
#define PYTHONIC_DISPATCH_UPDATE_HPP

#include "pythonic/include/__dispatch__/update.hpp"

#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace __dispatch__
{

  template <class Any, class... Arg0>
  auto update(Any &&any, Arg0 &&...arg0) -> decltype(any.update(std::forward<Arg0>(arg0)...))
  {
    return any.update(std::forward<Arg0>(arg0)...);
  }
} // namespace __dispatch__
PYTHONIC_NS_END

#endif
