#ifndef PYTHONIC_DISPATCH_REMOVE_HPP
#define PYTHONIC_DISPATCH_REMOVE_HPP

#include "pythonic/include/__dispatch__/remove.hpp"

#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace __dispatch__
{
  template <class Any, class Arg0>
  auto remove(Any &any, Arg0 const &arg0) -> decltype(any.remove(arg0))
  {
    return any.remove(arg0);
  }
} // namespace __dispatch__
PYTHONIC_NS_END

#endif
