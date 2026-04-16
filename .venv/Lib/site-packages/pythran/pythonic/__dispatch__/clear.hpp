#ifndef PYTHONIC_DISPATCH_CLEAR_HPP
#define PYTHONIC_DISPATCH_CLEAR_HPP

#include "pythonic/include/__dispatch__/clear.hpp"

#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace __dispatch__
{

  template <class Any>
  auto clear(Any &&any) -> decltype(any.clear());
}
PYTHONIC_NS_END

#endif
