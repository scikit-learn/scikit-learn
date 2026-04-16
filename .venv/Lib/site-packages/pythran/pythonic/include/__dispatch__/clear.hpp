#ifndef PYTHONIC_INCLUDE_DISPATCH_CLEAR_HPP
#define PYTHONIC_INCLUDE_DISPATCH_CLEAR_HPP

#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace __dispatch__
{
  template <class Any>
  auto clear(Any &&any) -> decltype(any.clear())
  {
    return any.clear();
  }

  DEFINE_FUNCTOR(pythonic::__dispatch__, clear);
} // namespace __dispatch__
PYTHONIC_NS_END

#endif
