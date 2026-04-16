#ifndef PYTHONIC_INCLUDE_DISPATCH_POP_HPP
#define PYTHONIC_INCLUDE_DISPATCH_POP_HPP

#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace __dispatch__
{
  template <class Any, class... Arg0>
  auto pop(Any &&any, Arg0 &&...arg0) -> decltype(any.pop(std::forward<Arg0>(arg0)...));

  DEFINE_FUNCTOR(pythonic::__dispatch__, pop);
} // namespace __dispatch__
PYTHONIC_NS_END

#endif
