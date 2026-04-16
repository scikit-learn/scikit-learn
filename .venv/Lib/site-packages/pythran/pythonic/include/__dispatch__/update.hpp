#ifndef PYTHONIC_INCLUDE_DISPATCH_UPDATE_HPP
#define PYTHONIC_INCLUDE_DISPATCH_UPDATE_HPP

#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace __dispatch__
{

  template <class Any, class... Arg0>
  auto update(Any &&any, Arg0 &&...arg0) -> decltype(any.update(std::forward<Arg0>(arg0)...));

  DEFINE_FUNCTOR(pythonic::__dispatch__, update);
} // namespace __dispatch__
PYTHONIC_NS_END

#endif
