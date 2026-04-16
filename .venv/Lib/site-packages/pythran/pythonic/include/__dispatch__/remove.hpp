#ifndef PYTHONIC_INCLUDE_DISPATCH_REMOVE_HPP
#define PYTHONIC_INCLUDE_DISPATCH_REMOVE_HPP

#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace __dispatch__
{
  template <class Any, class Arg0>
  auto remove(Any &any, Arg0 const &arg0) -> decltype(any.remove(arg0));

  DEFINE_FUNCTOR(pythonic::__dispatch__, remove);
} // namespace __dispatch__
PYTHONIC_NS_END

#endif
