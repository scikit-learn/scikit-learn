#ifndef PYTHONIC_INCLUDE_DISPATCH_EXTEND_HPP
#define PYTHONIC_INCLUDE_DISPATCH_EXTEND_HPP

#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace __dispatch__
{

  template <class Any, class Arg>
  types::none_type extend(Any &&any, Arg &&arg0);

  DEFINE_FUNCTOR(pythonic::__dispatch__, extend);
} // namespace __dispatch__
PYTHONIC_NS_END

#endif
