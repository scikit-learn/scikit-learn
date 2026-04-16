#ifndef PYTHONIC_INCLUDE_DISPATCH_APPEND_HPP
#define PYTHONIC_INCLUDE_DISPATCH_APPEND_HPP

#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace __dispatch__
{

  template <class Any, class Arg>
  types::none_type append(Any &&any, Arg &&arg0);

  DEFINE_FUNCTOR(pythonic::__dispatch__, append);
} // namespace __dispatch__
PYTHONIC_NS_END

#endif
