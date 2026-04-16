#ifndef PYTHONIC_DISPATCH_EXTEND_HPP
#define PYTHONIC_DISPATCH_EXTEND_HPP

#include "pythonic/include/__dispatch__/extend.hpp"

#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace __dispatch__
{

  template <class Any, class Arg>
  types::none_type extend(Any &&any, Arg &&arg)
  {
    std::forward<Any>(any) += std::forward<Arg>(arg);
    return {};
  }
} // namespace __dispatch__
PYTHONIC_NS_END

#endif
