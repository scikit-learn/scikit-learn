#ifndef PYTHONIC_DISPATCH_APPEND_HPP
#define PYTHONIC_DISPATCH_APPEND_HPP

#include "pythonic/include/__dispatch__/append.hpp"

#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace __dispatch__
{

  template <class Any, class Arg>
  types::none_type append(Any &&any, Arg &&arg)
  {
    std::forward<Any>(any).push_back(std::forward<Arg>(arg));
    return {};
  }
} // namespace __dispatch__
PYTHONIC_NS_END

#endif
