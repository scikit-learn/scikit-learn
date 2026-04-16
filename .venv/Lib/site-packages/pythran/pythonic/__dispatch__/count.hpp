#ifndef PYTHONIC_DISPATCH_COUNT_HPP
#define PYTHONIC_DISPATCH_COUNT_HPP

#include "pythonic/include/__dispatch__/count.hpp"

#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace __dispatch__
{

  template <class Any, class Value>
  auto count(Any &&any, Value &&value) -> decltype(any.count(std::forward<Value>(value)))
  {
    return any.count(std::forward<Value>(value));
  }
} // namespace __dispatch__
PYTHONIC_NS_END

#endif
