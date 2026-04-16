#ifndef PYTHONIC_DISPATCH_REVERSE_HPP
#define PYTHONIC_DISPATCH_REVERSE_HPP

#include "pythonic/include/__dispatch__/reverse.hpp"

#include "pythonic/utils/functor.hpp"

#include <algorithm>

PYTHONIC_NS_BEGIN

namespace __dispatch__
{

  template <class Any>
  types::none_type reverse(Any &&any)
  {
    std::reverse(any.begin(), any.end());
    return {};
  }
} // namespace __dispatch__
PYTHONIC_NS_END

#endif
