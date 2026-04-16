#ifndef PYTHONIC_INCLUDE_DISPATCH_REVERSE_HPP
#define PYTHONIC_INCLUDE_DISPATCH_REVERSE_HPP

#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace __dispatch__
{

  template <class Any>
  types::none_type reverse(Any &&any);

  DEFINE_FUNCTOR(pythonic::__dispatch__, reverse);
} // namespace __dispatch__
PYTHONIC_NS_END

#endif
