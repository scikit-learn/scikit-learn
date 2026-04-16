#ifndef PYTHONIC_INCLUDE_DISPATCH_COPY_HPP
#define PYTHONIC_INCLUDE_DISPATCH_COPY_HPP

#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace __dispatch__
{
  template <class Any>
  auto copy(Any const &any) -> decltype(any.copy());

  DEFINE_FUNCTOR(pythonic::__dispatch__, copy);
} // namespace __dispatch__
PYTHONIC_NS_END

#endif
