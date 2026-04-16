#ifndef PYTHONIC_INCLUDE_DISPATCH_CONJUGATE_HPP
#define PYTHONIC_INCLUDE_DISPATCH_CONJUGATE_HPP

#include "pythonic/include/numpy/conjugate.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace __dispatch__
{
  template <class Any>
  auto conjugate(Any const &any) -> decltype(numpy::functor::conjugate{}(any));

  DEFINE_FUNCTOR(pythonic::__dispatch__, conjugate);
} // namespace __dispatch__
PYTHONIC_NS_END

#endif
