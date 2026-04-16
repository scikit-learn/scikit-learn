#ifndef PYTHONIC_OPERATOR_RSHIFT_HPP
#define PYTHONIC_OPERATOR_RSHIFT_HPP

#include "pythonic/include/operator_/rshift.hpp"
#include "pythonic/operator_/overloads.hpp"

#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{
  template <class A, class B>
  auto rshift(A &&a, B &&b) -> decltype(std::forward<A>(a) >> std::forward<B>(b))
  {
    return std::forward<A>(a) >> std::forward<B>(b);
  }

  DEFINE_ALL_OPERATOR_OVERLOADS_IMPL(rshift, >>, true)
} // namespace operator_
PYTHONIC_NS_END

#endif
