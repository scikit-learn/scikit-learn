#ifndef PYTHONIC_OPERATOR_TRUEDIV_HPP
#define PYTHONIC_OPERATOR_TRUEDIV_HPP

#include "pythonic/include/operator_/truediv.hpp"

#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{
  template <class A, class B>
  auto truediv(A &&a, B &&b) -> decltype(std::forward<A>(a) / (double)std::forward<B>(b))
  {
    return std::forward<A>(a) / ((double)std::forward<B>(b));
  }
} // namespace operator_
PYTHONIC_NS_END

#endif
