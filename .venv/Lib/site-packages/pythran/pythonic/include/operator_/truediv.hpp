#ifndef PYTHONIC_INCLUDE_OPERATOR_TRUEDIV_HPP
#define PYTHONIC_INCLUDE_OPERATOR_TRUEDIV_HPP

#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{
  template <class A, class B>
  auto truediv(A &&a, B &&b) -> decltype(std::forward<A>(a) / (double)std::forward<B>(b));

  DEFINE_FUNCTOR(pythonic::operator_, truediv);
} // namespace operator_
PYTHONIC_NS_END

#endif
