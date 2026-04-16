#ifndef PYTHONIC_INCLUDE_OPERATOR_IFLOORDIV_HPP
#define PYTHONIC_INCLUDE_OPERATOR_IFLOORDIV_HPP

#include "pythonic/include/operator_/mod.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{
  template <class A, class B>
  A ifloordiv(A &&a, B &&b);

  template <class A, class B>
  auto ifloordiv(A const &a, B &&b) -> decltype((a - mod(a, b)) / b);

  DEFINE_FUNCTOR(pythonic::operator_, ifloordiv);
} // namespace operator_
PYTHONIC_NS_END

#endif
