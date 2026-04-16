#ifndef PYTHONIC_INCLUDE_OPERATOR_EQ_HPP
#define PYTHONIC_INCLUDE_OPERATOR_EQ_HPP

#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{
  template <class A, class B>
  auto eq(A &&a, B &&b) -> decltype(std::forward<A>(a) == std::forward<B>(b));

  bool eq(char const *a, char const *b);

  DEFINE_FUNCTOR(pythonic::operator_, eq);
} // namespace operator_
PYTHONIC_NS_END

#endif
