#ifndef PYTHONIC_INCLUDE_OPERATOR_GE_HPP
#define PYTHONIC_INCLUDE_OPERATOR_GE_HPP

#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{
  template <class A, class B>
  auto ge(A &&a, B &&b) -> decltype(std::forward<A>(a) >= std::forward<B>(b));

  bool ge(char const *, char const *);

  DEFINE_FUNCTOR(pythonic::operator_, ge);
} // namespace operator_
PYTHONIC_NS_END

#endif
