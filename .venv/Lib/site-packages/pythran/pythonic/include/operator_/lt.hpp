#ifndef PYTHONIC_INCLUDE_OPERATOR_LT_HPP
#define PYTHONIC_INCLUDE_OPERATOR_LT_HPP

#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{

  template <class A, class B>
  auto lt(A &&a, B &&b) -> decltype(std::forward<A>(a) < std::forward<B>(b));
  bool lt(char const *self, char const *other);

  DEFINE_FUNCTOR(pythonic::operator_, lt);
} // namespace operator_
PYTHONIC_NS_END

#endif
