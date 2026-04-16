#ifndef PYTHONIC_INCLUDE_OPERATOR_GT_HPP
#define PYTHONIC_INCLUDE_OPERATOR_GT_HPP

#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{
  template <class A, class B>
  auto gt(A &&a, B &&b) -> decltype(std::forward<A>(a) > std::forward<B>(b));

  bool gt(char const *, char const *);

  DEFINE_FUNCTOR(pythonic::operator_, gt);
} // namespace operator_
PYTHONIC_NS_END

#endif
