#ifndef PYTHONIC_INCLUDE_OPERATOR_SUB_HPP
#define PYTHONIC_INCLUDE_OPERATOR_SUB_HPP

#include "pythonic/include/operator_/overloads.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{

  template <class A, class B>
  auto sub(A &&a, B &&b) -> decltype(std::forward<A>(a) - std::forward<B>(b));

  DEFINE_ALL_OPERATOR_OVERLOADS_DECL(sub, -)

  DEFINE_FUNCTOR(pythonic::operator_, sub);
} // namespace operator_
PYTHONIC_NS_END

#endif
