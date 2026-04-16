#ifndef PYTHONIC_INCLUDE_OPERATOR_OR_HPP
#define PYTHONIC_INCLUDE_OPERATOR_OR_HPP

#include "pythonic/include/operator_/overloads.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{

  template <class A, class B>
  auto or_(A &&a, B &&b) -> decltype(std::forward<A>(a) | std::forward<B>(b));

  DEFINE_ALL_OPERATOR_OVERLOADS_DECL(or_, |)

  DEFINE_FUNCTOR(pythonic::operator_, or_);
} // namespace operator_
PYTHONIC_NS_END

#endif
