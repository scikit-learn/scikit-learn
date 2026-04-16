#ifndef PYTHONIC_INCLUDE_OPERATOR_RSHIFT_HPP
#define PYTHONIC_INCLUDE_OPERATOR_RSHIFT_HPP

#include "pythonic/include/operator_/overloads.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{
  template <class A, class B>
  auto rshift(A &&a, B &&b) -> decltype(std::forward<A>(a) >> std::forward<B>(b));

  DEFINE_ALL_OPERATOR_OVERLOADS_DECL(rshift, >>)

  DEFINE_FUNCTOR(pythonic::operator_, rshift);
} // namespace operator_
PYTHONIC_NS_END

#endif
