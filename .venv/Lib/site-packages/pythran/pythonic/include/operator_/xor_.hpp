#ifndef PYTHONIC_INCLUDE_OPERATOR_XOR_HPP
#define PYTHONIC_INCLUDE_OPERATOR_XOR_HPP

#include "pythonic/include/operator_/overloads.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{
  template <class A, class B>
  auto xor_(A &&a, B &&b) -> decltype(std::forward<A>(a) ^ std::forward<B>(b));

  DEFINE_ALL_OPERATOR_OVERLOADS_DECL(xor_, ^)

  DEFINE_FUNCTOR(pythonic::operator_, xor_);
} // namespace operator_
PYTHONIC_NS_END

#endif
