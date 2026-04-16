#ifndef PYTHONIC_INCLUDE_OPERATOR_CONTAINS_HPP
#define PYTHONIC_INCLUDE_OPERATOR_CONTAINS_HPP

#include "pythonic/include/builtins/in.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{
  template <class A, class B>
  auto contains(A &&a, B &&b) -> decltype(in(std::forward<A>(a), std::forward<B>(b)));

  DEFINE_FUNCTOR(pythonic::operator_, contains);
} // namespace operator_
PYTHONIC_NS_END

#endif
