#ifndef PYTHONIC_INCLUDE_OPERATOR_NEG_HPP
#define PYTHONIC_INCLUDE_OPERATOR_NEG_HPP

#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{

  template <class A>
  auto neg(A &&a) -> decltype(-std::forward<A>(a));

  DEFINE_FUNCTOR(pythonic::operator_, neg);
} // namespace operator_
PYTHONIC_NS_END

#endif
