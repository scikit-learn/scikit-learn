#ifndef PYTHONIC_INCLUDE_OPERATOR_INVERT_HPP
#define PYTHONIC_INCLUDE_OPERATOR_INVERT_HPP

#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{

  template <class A>
  auto invert(A &&a) -> decltype(~std::forward<A>(a));

  DEFINE_FUNCTOR(pythonic::operator_, invert);
} // namespace operator_
PYTHONIC_NS_END

#endif
