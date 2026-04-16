#ifndef PYTHONIC_INCLUDE_OPERATOR_ISNOT_HPP
#define PYTHONIC_INCLUDE_OPERATOR_ISNOT_HPP

#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{

  template <class A, class B>
  auto is_not(A &&a, B &&b)
      -> decltype(builtins::id(std::forward<A>(a)) != builtins::id(std::forward<B>(b)));

  DEFINE_FUNCTOR(pythonic::operator_, is_not);
} // namespace operator_
PYTHONIC_NS_END

#endif
