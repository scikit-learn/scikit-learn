#ifndef PYTHONIC_OPERATOR_OR_HPP
#define PYTHONIC_OPERATOR_OR_HPP

#include "pythonic/include/operator_/or_.hpp"
#include "pythonic/operator_/overloads.hpp"

#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{

  template <class A, class B>
  auto or_(A &&a, B &&b) -> decltype(std::forward<A>(a) | std::forward<B>(b))
  {
    return std::forward<A>(a) | std::forward<B>(b);
  }

  DEFINE_ALL_OPERATOR_OVERLOADS_IMPL(or_, |, true)
} // namespace operator_
PYTHONIC_NS_END

#endif
