#ifndef PYTHONIC_OPERATOR_SUB_HPP
#define PYTHONIC_OPERATOR_SUB_HPP

#include "pythonic/include/operator_/sub.hpp"

#include "pythonic/operator_/overloads.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{

  template <class A, class B>
  auto sub(A &&a, B &&b) -> decltype(std::forward<A>(a) - std::forward<B>(b))
  {
    return std::forward<A>(a) - std::forward<B>(b);
  }

  DEFINE_ALL_OPERATOR_OVERLOADS_IMPL(sub, -,
                                     (((b < 0)
                                           ? (a <= std::numeric_limits<decltype(b)>::max() + b)
                                           : (std::numeric_limits<decltype(b)>::min() + b <= a))))
} // namespace operator_
PYTHONIC_NS_END

#endif
