#ifndef PYTHONIC_OPERATOR_CONCAT_HPP
#define PYTHONIC_OPERATOR_CONCAT_HPP

#include "pythonic/include/operator_/concat.hpp"

#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{

  template <class A, class B>
  auto concat(A &&a, B &&b) -> decltype(std::forward<A>(a) + std::forward<B>(b))
  {
    return std::forward<A>(a) + std::forward<B>(b);
  }
} // namespace operator_
PYTHONIC_NS_END

#endif
