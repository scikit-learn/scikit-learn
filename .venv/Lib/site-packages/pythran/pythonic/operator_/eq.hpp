#ifndef PYTHONIC_OPERATOR_EQ_HPP
#define PYTHONIC_OPERATOR_EQ_HPP

#include "pythonic/include/operator_/eq.hpp"

#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{
  template <class A, class B>
  auto eq(A &&a, B &&b) -> decltype(std::forward<A>(a) == std::forward<B>(b))
  {
    return std::forward<A>(a) == std::forward<B>(b);
  }

  inline bool eq(char const *a, char const *b)
  {
    return strcmp(a, b) == 0;
  }
} // namespace operator_
PYTHONIC_NS_END

#endif
