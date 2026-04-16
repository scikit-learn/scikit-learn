#ifndef PYTHONIC_OPERATOR_NE_HPP
#define PYTHONIC_OPERATOR_NE_HPP

#include "pythonic/include/operator_/ne.hpp"

#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{

  template <class A, class B>
  auto ne(A &&a, B &&b) -> decltype(std::forward<A>(a) != std::forward<B>(b))
  {
    return std::forward<A>(a) != std::forward<B>(b);
  }

  inline bool ne(char const *a, char const *b)
  {
    return strcmp(a, b) != 0;
  }
} // namespace operator_
PYTHONIC_NS_END

#endif
