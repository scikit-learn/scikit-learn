#ifndef PYTHONIC_OPERATOR_LE_HPP
#define PYTHONIC_OPERATOR_LE_HPP

#include "pythonic/include/operator_/le.hpp"

#include "pythonic/utils/functor.hpp"
#include <cstring>

PYTHONIC_NS_BEGIN

namespace operator_
{
  template <class A, class B>
  auto le(A &&a, B &&b) -> decltype(std::forward<A>(a) <= std::forward<B>(b))
  {
    return std::forward<A>(a) <= std::forward<B>(b);
  }
  inline bool le(char const *self, char const *other)
  {
    return strcmp(self, other) <= 0;
  }
} // namespace operator_
PYTHONIC_NS_END

#endif
