#ifndef PYTHONIC_OPERATOR_NEG_HPP
#define PYTHONIC_OPERATOR_NEG_HPP

#include "pythonic/include/operator_/neg.hpp"

#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{

  template <class A>
  auto neg(A &&a) -> decltype(-std::forward<A>(a))
  {
    return -std::forward<A>(a);
  }
} // namespace operator_
PYTHONIC_NS_END

#endif
