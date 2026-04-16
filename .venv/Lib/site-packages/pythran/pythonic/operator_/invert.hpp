#ifndef PYTHONIC_OPERATOR_INVERT_HPP
#define PYTHONIC_OPERATOR_INVERT_HPP

#include "pythonic/include/operator_/invert.hpp"

#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{

  template <class A>
  auto invert(A &&a) -> decltype(~std::forward<A>(a))
  {
    return ~std::forward<A>(a);
  }
} // namespace operator_
PYTHONIC_NS_END

#endif
