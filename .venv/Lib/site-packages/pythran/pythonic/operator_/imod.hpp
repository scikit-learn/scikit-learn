#ifndef PYTHONIC_OPERATOR_IMOD_HPP
#define PYTHONIC_OPERATOR_IMOD_HPP

#include "pythonic/include/operator_/imod.hpp"

#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{
  template <class A, class B>
  A imod(A const &a, B &&b)
  {
    return a % std::forward<B>(b);
  }

  template <class A, class B>
  A &imod(A &a, B &&b)
  {
    return a %= std::forward<B>(b);
  }
} // namespace operator_
PYTHONIC_NS_END

#endif
