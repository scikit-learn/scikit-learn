#ifndef PYTHONIC_OPERATOR_IPOW_HPP
#define PYTHONIC_OPERATOR_IPOW_HPP

#include "pythonic/include/operator_/ipow.hpp"

#include "pythonic/builtins/pow.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{

  template <class A, class B>
  A ipow(A const &a, B &&b)
  {
    return builtins::pow(a, std::forward<B>(b));
  }
  template <class A, class B>
  A &ipow(A &a, B &&b)
  {
    return a = builtins::pow(a, std::forward<B>(b));
  }
} // namespace operator_
PYTHONIC_NS_END

#endif
