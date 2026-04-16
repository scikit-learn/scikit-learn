#ifndef PYTHONIC_OPERATOR_ITRUEDIV_HPP
#define PYTHONIC_OPERATOR_ITRUEDIV_HPP

#include "pythonic/include/operator_/itruediv.hpp"

#include "pythonic/operator_/truediv.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{
  template <class A, class B>
  auto itruediv(A const &a, B &&b) -> decltype(truediv(a, std::forward<B>(b)))
  {
    return truediv(a, std::forward<B>(b));
  }
  template <class A, class B>
  auto itruediv(A &a, B &&b)
      -> std::enable_if_t<std::is_same<A, decltype(truediv(a, std::forward<B>(b)))>::value, A &>
  {
    return a = truediv(a, std::forward<B>(b));
  }
  template <class A, class B>
  auto itruediv(A &a, B &&b)
      -> std::enable_if_t<!std::is_same<A, decltype(truediv(a, std::forward<B>(b)))>::value,
                          decltype(truediv(a, std::forward<B>(b)))>
  {
    return truediv(a, std::forward<B>(b));
  }
} // namespace operator_
PYTHONIC_NS_END

#endif
