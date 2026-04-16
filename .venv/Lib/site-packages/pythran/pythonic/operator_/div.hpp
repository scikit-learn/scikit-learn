#ifndef PYTHONIC_OPERATOR_DIV_HPP
#define PYTHONIC_OPERATOR_DIV_HPP

#include "pythonic/include/operator_/div.hpp"

#include "pythonic/operator_/overloads.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{

  template <class A, class B>
  auto div(A &&a, B &&b) // for ndarrays
      -> std::enable_if_t<!std::is_fundamental<std::decay_t<A>>::value ||
                              !std::is_fundamental<std::decay_t<B>>::value,
                          decltype(std::forward<A>(a) / std::forward<B>(b))>
  {
    return std::forward<A>(a) / std::forward<B>(b);
  }

  inline double div(double a, double b)
  {
    assert(b != 0 && "divide by zero");
    return a / b;
  }
} // namespace operator_
PYTHONIC_NS_END

#endif
