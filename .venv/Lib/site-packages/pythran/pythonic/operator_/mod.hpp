#ifndef PYTHONIC_OPERATOR_MOD_HPP
#define PYTHONIC_OPERATOR_MOD_HPP

#include "pythonic/include/operator_/mod.hpp"

#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{

  template <class A, class B>
  auto mod(A &&a, B &&b) -> std::enable_if_t<std::is_fundamental<std::decay_t<A>>::value &&
                                                 std::is_fundamental<std::decay_t<B>>::value,
                                             decltype(std::forward<A>(a) % std::forward<B>(b))>
  {
    auto t = std::forward<A>(a) % b;
    return t < 0 ? (t + b) : t;
  }

  inline double mod(double a, long b)
  {
    auto t = std::fmod(a, double(b));
    return t < 0 ? (t + b) : t;
  }

  inline double mod(double a, double b)
  {
    auto t = std::fmod(a, b);
    return t < 0 ? (t + b) : t;
  }

  template <class A, class B>
  auto mod(A &&a, B &&b) // for ndarrays
      -> std::enable_if_t<!std::is_fundamental<std::decay_t<A>>::value ||
                              !std::is_fundamental<std::decay_t<B>>::value,
                          decltype(std::forward<A>(a) % std::forward<B>(b))>
  {
    return std::forward<A>(a) % std::forward<B>(b);
  }
} // namespace operator_
PYTHONIC_NS_END

#endif
