#ifndef PYTHONIC_INCLUDE_OPERATOR_MOD_HPP
#define PYTHONIC_INCLUDE_OPERATOR_MOD_HPP

#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{

  template <class A, class B>
  auto mod(A &&a, B &&b) -> std::enable_if_t<std::is_fundamental<std::decay_t<A>>::value &&
                                                 std::is_fundamental<std::decay_t<B>>::value,
                                             decltype(std::forward<A>(a) % std::forward<B>(b))>;

  inline double mod(double a, long b);

  inline double mod(double a, double b);

  template <class A, class B>
  auto mod(A &&a, B &&b) // for ndarrays
      -> std::enable_if_t<!std::is_fundamental<std::decay_t<A>>::value ||
                              !std::is_fundamental<std::decay_t<B>>::value,
                          decltype(std::forward<A>(a) % std::forward<B>(b))>;

  DEFINE_FUNCTOR(pythonic::operator_, mod);
} // namespace operator_
PYTHONIC_NS_END

#endif
