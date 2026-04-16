#ifndef PYTHONIC_INCLUDE_OPERATOR_IMIN_HPP
#define PYTHONIC_INCLUDE_OPERATOR_IMIN_HPP

#include "pythonic/include/numpy/minimum.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{
  template <class A, class B>
  auto imin(A &&a, B &&b)
      -> std::enable_if_t<std::is_const<A>::value || !std::is_assignable<A, B>::value,
                          decltype(numpy::functor::minimum{}(std::forward<A>(a),
                                                             std::forward<B>(b)))>;

  template <class A, class B>
  auto imin(A &&a, B &&b)
      -> std::enable_if_t<!std::is_const<A>::value && std::is_assignable<A, B>::value,
                          decltype(a = numpy::functor::minimum{}(std::forward<A>(a),
                                                                 std::forward<B>(b)))>;

  DEFINE_FUNCTOR(pythonic::operator_, imin);
} // namespace operator_
PYTHONIC_NS_END

#endif
