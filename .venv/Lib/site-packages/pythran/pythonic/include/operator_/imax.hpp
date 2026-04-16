#ifndef PYTHONIC_INCLUDE_OPERATOR_IMAX_HPP
#define PYTHONIC_INCLUDE_OPERATOR_IMAX_HPP

#include "pythonic/include/numpy/maximum.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{
  template <class A, class B>
  auto imax(A &&a, B &&b)
      -> std::enable_if_t<std::is_const<A>::value || !std::is_assignable<A, B>::value,
                          decltype(numpy::functor::maximum{}(std::forward<A>(a),
                                                             std::forward<B>(b)))>;

  template <class A, class B>
  auto imax(A &&a, B &&b)
      -> std::enable_if_t<!std::is_const<A>::value && std::is_assignable<A, B>::value,
                          decltype(a = numpy::functor::maximum{}(std::forward<A>(a),
                                                                 std::forward<B>(b)))>;

  DEFINE_FUNCTOR(pythonic::operator_, imax);
} // namespace operator_
PYTHONIC_NS_END

#endif
