#ifndef PYTHONIC_OPERATOR_IMAX_HPP
#define PYTHONIC_OPERATOR_IMAX_HPP

#include "pythonic/include/operator_/imax.hpp"

#include "pythonic/numpy/maximum.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{
  template <class A, class B>
  auto imax(A &&a, B &&b)
      -> std::enable_if_t<std::is_const<A>::value || !std::is_assignable<A, B>::value,
                          decltype(numpy::functor::maximum{}(std::forward<A>(a),
                                                             std::forward<B>(b)))>
  {
    return numpy::functor::maximum{}(std::forward<A>(a), std::forward<B>(b));
  }

  template <class A, class B>
  auto imax(A &&a, B &&b)
      -> std::enable_if_t<!std::is_const<A>::value && std::is_assignable<A, B>::value,
                          decltype(a = numpy::functor::maximum{}(std::forward<A>(a),
                                                                 std::forward<B>(b)))>
  {
    return a = numpy::functor::maximum{}(a, std::forward<B>(b));
  }
} // namespace operator_
PYTHONIC_NS_END

#endif
