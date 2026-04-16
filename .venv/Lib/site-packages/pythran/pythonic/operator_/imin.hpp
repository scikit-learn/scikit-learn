#ifndef PYTHONIC_OPERATOR_IMIN_HPP
#define PYTHONIC_OPERATOR_IMIN_HPP

#include "pythonic/include/operator_/imin.hpp"

#include "pythonic/numpy/minimum.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{

  template <class A, class B>
  auto imin(A &&a, B &&b)
      -> std::enable_if_t<std::is_const<A>::value || !std::is_assignable<A, B>::value,
                          decltype(numpy::functor::minimum{}(std::forward<A>(a),
                                                             std::forward<B>(b)))>
  {
    return numpy::functor::minimum{}(std::forward<A>(a), std::forward<B>(b));
  }

  template <class A, class B>
  auto imin(A &&a, B &&b)
      -> std::enable_if_t<!std::is_const<A>::value && std::is_assignable<A, B>::value,
                          decltype(a = numpy::functor::minimum{}(std::forward<A>(a),
                                                                 std::forward<B>(b)))>
  {
    return a = numpy::functor::minimum{}(a, std::forward<B>(b));
  }
} // namespace operator_
PYTHONIC_NS_END

#endif
