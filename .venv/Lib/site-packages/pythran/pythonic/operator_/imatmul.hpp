#ifndef PYTHONIC_OPERATOR_IMATMUL_HPP
#define PYTHONIC_OPERATOR_IMATMUL_HPP

#include "pythonic/include/operator_/imatmul.hpp"

#include "pythonic/numpy/dot.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{
  template <class A, class B>
  A imatmul(A const &a, B &&b)
  {
    return numpy::functor::dot{}(a, std::forward<B>(b));
  }

  template <class A, class B>
  A &imatmul(A &a, B &&b)
  {
    return a = numpy::functor::dot{}(a,
                                     std::forward<B>(b)); // FIXME: improve that
  }
} // namespace operator_
PYTHONIC_NS_END

#endif
