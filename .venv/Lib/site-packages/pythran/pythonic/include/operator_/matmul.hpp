#ifndef PYTHONIC_INCLUDE_OPERATOR_MATMUL_HPP
#define PYTHONIC_INCLUDE_OPERATOR_MATMUL_HPP

#include "pythonic/include/numpy/dot.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{

  template <class A, class B>
  auto matmul(A &&a, B &&b)
      -> decltype(numpy::functor::dot{}(std::forward<A>(a), std::forward<B>(b)));

  DEFINE_FUNCTOR(pythonic::operator_, matmul);
} // namespace operator_
PYTHONIC_NS_END

#endif
