#ifndef PYTHONIC_INCLUDE_OPERATOR_COUNTOF_HPP
#define PYTHONIC_INCLUDE_OPERATOR_COUNTOF_HPP

#include "pythonic/include/utils/functor.hpp"
#include <algorithm>

PYTHONIC_NS_BEGIN

namespace operator_
{
  template <class A, class B>
  long countOf(A &&a, B &&b);

  DEFINE_FUNCTOR(pythonic::operator_, countOf);
} // namespace operator_
PYTHONIC_NS_END

#endif
