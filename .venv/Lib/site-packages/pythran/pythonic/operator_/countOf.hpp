#ifndef PYTHONIC_OPERATOR_COUNTOF_HPP
#define PYTHONIC_OPERATOR_COUNTOF_HPP

#include "pythonic/include/operator_/countOf.hpp"

#include "pythonic/utils/functor.hpp"
#include <algorithm>

PYTHONIC_NS_BEGIN

namespace operator_
{
  template <class A, class B>
  long countOf(A &&a, B &&b)
  {
    return std::count(a.begin(), a.end(), std::forward<B>(b));
  }
} // namespace operator_
PYTHONIC_NS_END

#endif
