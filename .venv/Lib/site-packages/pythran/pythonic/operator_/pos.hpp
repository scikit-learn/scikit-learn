#ifndef PYTHONIC_OPERATOR_POS_HPP
#define PYTHONIC_OPERATOR_POS_HPP

#include "pythonic/include/operator_/pos.hpp"

#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{

  template <class A>
  A pos(A const &a)
  {
    return a;
  }
} // namespace operator_
PYTHONIC_NS_END

#endif
