#ifndef PYTHONIC_INCLUDE_OPERATOR_POS_HPP
#define PYTHONIC_INCLUDE_OPERATOR_POS_HPP

#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{

  template <class A>
  A pos(A const &a);

  DEFINE_FUNCTOR(pythonic::operator_, pos);
} // namespace operator_
PYTHONIC_NS_END

#endif
