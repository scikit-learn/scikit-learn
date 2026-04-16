#ifndef PYTHONIC_INCLUDE_OPERATOR_IMOD_HPP
#define PYTHONIC_INCLUDE_OPERATOR_IMOD_HPP

#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{

  template <class A, class B>
  A &imod(A &a, B &&b);
  template <class A, class B>
  A imod(A const &a, B &&b);

  DEFINE_FUNCTOR(pythonic::operator_, imod);
} // namespace operator_
PYTHONIC_NS_END

#endif
