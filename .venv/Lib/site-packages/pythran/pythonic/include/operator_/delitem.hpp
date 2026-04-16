#ifndef PYTHONIC_INCLUDE_OPERATOR_DELITEM_HPP
#define PYTHONIC_INCLUDE_OPERATOR_DELITEM_HPP

#include "pythonic/include/builtins/None.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{
  template <class A, class B>
  types::none_type delitem(A &&a, B &&b);

  DEFINE_FUNCTOR(pythonic::operator_, delitem);
} // namespace operator_
PYTHONIC_NS_END

#endif
