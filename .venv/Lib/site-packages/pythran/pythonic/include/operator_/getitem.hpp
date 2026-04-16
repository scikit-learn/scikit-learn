#ifndef PYTHONIC_INCLUDE_OPERATOR_GETITEM_HPP
#define PYTHONIC_INCLUDE_OPERATOR_GETITEM_HPP

#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{
  template <class A, class B>
  auto getitem(A &&a, B &&b) -> decltype(std::forward<A>(a)[std::forward<B>(b)]);

  DEFINE_FUNCTOR(pythonic::operator_, getitem);
} // namespace operator_
PYTHONIC_NS_END

#endif
