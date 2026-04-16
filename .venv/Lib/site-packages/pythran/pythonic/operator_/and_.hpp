#ifndef PYTHONIC_OPERATOR_AND_HPP
#define PYTHONIC_OPERATOR_AND_HPP

#include "pythonic/include/operator_/and_.hpp"
#include "pythonic/operator_/overloads.hpp"

#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{
  template <class A, class B>
  auto and_(A &&a, B &&b) -> decltype(std::forward<A>(a) & std::forward<B>(b))
  {
    return std::forward<A>(a) & std::forward<B>(b);
  }

  DEFINE_ALL_OPERATOR_OVERLOADS_IMPL(and_, &, true)
} // namespace operator_
PYTHONIC_NS_END

#endif
