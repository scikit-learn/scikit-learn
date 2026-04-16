#ifndef PYTHONIC_INCLUDE_OPERATOR_IS_HPP
#define PYTHONIC_INCLUDE_OPERATOR_IS_HPP

#include "pythonic/include/builtins/id.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{

  template <class A, class B>
  auto is_(A &&a, B &&b)
      -> decltype(builtins::id(std::forward<A>(a)) == builtins::id(std::forward<B>(b)));

  DEFINE_FUNCTOR(pythonic::operator_, is_);
} // namespace operator_
PYTHONIC_NS_END

#endif
