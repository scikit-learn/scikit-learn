#ifndef PYTHONIC_OPERATOR_ISNOT_HPP
#define PYTHONIC_OPERATOR_ISNOT_HPP

#include "pythonic/include/operator_/is_not.hpp"

#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{

  template <class A, class B>
  auto is_not(A &&a, B &&b)
      -> decltype(builtins::id(std::forward<A>(a)) != builtins::id(std::forward<B>(b)))
  {
    return builtins::id(std::forward<A>(a)) != builtins::id(std::forward<B>(b));
  }
} // namespace operator_
PYTHONIC_NS_END

#endif
