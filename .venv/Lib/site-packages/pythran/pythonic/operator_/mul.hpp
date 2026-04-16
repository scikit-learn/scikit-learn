#ifndef PYTHONIC_OPERATOR_MUL_HPP
#define PYTHONIC_OPERATOR_MUL_HPP

#include "pythonic/include/operator_/mul.hpp"

#include "pythonic/operator_/overloads.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{

  template <class A, class B>
  auto mul(A &&a, B &&b) -> decltype(std::forward<A>(a) * std::forward<B>(b))
  {
    return std::forward<A>(a) * std::forward<B>(b);
  }

  // Just here for the sake of completeness, and has a specific definition to
  // avoid some warnings.
  inline bool mul(bool a, bool b)
  {
    return a & b;
  }

  DEFINE_ALL_OPERATOR_OVERLOADS_NO_BOOL_IMPL(
      mul, *,
      (b == 0 ||
       (a * b >= 0 && std::abs(a) <= std::numeric_limits<decltype(b)>::max() / std::abs(b)) ||
       (a * b <= 0 && std::abs(a) >= std::numeric_limits<decltype(b)>::min() / std::abs(b))))
} // namespace operator_
PYTHONIC_NS_END

#endif
