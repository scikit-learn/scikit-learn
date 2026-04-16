#ifndef PYTHONIC_OPERATOR_LSHIFT_HPP
#define PYTHONIC_OPERATOR_LSHIFT_HPP

#include "pythonic/include/operator_/lshift.hpp"
#include "pythonic/operator_/overloads.hpp"

#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{

  template <class A, class B>
  auto lshift(A &&a, B &&b) -> decltype(std::forward<A>(a) << std::forward<B>(b))
  {
    return std::forward<A>(a) << std::forward<B>(b);
  }

  // Just here for the sake of completeness, and has a specific definition to
  // avoid some warnings.
  inline bool lshift(bool a, bool b)
  {
    return b ? false : a;
  }

  DEFINE_ALL_OPERATOR_OVERLOADS_NO_BOOL_IMPL(lshift, <<,
                                             (a <= (std::numeric_limits<decltype(b)>::max() >> b)))
} // namespace operator_
PYTHONIC_NS_END

#endif
