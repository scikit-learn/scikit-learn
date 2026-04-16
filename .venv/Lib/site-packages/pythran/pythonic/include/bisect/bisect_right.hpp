#ifndef PYTHONIC_INCLUDE_BISECT_BISECTRIGHT_HPP
#define PYTHONIC_INCLUDE_BISECT_BISECTRIGHT_HPP

#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace bisect
{
  template <class X, class A>
  long bisect_right(X const &x, A const &a, long lo = 0);

  template <class X, class A>
  long bisect_right(X const &x, A const &a, long lo, long hi);

  DEFINE_FUNCTOR(pythonic::bisect, bisect_right);
} // namespace bisect
PYTHONIC_NS_END

#endif
