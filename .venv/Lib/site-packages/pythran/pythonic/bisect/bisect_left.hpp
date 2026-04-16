#ifndef PYTHONIC_BISECT_BISECTLEFT_HPP
#define PYTHONIC_BISECT_BISECTLEFT_HPP

#include "pythonic/include/bisect/bisect_left.hpp"

#include "pythonic/bisect/bisect.hpp"
#include "pythonic/utils/functor.hpp"

#include <algorithm>

PYTHONIC_NS_BEGIN

namespace bisect
{
  template <class X, class A>
  long bisect_left(X const &x, A const &a, long lo)
  {
    return bisect(x, a, lo, std::lower_bound<typename X::const_iterator, A>);
  }

  template <class X, class A>
  long bisect_left(X const &x, A const &a, long lo, long hi)
  {
    return bisect(x, a, lo, hi, std::lower_bound<typename X::const_iterator, A>);
  }
} // namespace bisect
PYTHONIC_NS_END

#endif
