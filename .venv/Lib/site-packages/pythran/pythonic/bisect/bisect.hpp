#ifndef PYTHONIC_BISECT_BISECT_HPP
#define PYTHONIC_BISECT_BISECT_HPP

#include "pythonic/builtins/ValueError.hpp"
#include "pythonic/include/bisect/bisect.hpp"

#include "pythonic/utils/functor.hpp"

#include <iterator>

PYTHONIC_NS_BEGIN

namespace bisect
{
  template <class X, class A>
  long bisect(X const &x, A const &a, long lo, details::bisect_fun<X, A> const &fun)
  {
    if (lo < 0)
      throw types::ValueError("lo must be non-negative");
    return std::distance(x.begin(), fun(x.begin() + lo, x.end(), a));
  }

  template <class X, class A>
  long bisect(X const &x, A const &a, long lo, long hi, details::bisect_fun<X, A> const &fun)
  {
    if (lo < 0)
      throw types::ValueError("lo must be non-negative");
    return std::distance(x.begin(), fun(x.begin() + lo, x.begin() + hi, a));
  }
} // namespace bisect
PYTHONIC_NS_END

#endif
