#ifndef PYTHONIC_BUILTIN_ROUND_HPP
#define PYTHONIC_BUILTIN_ROUND_HPP

#include "pythonic/include/builtins/round.hpp"

#include "pythonic/builtins/pow.hpp"
#include "pythonic/utils/functor.hpp"

#include <cmath>

PYTHONIC_NS_BEGIN

namespace builtins
{

  template <class T>
  double round(T const &v, size_t n)
  {
    T p = functor::pow()(10, n);
    return std::lround(v * p) / p;
  }

  template <class T>
  double round(T const &v)
  {
    return std::lround(v);
  }
} // namespace builtins
PYTHONIC_NS_END

#endif
