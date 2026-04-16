#ifndef PYTHONIC_BUILTIN_POW_HPP
#define PYTHONIC_BUILTIN_POW_HPP

#include "pythonic/include/builtins/pow.hpp"

#include "pythonic/numpy/power.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{
  inline double pow(long x, long y)
  {
    return std::pow((double)x, (double)y);
  }

  template <long N>
  long pow(long x, std::integral_constant<long, N>)
  {
    if (N == 0)
      return 1;
    if (N == 1)
      return x;
    long tmp = pow(x, std::integral_constant<long, N / 2>{});
    if (N % 2 == 0)
      return tmp * tmp;
    else
      return tmp * tmp * x;
  }

  template <class... Types>
  auto pow(Types &&...args) -> decltype(numpy::functor::power{}(std::forward<Types>(args)...))
  {
    return numpy::functor::power{}(std::forward<Types>(args)...);
  }
} // namespace builtins
PYTHONIC_NS_END

#endif
