#ifndef PYTHONIC_INCLUDE_BUILTIN_POW_HPP
#define PYTHONIC_INCLUDE_BUILTIN_POW_HPP

#include "pythonic/include/numpy/power.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  // this is only the case in python if the exponent is negative
  double pow(long, long);
  // in that case we are sure we have a positive exponent
  template <long N>
  long pow(long, std::integral_constant<long, N>);

  template <class... Types>
  auto pow(Types &&...args) -> decltype(numpy::functor::power{}(std::forward<Types>(args)...));

  DEFINE_FUNCTOR(pythonic::builtins, pow);
} // namespace builtins
PYTHONIC_NS_END

#endif
