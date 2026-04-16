#ifndef PYTHONIC_INCLUDE_BUILTIN_ROUND_HPP
#define PYTHONIC_INCLUDE_BUILTIN_ROUND_HPP

#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{
  template <class T>
  double round(T const &v, size_t n);

  template <class T>
  double round(T const &v);

  DEFINE_FUNCTOR(pythonic::builtins, round);
} // namespace builtins
PYTHONIC_NS_END

#endif
