#ifndef PYTHONIC_INCLUDE_BUILTIN_BIN_HPP
#define PYTHONIC_INCLUDE_BUILTIN_BIN_HPP

#include "pythonic/include/utils/functor.hpp"

#include "pythonic/include/types/str.hpp"

#include <type_traits>

PYTHONIC_NS_BEGIN

namespace builtins
{
  template <class T>
  std::enable_if_t<std::is_scalar<T>::value, types::str> bin(T const &v);

  DEFINE_FUNCTOR(pythonic::builtins, bin);
} // namespace builtins
PYTHONIC_NS_END

#endif
