#ifndef PYTHONIC_INCLUDE_BUILTIN_TYPE_HPP
#define PYTHONIC_INCLUDE_BUILTIN_TYPE_HPP

#include "pythonic/include/types/type.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{
  template <class T>
  typename types::type_functor<T>::type type(T const &t);
  DEFINE_FUNCTOR(pythonic::builtins, type);
} // namespace builtins
PYTHONIC_NS_END

#endif
