#ifndef PYTHONIC_BUILTIN_TYPE_HPP
#define PYTHONIC_BUILTIN_TYPE_HPP

#include "pythonic/include/builtins/type.hpp"
#include "pythonic/types/type.hpp"

#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{
  template <class T>
  typename types::type_functor<T>::type type(T const &)
  {
    return {};
  }
} // namespace builtins
PYTHONIC_NS_END

#endif
