#ifndef PYTHONIC_INCLUDE_BUILTIN_PYTHRAN_STATICIFRETURN_HPP
#define PYTHONIC_INCLUDE_BUILTIN_PYTHRAN_STATICIFRETURN_HPP

#include "pythonic/include/types/static_if.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace pythran
  {
    template <class T>
    types::StaticIfReturn<T> StaticIfReturn(T const &arg);

    DEFINE_FUNCTOR(pythonic::builtins::pythran, StaticIfReturn);
  } // namespace pythran
} // namespace builtins
PYTHONIC_NS_END

#endif
