#ifndef PYTHONIC_BUILTIN_PYTHRAN_STATICIFRETURN_HPP
#define PYTHONIC_BUILTIN_PYTHRAN_STATICIFRETURN_HPP

#include "pythonic/include/builtins/pythran/StaticIfReturn.hpp"
#include "pythonic/types/static_if.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{
  namespace pythran
  {

    template <class T>
    types::StaticIfReturn<T> StaticIfReturn(T const &arg)
    {
      return {arg};
    }
  } // namespace pythran
} // namespace builtins
PYTHONIC_NS_END

#endif
