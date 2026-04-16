#ifndef PYTHONIC_BUILTIN_PYTHRAN_STATICIFCONT_HPP
#define PYTHONIC_BUILTIN_PYTHRAN_STATICIFCONT_HPP

#include "pythonic/include/builtins/pythran/StaticIfCont.hpp"
#include "pythonic/types/static_if.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace pythran
  {
    template <class T>
    types::StaticIfCont<T> StaticIfCont(T const &arg)
    {
      return {arg};
    }
  } // namespace pythran
} // namespace builtins
PYTHONIC_NS_END

#endif
