#ifndef PYTHONIC_INCLUDE_BUILTIN_PYTHRAN_STATICIFCONT_HPP
#define PYTHONIC_INCLUDE_BUILTIN_PYTHRAN_STATICIFCONT_HPP

#include "pythonic/include/types/static_if.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace pythran
  {
    template <class T>
    types::StaticIfCont<T> StaticIfCont(T const &arg);

    DEFINE_FUNCTOR(pythonic::builtins::pythran, StaticIfCont);
  } // namespace pythran
} // namespace builtins

PYTHONIC_NS_END

#endif
