#ifndef PYTHONIC_BUILTIN_PYTHRAN_STATIC_IF_HPP
#define PYTHONIC_BUILTIN_PYTHRAN_STATIC_IF_HPP

#include "pythonic/builtins/pythran/is_none.hpp"
#include "pythonic/include/builtins/pythran/static_if.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace pythran
  {

    template <class T, class F0, class F1>
    auto static_if(T const &cond, F0 f0, F1 f1) -> decltype(details::static_if<T>{cond}(f0, f1))
    {
      return details::static_if<T>{cond}(f0, f1);
    }
  } // namespace pythran
} // namespace builtins
PYTHONIC_NS_END

#endif
