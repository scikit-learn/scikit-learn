#ifndef PYTHONIC_BUILTIN_PYTHRAN_AND_HPP
#define PYTHONIC_BUILTIN_PYTHRAN_AND_HPP

#include "pythonic/include/builtins/pythran/and_.hpp"

#include "pythonic/types/combined.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace pythran
  {

    template <class T0, class T1>
    types::lazy_combined_t<T0, T1> and_(T0 &&v0, T1 &&v1)
    {
      auto &&val0 = std::forward<T0>(v0)();
      if (val0)
        return (types::lazy_combined_t<T0, T1>)std::forward<T1>(v1)();
      else
        return (types::lazy_combined_t<T0, T1>)val0;
    }
  } // namespace pythran
} // namespace builtins
PYTHONIC_NS_END

#endif
