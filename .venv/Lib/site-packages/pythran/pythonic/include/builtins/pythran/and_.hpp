#ifndef PYTHONIC_INCLUDE_BUILTIN_PYTHRAN_AND_HPP
#define PYTHONIC_INCLUDE_BUILTIN_PYTHRAN_AND_HPP

#include "pythonic/include/types/combined.hpp"
#include "pythonic/include/types/lazy.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace pythran
  {

    template <class T0, class T1>
    types::lazy_combined_t<T0, T1> and_(T0 &&, T1 &&);

    DEFINE_FUNCTOR(pythonic::builtins::pythran, and_);
  } // namespace pythran
} // namespace builtins
PYTHONIC_NS_END

#endif
