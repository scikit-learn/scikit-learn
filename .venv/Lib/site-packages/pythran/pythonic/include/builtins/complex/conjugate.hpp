#ifndef PYTHONIC_INCLUDE_BUILTIN_COMPLEX_CONJUGATE_HPP
#define PYTHONIC_INCLUDE_BUILTIN_COMPLEX_CONJUGATE_HPP

#include "pythonic/include/numpy/conjugate.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN
namespace builtins
{
  namespace complex
  {
    USING_FUNCTOR(conjugate, numpy::functor::conjugate);
  }
} // namespace builtins
PYTHONIC_NS_END
#endif
