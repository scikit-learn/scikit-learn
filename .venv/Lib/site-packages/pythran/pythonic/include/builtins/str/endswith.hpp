#ifndef PYTHONIC_INCLUDE_BUILTIN_STR_ENDSWITH_HPP
#define PYTHONIC_INCLUDE_BUILTIN_STR_ENDSWITH_HPP

#include "pythonic/include/types/str.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace str
  {

    // TODO : Add implementation for tuple as first argument.
    bool endswith(types::str const &s, types::str const &suffix, long start = 0, long end = -1);

    DEFINE_FUNCTOR(pythonic::builtins::str, endswith);
  } // namespace str
} // namespace builtins
PYTHONIC_NS_END
#endif
