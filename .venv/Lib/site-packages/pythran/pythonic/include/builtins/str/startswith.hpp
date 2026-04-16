#ifndef PYTHONIC_INCLUDE_BUILTIN_STR_STARTSWITH_HPP
#define PYTHONIC_INCLUDE_BUILTIN_STR_STARTSWITH_HPP

#include "pythonic/include/types/str.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace str
  {

    bool startswith(types::str const &s, types::str const &prefix, long start = 0, long end = -1);

    DEFINE_FUNCTOR(pythonic::builtins::str, startswith);
  } // namespace str
} // namespace builtins
PYTHONIC_NS_END
#endif
