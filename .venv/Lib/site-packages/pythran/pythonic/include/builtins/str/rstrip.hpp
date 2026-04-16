#ifndef PYTHONIC_INCLUDE_BUILTIN_STR_RSTRIP_HPP
#define PYTHONIC_INCLUDE_BUILTIN_STR_RSTRIP_HPP

#include "pythonic/include/types/str.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace str
  {

    types::str rstrip(types::str const &self, types::str const &to_del = " ");

    DEFINE_FUNCTOR(pythonic::builtins::str, rstrip);
  } // namespace str
} // namespace builtins
PYTHONIC_NS_END
#endif
