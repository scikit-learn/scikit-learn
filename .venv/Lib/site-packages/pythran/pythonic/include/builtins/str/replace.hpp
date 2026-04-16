#ifndef PYTHONIC_INCLUDE_BUILTIN_STR_REPLACE_HPP
#define PYTHONIC_INCLUDE_BUILTIN_STR_REPLACE_HPP

#include "pythonic/include/types/str.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace str
  {

    types::str replace(types::str const &self, types::str const &old_pattern,
                       types::str const &new_pattern,
                       long count = std::numeric_limits<long>::max());

    DEFINE_FUNCTOR(pythonic::builtins::str, replace);
  } // namespace str
} // namespace builtins
PYTHONIC_NS_END
#endif
