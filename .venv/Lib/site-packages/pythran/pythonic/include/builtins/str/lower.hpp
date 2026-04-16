#ifndef PYTHONIC_INCLUDE_BUILTIN_STR_LOWER_HPP
#define PYTHONIC_INCLUDE_BUILTIN_STR_LOWER_HPP

#include "pythonic/include/types/str.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace str
  {

    types::str lower(types::str const &s);

    DEFINE_FUNCTOR(pythonic::builtins::str, lower);
  } // namespace str
} // namespace builtins
PYTHONIC_NS_END
#endif
