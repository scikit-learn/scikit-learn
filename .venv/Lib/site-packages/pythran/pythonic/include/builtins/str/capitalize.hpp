#ifndef PYTHONIC_INCLUDE_BUILTIN_STR_CAPITALIZE_HPP
#define PYTHONIC_INCLUDE_BUILTIN_STR_CAPITALIZE_HPP

#include "pythonic/include/types/str.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace str
  {

    types::str capitalize(types::str const &s);

    DEFINE_FUNCTOR(pythonic::builtins::str, capitalize);
  } // namespace str
} // namespace builtins
PYTHONIC_NS_END
#endif
