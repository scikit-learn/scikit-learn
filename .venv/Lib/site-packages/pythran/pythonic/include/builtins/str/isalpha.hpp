#ifndef PYTHONIC_INCLUDE_BUILTIN_STR_ISALPHA_HPP
#define PYTHONIC_INCLUDE_BUILTIN_STR_ISALPHA_HPP

#include "pythonic/include/types/str.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace str
  {

    bool isalpha(types::str const &s);

    DEFINE_FUNCTOR(pythonic::builtins::str, isalpha);
  } // namespace str
} // namespace builtins
PYTHONIC_NS_END
#endif
