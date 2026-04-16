#ifndef PYTHONIC_INCLUDE_BUILTIN_STR_UPPER_HPP
#define PYTHONIC_INCLUDE_BUILTIN_STR_UPPER_HPP

#include "pythonic/include/types/str.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace str
  {

    types::str upper(types::str const &s);

    DEFINE_FUNCTOR(pythonic::builtins::str, upper);
  } // namespace str
} // namespace builtins
PYTHONIC_NS_END
#endif
