#ifndef PYTHONIC_INCLUDE_BUILTIN_STR_STRIP_HPP
#define PYTHONIC_INCLUDE_BUILTIN_STR_STRIP_HPP

#include "pythonic/include/types/str.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace str
  {

    types::str strip(types::str const &self, types::str const &to_del = " \n");

    DEFINE_FUNCTOR(pythonic::builtins::str, strip);
  } // namespace str
} // namespace builtins
PYTHONIC_NS_END
#endif
