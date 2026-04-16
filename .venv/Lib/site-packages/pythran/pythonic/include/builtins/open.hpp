#ifndef PYTHONIC_INCLUDE_BUILTIN_OPEN_HPP
#define PYTHONIC_INCLUDE_BUILTIN_OPEN_HPP

#include "pythonic/include/types/file.hpp"
#include "pythonic/include/types/str.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  types::file open(types::str const &filename, types::str const &strmode = "r");

  DEFINE_FUNCTOR(pythonic::builtins, open);
} // namespace builtins
PYTHONIC_NS_END

#endif
