#ifndef PYTHONIC_BUILTIN_OPEN_HPP
#define PYTHONIC_BUILTIN_OPEN_HPP

#include "pythonic/include/builtins/open.hpp"

#include "pythonic/types/file.hpp"
#include "pythonic/types/str.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  types::file open(types::str const &filename, types::str const &strmode)
  {
    return {filename, strmode};
  }
} // namespace builtins
PYTHONIC_NS_END

#endif
