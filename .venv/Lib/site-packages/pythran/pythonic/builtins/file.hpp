#ifndef PYTHONIC_BUILTIN_FILE_HPP
#define PYTHONIC_BUILTIN_FILE_HPP

#include "pythonic/include/builtins/file.hpp"

#include "pythonic/types/file.hpp"
#include "pythonic/types/str.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace anonymous
  {
    types::file file(types::str const &filename, types::str const &strmode)
    {
      return {filename, strmode};
    }
  } // namespace anonymous
} // namespace builtins
PYTHONIC_NS_END

#endif
