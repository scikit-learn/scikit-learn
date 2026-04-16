#ifndef PYTHONIC_BUILTIN_FILE_FILENO_HPP
#define PYTHONIC_BUILTIN_FILE_FILENO_HPP

#include "pythonic/include/builtins/file/fileno.hpp"

#include "pythonic/types/file.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace file
  {

    long fileno(types::file const &f)
    {
      return f.fileno();
    }
  } // namespace file
} // namespace builtins
PYTHONIC_NS_END
#endif
