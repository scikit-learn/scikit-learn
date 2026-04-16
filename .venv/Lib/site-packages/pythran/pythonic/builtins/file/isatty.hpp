#ifndef PYTHONIC_BUILTIN_FILE_ISATTY_HPP
#define PYTHONIC_BUILTIN_FILE_ISATTY_HPP

#include "pythonic/include/builtins/file/isatty.hpp"

#include "pythonic/types/file.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace file
  {

    bool isatty(types::file const &f)
    {
      return f.isatty();
    }
  } // namespace file
} // namespace builtins
PYTHONIC_NS_END
#endif
