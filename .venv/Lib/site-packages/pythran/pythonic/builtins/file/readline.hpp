#ifndef PYTHONIC_BUILTIN_FILE_READLINE_HPP
#define PYTHONIC_BUILTIN_FILE_READLINE_HPP

#include "pythonic/include/builtins/file/readline.hpp"

#include "pythonic/types/file.hpp"
#include "pythonic/types/str.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace file
  {

    inline types::str readline(types::file &f, long size)
    {
      return size < 0 ? f.readline() : f.readline(size);
    }

    inline types::str readline(types::file &&f, long size)
    {
      return size < 0 ? f.readline() : f.readline(size);
    }
  } // namespace file
} // namespace builtins
PYTHONIC_NS_END
#endif
