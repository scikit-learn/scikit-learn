#ifndef PYTHONIC_BUILTIN_FILE_SEEK_HPP
#define PYTHONIC_BUILTIN_FILE_SEEK_HPP

#include "pythonic/include/builtins/file/seek.hpp"

#include "pythonic/types/file.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace file
  {

    inline void seek(types::file &f, long offset)
    {
      f.seek(offset);
    }

    inline void seek(types::file &&f, long offset)
    {
      // Nothing have to be done as it is a lvalue
    }

    inline void seek(types::file &f, long offset, long whence)
    {
      f.seek(offset, whence);
    }

    inline void seek(types::file &&f, long offset, long whence)
    {
      // Nothing have to be done as it is a lvalue
    }
  } // namespace file
} // namespace builtins
PYTHONIC_NS_END
#endif
