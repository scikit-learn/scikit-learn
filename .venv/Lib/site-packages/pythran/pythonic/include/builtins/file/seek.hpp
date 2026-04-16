#ifndef PYTHONIC_INCLUDE_BUILTIN_FILE_SEEK_HPP
#define PYTHONIC_INCLUDE_BUILTIN_FILE_SEEK_HPP

#include "pythonic/include/types/file.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace file
  {

    void seek(types::file &f, long offset);
    void seek(types::file &&f, long offset);
    void seek(types::file &f, long offset, long whence);
    void seek(types::file &&f, long offset, long whence);

    DEFINE_FUNCTOR(pythonic::builtins::file, seek);
  } // namespace file
} // namespace builtins
PYTHONIC_NS_END
#endif
