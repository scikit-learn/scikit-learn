#ifndef PYTHONIC_INCLUDE_BUILTIN_FILE_READLINE_HPP
#define PYTHONIC_INCLUDE_BUILTIN_FILE_READLINE_HPP

#include "pythonic/include/types/file.hpp"
#include "pythonic/include/types/str.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace file
  {

    types::str readline(types::file &f, long size = -1);
    types::str readline(types::file &&f, long size = -1);

    DEFINE_FUNCTOR(pythonic::builtins::file, readline);
  } // namespace file
} // namespace builtins
PYTHONIC_NS_END
#endif
