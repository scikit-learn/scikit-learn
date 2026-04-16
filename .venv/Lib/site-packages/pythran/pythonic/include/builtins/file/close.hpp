#ifndef PYTHONIC_INCLUDE_BUILTIN_FILE_CLOSE_HPP
#define PYTHONIC_INCLUDE_BUILTIN_FILE_CLOSE_HPP

#include "pythonic/include/types/file.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace file
  {

    void close(types::file &f);
    void close(types::file &&f);

    DEFINE_FUNCTOR(pythonic::builtins::file, close);
  } // namespace file
} // namespace builtins
PYTHONIC_NS_END
#endif
