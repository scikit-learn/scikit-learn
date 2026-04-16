#ifndef PYTHONIC_INCLUDE_BUILTIN_FILE_WRITE_HPP
#define PYTHONIC_INCLUDE_BUILTIN_FILE_WRITE_HPP

#include "pythonic/include/types/file.hpp"
#include "pythonic/include/types/str.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace file
  {

    long write(types::file &f, types::str const &str);
    long write(types::file &&f, types::str const &str);

    DEFINE_FUNCTOR(pythonic::builtins::file, write);
  } // namespace file
} // namespace builtins
PYTHONIC_NS_END
#endif
