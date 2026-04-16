#ifndef PYTHONIC_INCLUDE_BUILTIN_FILE_TRUNCATE_HPP
#define PYTHONIC_INCLUDE_BUILTIN_FILE_TRUNCATE_HPP

#include "pythonic/include/types/file.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace file
  {

    void truncate(types::file &f);
    void truncate(types::file &&f);
    void truncate(types::file &f, long size);
    void truncate(types::file &&f, long size);

    DEFINE_FUNCTOR(pythonic::builtins::file, truncate);
  } // namespace file
} // namespace builtins
PYTHONIC_NS_END
#endif
