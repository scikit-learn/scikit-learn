#ifndef PYTHONIC_BUILTIN_FILE_CLOSE_HPP
#define PYTHONIC_BUILTIN_FILE_CLOSE_HPP

#include "pythonic/include/builtins/file/close.hpp"

#include "pythonic/types/file.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace file
  {

    inline void close(types::file &f)
    {
      f.close();
    }

    inline void close(types::file &&f)
    {
      f.close();
    }
  } // namespace file
} // namespace builtins
PYTHONIC_NS_END
#endif
