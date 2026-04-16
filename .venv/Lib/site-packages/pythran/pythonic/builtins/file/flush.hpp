#ifndef PYTHONIC_BUILTIN_FILE_FLUSH_HPP
#define PYTHONIC_BUILTIN_FILE_FLUSH_HPP

#include "pythonic/include/builtins/file/flush.hpp"

#include "pythonic/types/file.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace file
  {

    inline void flush(types::file &f)
    {
      f.flush();
    }

    inline void flush(types::file &&f)
    {
      f.flush();
    }
  } // namespace file
} // namespace builtins
PYTHONIC_NS_END
#endif
