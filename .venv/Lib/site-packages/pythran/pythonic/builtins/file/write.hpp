#ifndef PYTHONIC_BUILTIN_FILE_WRITE_HPP
#define PYTHONIC_BUILTIN_FILE_WRITE_HPP

#include "pythonic/include/builtins/file/write.hpp"

#include "pythonic/types/file.hpp"
#include "pythonic/types/str.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace file
  {

    inline long write(types::file &f, types::str const &str)
    {
      return f.write(str);
    }

    inline long write(types::file &&f, types::str const &str)
    {
      return f.write(str);
    }
  } // namespace file
} // namespace builtins
PYTHONIC_NS_END
#endif
