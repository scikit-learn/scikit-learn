#ifndef PYTHONIC_BUILTIN_FILE_TELL_HPP
#define PYTHONIC_BUILTIN_FILE_TELL_HPP

#include "pythonic/include/builtins/file/tell.hpp"

#include "pythonic/types/file.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace file
  {

    long tell(types::file const &f)
    {
      return f.tell();
    }
  } // namespace file
} // namespace builtins
PYTHONIC_NS_END
#endif
