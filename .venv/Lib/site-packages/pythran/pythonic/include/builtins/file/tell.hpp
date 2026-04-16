#ifndef PYTHONIC_INCLUDE_BUILTIN_FILE_TELL_HPP
#define PYTHONIC_INCLUDE_BUILTIN_FILE_TELL_HPP

#include "pythonic/include/types/file.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace file
  {

    long tell(types::file const &f);

    DEFINE_FUNCTOR(pythonic::builtins::file, tell);
  } // namespace file
} // namespace builtins
PYTHONIC_NS_END
#endif
