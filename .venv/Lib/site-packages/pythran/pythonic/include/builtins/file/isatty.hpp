#ifndef PYTHONIC_INCLUDE_BUILTIN_FILE_ISATTY_HPP
#define PYTHONIC_INCLUDE_BUILTIN_FILE_ISATTY_HPP

#include "pythonic/include/types/file.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace file
  {

    bool isatty(types::file const &f);

    DEFINE_FUNCTOR(pythonic::builtins::file, isatty);
  } // namespace file
} // namespace builtins
PYTHONIC_NS_END
#endif
