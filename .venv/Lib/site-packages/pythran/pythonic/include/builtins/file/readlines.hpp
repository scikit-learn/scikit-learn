#ifndef PYTHONIC_INCLUDE_BUILTIN_FILE_READLINES_HPP
#define PYTHONIC_INCLUDE_BUILTIN_FILE_READLINES_HPP

#include "pythonic/include/types/file.hpp"
#include "pythonic/include/types/list.hpp"
#include "pythonic/include/types/str.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace file
  {

    template <class F>
    types::list<types::str> readlines(F &&f);
    template <class F>
    types::list<types::str> readlines(F &&f, long sizehint);

    DEFINE_FUNCTOR(pythonic::builtins::file, readlines);
  } // namespace file
} // namespace builtins
PYTHONIC_NS_END
#endif
