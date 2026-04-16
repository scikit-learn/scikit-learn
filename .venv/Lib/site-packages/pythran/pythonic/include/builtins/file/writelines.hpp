#ifndef PYTHONIC_INCLUDE_BUILTIN_FILE_WRITELINES_HPP
#define PYTHONIC_INCLUDE_BUILTIN_FILE_WRITELINES_HPP

#include "pythonic/include/types/file.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace file
  {

    template <class F, class T>
    void writelines(F &&f, T const &sequence);

    DEFINE_FUNCTOR(pythonic::builtins::file, writelines);
  } // namespace file
} // namespace builtins
PYTHONIC_NS_END
#endif
